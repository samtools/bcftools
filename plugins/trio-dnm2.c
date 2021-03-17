/* The MIT License

   Copyright (c) 2018-2021 Genome Research Ltd.

   Author: Petr Danecek <pd3@sanger.ac.uk>
   
   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:
   
   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.
   
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <unistd.h>     // for isatty
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/kfunc.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <assert.h>
#include <errno.h>
#include "bcftools.h"
#include "regidx.h"
#include "filter.h"

#define USE_DNG    1     // DeNovoGear model
#define USE_ACM    2     // the new "allele-centric model" which combines fixed DNG priors with allele centric approach

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define iFATHER 0   // don't modify, QS calculations depend on this order!
#define iMOTHER 1
#define iCHILD  2

typedef struct
{
    int idx[3];     // VCF sample index for child, father, mother
    int pass,       // do all three pass the filters?
        is_male;    // male pattern of chrX inheritance?
}
trio_t;

typedef struct
{
    // combines priors, mutation rates, genotype transmission probability; see init_priors()
    double pprob[10][10][10];           // prior probability; the order is father,mother,child
    uint8_t denovo[10][10][10];         // is the GT combination not compatible with normal inheritence (0) or is de novo (1)
    uint8_t denovo_allele[10][10][10];  // which of the alleles is de novo for this configuration
}
priors_t;

typedef struct
{
    int argc, filter_logic, regions_is_file, targets_is_file, output_type;
    char *filter_str;
    filter_t *filter;
    char **argv, *ped_fname, *pfm, *output_fname, *fname, *regions, *targets;
    htsFile *out_fh;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr, *hdr_out;
    char *chrX_list_str;
    regidx_t *chrX_idx;
    trio_t *trio;
    int has_fmt_ad;
    int ntrio, mtrio;
    int32_t *pl, *ad, *qs, *dnm_qual_int, *dnm_allele, *vaf;    // input FMT/PL, AD, QS values, output DNM and VAF
    float *dnm_qual_float;
    int mpl, mad, mqs;
    double min_score;
    double *aprob;  // proband's allele probabilities
    double *pl3;    // normalized PLs converted to probs for iFATHER,iMOTHER,iCHILD
    double *qs3;    // QS converted to probs for iFATHER,iMOTHER,iCHILD
    int maprob, mpl3, mqs3, midx, *idx, force_ad, use_model;
    char *dnm_score_tag,            // the argument of --use tag, by default DNM:log
         *dnm_vaf_tag,
         *dnm_allele_tag;
    int dnm_score_is_float;         // given by e.g. --use tag DNM:float
    double mrate;                   // --use mrate, mutation rate
    double pnoise_abs,pnoise_frac;  // --use pn|pnoise or --use pns
    int pnoise_strict;              // set to 1 if pns was used or 0 if pn
    int use_ppl, use_ppl_qs;        // --use ppl and --use ppl-qs
    int use_dng_priors;             // --use dng-priors
    priors_t priors, priors_X, priors_XX;
}
args_t;

args_t args;

const char *about(void)
{
    return "Screen variants for possible de-novo mutations in trios.\n";
}

static const char *usage_text(void)
{
    return 
        "\n"
        "About: Screen variants for possible de-novo mutations in trios\n"
        "Usage: bcftools +trio-dnm2 [Plugin Options]\n"
        "Plugin options:\n"
        "   -e, --exclude EXPR              Exclude trios for which the expression is true (one matching sample invalidates a trio)\n"
        "       --force-AD                  Calculate VAF even if the number of FMT/AD fields is incorrect. Use at your own risk!\n"
        "   -i, --include EXPR              Include trios for which the expression is true (one failing samples invalidates a trio)\n"
        "   -m, --min-score NUM             Do not add FMT/DNM annotation if the score is smaller than NUM\n"
        "   -o, --output FILE               Output file name [stdout]\n"
        "   -O, --output-type <b|u|z|v>     b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n"
        "   -p, --pfm [1X:|2X:]P,F,M        Sample names of child (the proband), father, mother; \"1X:\" for male pattern of chrX inheritance [2X:]\n"
        "   -P, --ped FILE                  PED file with the columns: <ignored>,proband,father,mother,sex(1:male,2:female)\n"
        "   -r, --regions REG               Restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE         Restrict to regions listed in a file\n"
        "   -t, --targets REG               Similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE         Similar to -R but streams rather than index-jumps\n"
        "   -u, --use OPTION[=VALUE]        Various options to tweak:\n"
        "          DNG                         Use the original DeNovoGear model, implies -u dng-priors\n"
        "          dng-priors                  Use the original DeNovoGear priors (including bugs in prior assignment)\n"
        "          mrate=NUM                   Mutation rate [-u mrate=1e-8]\n"
        "          pn=FRAC[,NUM]               Tolerance to parental noise or mosaicity, given as fraction of QS or number of reads [-u pn=0,0]\n"
        "          pns=FRAC[,NUM]              Same as `pn` but is not applied to alleles observed in both parents [-u pns=0.045,0]\n"
        "          ppl                         Use parental genotype likelihoods (FMT/PL rather than FMT/QS)\n"
        "          tag=TAG[:phred|log]         Annotation to add, either as phred quality (int) or log-scaled (float) [-u tag=DNM:log]\n"
        "          vaf=TAG                     The tag name for variant allele fraction annotation to add [VAF]\n"
        "          va=TAG                      The tag name for variant allele annotation [VA]\n"
        "   -X, --chrX LIST                 Regions with the chr X inheritance pattern or one of the predefined lists, exclude PARs [GRCh37]\n"
        "                                      GRCh37 .. X:1-60000,chrX:1-60000,X:2699521-154931043,chrX:2699521-154931043\n"
        "                                      GRCh38 .. X:1-9999,chrX:1-9999,X:2781480-155701381,chrX:2781480-155701381\n"
        "\n"
        "Example:\n"
        "   # Annotate VCF with FORMAT/DNM, run for a single trio\n"
        "   bcftools +trio-dnm2 -p proband,father,mother file.bcf\n"
        "\n"
        "   # Same as above, but read the trio(s) from a PED file\n"
        "   bcftools +trio-dnm2 -P file.ped file.bcf\n"
        "\n"
        "   # Same as above plus extract a list of significant DNMs using the bcftools/query command\n"
        "   bcftools +trio-dnm2 -P file.ped file.bcf -Ou | bcftools query -i'DNM>10' -f'[%CHROM:%POS %SAMPLE %DNM\\n]'\n"
        "\n"
        "   # A complete example with a variant calling step. Note that this is one long\n"
        "   # command and should be on a single line. Also note that a filtering step is\n"
        "   # recommended, e.g. by depth and VAF (not shown here):\n"
        "   bcftools mpileup -a AD,QS -f ref.fa -Ou proband.bam father.bam mother.bam |\n"
        "     bcftools call -mv -Ou |\n"
        "       bcftools +trio-dnm2 -p proband,father,mother -Oz -o output.vcf.gz\n"
        "\n";
}

static int cmp_trios(const void *_a, const void *_b)
{
    trio_t *a = (trio_t *) _a;
    trio_t *b = (trio_t *) _b;
    int i;
    int amin = a->idx[0];
    for (i=1; i<3; i++)
        if ( amin > a->idx[i] ) amin = a->idx[i];
    int bmin = b->idx[0];
    for (i=1; i<3; i++)
        if ( bmin > b->idx[i] ) bmin = b->idx[i];
    if ( amin < bmin ) return -1;
    if ( amin > bmin ) return 1;
    return 0;
}
static void parse_ped(args_t *args, char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);

    kstring_t str = {0,0,0};
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Empty file: %s\n", fname);

    int moff = 0, *off = NULL;
    do
    {
        // familyID    sampleID paternalID maternalID sex   phenotype   population relationship   siblings   secondOrder   thirdOrder   children    comment
        // BB03    HG01884 HG01885 HG01956 2   0   ACB child   0   0   0   0
        int ncols = ksplit_core(str.s,0,&moff,&off);
        if ( ncols<4 ) error("Could not parse the ped file: %s\n", str.s);

        int father = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[2]]);
        if ( father<0 ) continue;
        int mother = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[3]]);
        if ( mother<0 ) continue;
        int child = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[1]]);
        if ( child<0 ) continue;

        int sex = 0;
        if ( ncols>=5 )
        {
            char *tmp;
            sex = strtol(&str.s[off[4]],&tmp,10);
            if ( tmp==&str.s[off[4]] || *tmp ) error("Could not parse the PED file, the 5th column should be numeric: %s\n",str.s);
            if ( sex!=1 && sex!=2 ) sex = 0;
        }

        args->ntrio++;
        hts_expand0(trio_t,args->ntrio,args->mtrio,args->trio);
        trio_t *trio = &args->trio[args->ntrio-1];
        trio->idx[iFATHER] = father;
        trio->idx[iMOTHER] = mother;
        trio->idx[iCHILD]  = child;
        trio->is_male = sex==1 ? 1 : 0;
    }
    while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );

    // sort the sample by index so that they are accessed more or less sequentially
    qsort(args->trio,args->ntrio,sizeof(trio_t),cmp_trios);

    // check for duplicates
    int i;
    for (i=1; i<args->ntrio; i++)
    {
        trio_t *ta = &args->trio[i-1];
        trio_t *tb = &args->trio[i];
        if ( ta->idx[0]==tb->idx[0] && ta->idx[1]==tb->idx[1] && ta->idx[2]==tb->idx[2] )
            error("Error: duplicate trio entries detected in the PED file: %s\n",fname);
    }

    fprintf(stderr,"Identified %d complete trio%s in the VCF file\n", args->ntrio,args->ntrio==1?"":"s");

    free(str.s);
    free(off);
    if ( hts_close(fp)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,fname);
}


static const uint8_t seq1[10] = {0,1,1,2,2,2,3,3,3,3};
static const uint8_t seq2[10] = {0,0,1,0,1,2,0,1,2,3};
typedef enum { include_ref, only_alts } count_unique_t;
static int count_unique_alleles(int ngt, int gt[3], count_unique_t count)
{
    int i, als[4] = {0,0,0,0};
    for (i=0; i<ngt; i++)
    {
        int igt = gt[i];
        als[seq1[igt]] = 1;
        als[seq2[igt]] = 1;
    }
    int nals = 0;
    int ibeg = count==include_ref ? 0 : 1;
    for (i=ibeg; i<4; i++) nals += als[i];
    return nals;
}

// Parent genotype probability L(GM,GF)
// The FIGL model from the supplement "Variation in genome-wide mutation rates within and between human families",
// see also the actual implementation in https://github.com/ultimatesource/denovogear/blob/develop/src/dnm/makeLookup.cc
// This is the original method, including bugs in prior assignment.
static double init_DNG_mf_priors(args_t *args, int fi, int mi, int ci)
{
    double gt_prior = 0;    // parent genotype probability L(GM,GF)
    int fa = seq1[fi];
    int fb = seq2[fi];
    int ma = seq1[mi];
    int mb = seq2[mi];
    int gts[3]; gts[0] = fi; gts[1] = mi; gts[2] = 0;
    int nals_mf = count_unique_alleles(2,gts,include_ref);
    int ca = seq1[ci];
    int cb = seq2[ci];
    gts[0] = fi; gts[1] = mi; gts[2] = ci;
    int nals_mfc = count_unique_alleles(3,gts,include_ref);
    int nref_mf  = (fa==0 ? 1 : 0) + (fb==0 ? 1 : 0) + (ma==0 ? 1 : 0) + (mb==0 ? 1 : 0);

    if ( nals_mfc>3 )                                               // 4 different alleles in the trio
        gt_prior = 1e-26;
    else if ( nals_mf>=3 )                                          // 3 different alleles in parents,
        gt_prior = 0.002 * 0.002 / 414;                             //      split equally amongst all triallelic cases
    else if ( nals_mfc==3 )                                         // 3rd allele in the child
        gt_prior = 1e-3 * 1e-3;                                     // This is what g_PolyRate evaluates in DNG code
    else if ( nref_mf==4 )
        gt_prior = 0.995 * 0.998;                                   // 4 copies of ref in parents
    else if ( nref_mf==3 )
        gt_prior = 0.995 * 0.002 * (3.0/5.0) * (4.0/5.0) * 0.5;     // 3 copies of ref in parents
    else if ( nref_mf==2 && fa==fb && ma==mb )
        gt_prior = 0.995 * 0.002 * (2.0/5.0) * (1.0/5.0) * 0.5;     // 2 copies of ref in parents, both homs
    else if ( nref_mf==2 )
        gt_prior = 0.995 * 0.002 * (2.0/5.0) * (2.0/5.0);           // 2 copies of ref in parents, both hets
    else if ( nref_mf==1 )
    {
        assert( nals_mf==2 && nals_mfc==2 );
        gt_prior = 0.995 * 0.002 * (2.0/5.0) * (2.0/5.0) * 0.5;     // 1 copy of ref in parents
    }
    else if ( nref_mf==0 )
    {
        if ( nals_mf==1 )
            gt_prior = 0.995 * 0.002 * (3.0/5.0) * (1.0/5.0);       // 1 alt allele in the trio
        else if ( nals_mf==2 )
        {
            assert( ca!=0 && cb!=0 );
            gt_prior = 0.002 * 0.002 / 414;                         // 2 alt alleles and 0 refs in the trio
        }
        else
            error("Fixme: %s:%d\n",__FILE__,__LINE__);
    }
    else
        error("Fixme: %s:%d\n",__FILE__,__LINE__);
    return gt_prior;
}
// Parent genotype probability L(GM,GF), with DNG bugs fixed
static double init_mf_priors(args_t *args, int fi, int mi)
{
    double gt_prior = 0;    // parent genotype probability L(GM,GF)
    int fa = seq1[fi];
    int fb = seq2[fi];
    int ma = seq1[mi];
    int mb = seq2[mi];
    int gts[3]; gts[0] = fi; gts[1] = mi; gts[2] = 0;
    int nalt_mf = count_unique_alleles(2,gts,only_alts);
    int nref_mf = (fa==0 ? 1 : 0) + (fb==0 ? 1 : 0) + (ma==0 ? 1 : 0) + (mb==0 ? 1 : 0);

    const double p_homref = 0.998;                                  // this assumes bi-allelic sites
    const double p_poly   = (1 - p_homref) * (1 - p_homref);        // p of this occuring twice for a different allele
    const double p_nonref = 1 - p_homref - p_poly;

    if ( nalt_mf>=3 )                                               // penalize heavily sites with 3 unique ALTs
        gt_prior = 1e-26;
    else if ( nalt_mf>=2 )                                          // 2 unique ALTs, 19*3 = 57 cases
        gt_prior = p_poly / 57.;
    else if ( nref_mf==4 )                                          // 0 ALTs; 00,00
        gt_prior = p_homref;
    else if ( nref_mf==3 )                                          // this and all remaining have 1 unique ALT allele; 00,0x 
        gt_prior = p_nonref * (4.0/15.0) * (1.0/3.0);
    else if ( nref_mf==2 && ma==mb )                                // hom alt; 00,xx
        gt_prior = p_nonref * (2.0/15.0) * (1.0/3.0);
    else if ( nref_mf==2 )                                          // two hets; 0x,0x
            gt_prior = p_nonref * (4.0/15.0) * (1.0/3.0);
    else if ( nref_mf==1 )                                          // single ref; 0x,xx
        gt_prior = p_nonref * (4.0/15.0) * (1.0/3.0);
    else if ( nref_mf==0 )                                          // no ref; xx,xx
        gt_prior = p_nonref * (1.0/15.0) * (1.0/3.0);
    else
        error("Fixme: %s:%d\n",__FILE__,__LINE__);
    return gt_prior;
}
static double init_mf_priors_chrX(args_t *args, int mi)
{
    double gt_prior = 0;    // parent genotype probability L(GM)
    int ma = seq1[mi];
    int mb = seq2[mi];
    int gts[3]; gts[0] = mi; gts[1] = 0; gts[2] = 0;
    int nalt_m = count_unique_alleles(1,gts,only_alts);
    int nref_m = (ma==0 ? 1 : 0) + (mb==0 ? 1 : 0);

    const double p_homref = 0.999;                                  // this assumes bi-allelic sites
    const double p_poly   = (1 - p_homref) * (1 - p_homref);        // p of this occuring twice for a different allele
    const double p_nonref = 1 - p_homref - p_poly;

    if ( nalt_m>=2 )                                               // 2 unique ALTs, 3 cases
        gt_prior = p_poly / 3.;
    else if ( nref_m==2 )                                          // 00
        gt_prior = p_homref;
    else if ( nref_m==1 )                                          // single ref; 0x and x0
        gt_prior = p_nonref * (2.0/3.0) * (1.0/3.0);
    else if ( nref_m==0 )                                          // no ref; xx,xx
        gt_prior = p_nonref * (1.0/3.0) * (1.0/3.0);
    else
        error("Fixme: %s:%d\n",__FILE__,__LINE__);
    return gt_prior;
}
static double init_mf_priors_chrXX(args_t *args, int fi, int mi)
{
    double gt_prior = 0;    // parent genotype probability L(GM)
    int fa = seq1[fi];
    int fb = seq2[fi];
    int ma = seq1[mi];
    int mb = seq2[mi];
    int gts[3]; gts[0] = fi; gts[1] = mi; gts[2] = 0;
    int nalt_mf = count_unique_alleles(2,gts,only_alts);
    int nref_mf = (fa==0 ? 1 : 0) + (fb==0 ? 1 : 0) + (ma==0 ? 1 : 0) + (mb==0 ? 1 : 0);
    if ( fa!=fb ) return 0;     // father can't be a het
    if ( fa==0 ) nref_mf--;
    else nalt_mf--;

    const double p_homref = 0.998;                                  // this assumes bi-allelic sites
    const double p_poly   = (1 - p_homref) * (1 - p_homref);        // p of this occuring twice for a different allele
    const double p_nonref = 1 - p_homref - p_poly;

    if ( nalt_mf>=3 )                                               // 3 unique ALTs
        gt_prior = 1e-26;
    else if ( nalt_mf>=2 )                                          // 2 unique ALTs
        gt_prior = p_poly * (1.0/9.0) * (1.0/3.0);
    else if ( nref_mf==3 )                                          // 00,0
        gt_prior = p_homref;
    else if ( nref_mf==2 )                                          // 00,x; 0x,0; x0,0
        gt_prior = p_nonref * (3.0/7.0) * (1.0/3.0);
    else if ( nref_mf==1 )                                          // 0x,x; x0,x; xx,0
        gt_prior = p_nonref * (3.0/7.0) * (1.0/3.0);
    else if ( nref_mf==0 )                                          // no ref; xx,x
        gt_prior = p_nonref * (1.0/7.0) * (1.0/3.0);
    else
        error("Fixme: %s:%d\n",__FILE__,__LINE__);
    return gt_prior;
}
static void init_DNG_tprob_mprob(args_t *args, int fi, int mi, int ci, double *tprob, double *mprob, int *denovo_allele)
{
    int fa = seq1[fi];
    int fb = seq2[fi];
    int ma = seq1[mi];
    int mb = seq2[mi];
    int gts[3]; gts[0] = fi; gts[1] = mi; gts[2] = 0;
    int ca = seq1[ci];
    int cb = seq2[ci];
    gts[0] = fi; gts[1] = mi; gts[2] = ci;
    int nals_mfc = count_unique_alleles(3,gts,include_ref);
    *tprob = 1;                   // genotype transmission likelihood L(GC|GM,GF), 0 if not compatible with Mendelian inheritance
    *mprob = 1 - args->mrate;     // probability of mutation
    *denovo_allele = ca!=fa && ca!=fb && ca!=ma && ca!=mb ? ca : cb;

    if ( nals_mfc==4 )
        *tprob = 0;                     // 4 unique alleles
    else if ( nals_mfc==3 )             // 3 alleles
    {
        if ( ((ca==fa || ca==fb) && (cb==ma || cb==mb)) ||
                ((cb==fa || cb==fb) && (ca==ma || ca==mb)) )
        {
            if ( ca==cb ) *tprob = 0.25;
            else if ( fa==fb || ma==mb ) *tprob = 0.5;   // one parent is homozygous
            else *tprob = 0.25;
        }
        else
        {
            if ( ca!=fa  && ca!=fb && ca!=ma  && ca!=mb && 
                    cb!=fa  && cb!=fb && cb!=ma  && cb!=mb ) *mprob = args->mrate * args->mrate;    // two mutations
            else
                *mprob = args->mrate;
            *tprob = 0;
        }
    }
    else if ( nals_mfc==2 )                      // 2 alleles
    {
        if ( fa!=fb && ma!=mb ) *tprob = 0.25;   // both parents are hets
        else if ( fa==fb && ma==mb )             // both parents are homs
        {
            if ( fa==ma && ca==cb ) *tprob = 0, *mprob = args->mrate * args->mrate;   // parents same homs, child a hom, two alleles mutated
            else if ( fa==ma ) *tprob = 0, *mprob = args->mrate;                      // parents same homs, child a het, one allele mutated
            else if ( ca==cb ) *tprob = 0, *mprob = args->mrate;                      // parents diff homs, child a hom, one allele mutated
        }
        else if ( ca==cb && ((fa==fb && fa!=ca) || (ma==mb && ma!=ca)) )
            *tprob = 0, *mprob = args->mrate;                                         // child is (wrong) hom and one parent is hom
        else
            *tprob = 0.5;
    }
}
static void init_tprob_mprob(args_t *args, int fi, int mi, int ci, double *tprob, double *mprob, int *denovo_allele)
{
    int fa = seq1[fi];
    int fb = seq2[fi];
    int ma = seq1[mi];
    int mb = seq2[mi];
    int ca = seq1[ci];
    int cb = seq2[ci];

    *denovo_allele = ca!=fa && ca!=fb && ca!=ma && ca!=mb ? ca : cb;

    // tprob .. genotype transmission probability L(GC|GM,GF), 0 if not compatible with Mendelian inheritance
    // mprob .. probability of mutation

    if ( ((ca==fa||ca==fb) && (cb==ma||cb==mb)) || ((ca==ma||ca==mb) && (cb==fa||cb==fb)) )
    {
        if ( fa==fb && ma==mb ) *tprob = 1;
        else if ( fa==fb || ma==mb ) *tprob = 0.5;
        else *tprob = 0.25;
        *mprob = 1 - args->mrate;
    }
    else
    {
        *tprob = 0;
        if ( (ca==fa||ca==fb) || (ca==ma||ca==mb) || (cb==fa||cb==fb) || (cb==ma||cb==mb) ) *mprob = args->mrate;
        else *mprob = args->mrate * args->mrate;
    }
}
static void init_tprob_mprob_chrX(args_t *args, int mi, int ci, double *tprob, double *mprob, int *denovo_allele)
{
    int ma = seq1[mi];
    int mb = seq2[mi];
    int ca = seq1[ci];
    int cb = seq2[ci];

    *denovo_allele = ca!=ma && ca!=mb ? ca : cb;

    if ( ca!=cb )                   // male cannot be heterozygous in X
        *mprob = 0, *tprob = 0;
    else if ( ca==ma || ca==mb )    // inherited
    {
        if ( ma==mb ) *tprob = 1;
        else *tprob = 0.5;
        *mprob = 1 - args->mrate;
    }
    else                            // de novo
        *mprob = args->mrate, *tprob = 0;
}
static void init_tprob_mprob_chrXX(args_t *args, int fi, int mi, int ci, double *tprob, double *mprob, int *denovo_allele)
{
    int fa = seq1[fi];
    int fb = seq2[fi];
    int ma = seq1[mi];
    int mb = seq2[mi];
    int ca = seq1[ci];
    int cb = seq2[ci];

    *denovo_allele = ca!=fa && ca!=fb && ca!=ma && ca!=mb ? ca : cb;

    if ( fa!=fb )                   // father cannot be heterozygous in X
        *mprob = 0, *tprob = 0;
    else if ( (ca==fa && (cb==ma||cb==mb)) || (cb==fa && (ca==ma||ca==mb)) )
    {
        if ( ma==mb ) *tprob = 1;
        else *tprob = 0.5;
        *mprob = 1 - args->mrate;
    }
    else
    {
        *tprob = 0;
        if ( (ca==fa || (ca==ma||ca==mb)) || (cb==fa || (cb==ma||cb==mb)) ) *mprob = args->mrate, *tprob = 0;
        else *mprob = args->mrate * args->mrate;
    }
}
typedef enum { autosomal, chrX, chrXX } init_priors_t;
static void init_priors(args_t *args, priors_t *priors, init_priors_t type)
{
    // Based on the FIGL model from the supplement "Variation in genome-wide mutation rates within and between human families"
    int fi,mi,ci;
    for (fi=0; fi<10; fi++)
    {
        for (mi=0; mi<10; mi++)
        {
            for (ci=0; ci<10; ci++)
            {
                double gt_prior;                // parent genotype probability L(GM,GF)
                double tprob;                   // genotype transmission likelihood L(GC|GM,GF), 0 if not compatible with Mendelian inheritance
                double mprob;                   // probability of mutation
                int allele;                     // which of the alleles is de novo
                if ( args->use_dng_priors )
                    gt_prior = init_DNG_mf_priors(args,fi,mi,ci);
                else if ( type==autosomal )
                    gt_prior = init_mf_priors(args,fi,mi);
                else if ( type==chrX )
                    gt_prior = init_mf_priors_chrX(args,mi);
                else if ( type==chrXX )
                    gt_prior = init_mf_priors_chrXX(args,fi,mi);
                else
                    error("Can't happen\n");

                if ( args->use_dng_priors )
                    init_DNG_tprob_mprob(args,fi,mi,ci,&tprob,&mprob,&allele);
                else if ( type==autosomal )
                    init_tprob_mprob(args,fi,mi,ci,&tprob,&mprob,&allele);
                else if ( type==chrX )
                    init_tprob_mprob_chrX(args,mi,ci,&tprob,&mprob,&allele);
                else if ( type==chrXX )
                    init_tprob_mprob_chrXX(args,fi,mi,ci,&tprob,&mprob,&allele);
                else
                    error("Can't happen\n");

                priors->denovo_allele[fi][mi][ci] = tprob==0 ? allele : UINT8_MAX;  // the latter should never happen, making it fail deliberately
                priors->denovo[fi][mi][ci] = tprob==0 ? 1 : 0;
                priors->pprob[fi][mi][ci]  = log(gt_prior * mprob * (tprob==0 ? 1 : tprob));
            }
        }
    }
}
static void init_data(args_t *args)
{
    char *ptr = strchr(args->dnm_score_tag,':');
    if ( ptr )
    {
        if ( ptr==args->dnm_score_tag ) error("Error: could not parse --use tag=%s\n",ptr);
        *ptr = 0;
        if ( !strcasecmp(ptr+1,"log") ) args->dnm_score_is_float = 1;
        else if ( strcasecmp(ptr+1,"phred") ) error("Error: the type \"%s\" is not supported --use tag\n",ptr+1);
    }

    args->sr = bcf_sr_init();
    if ( args->regions )
    {
        args->sr->require_index = 1;
        if ( bcf_sr_set_regions(args->sr, args->regions, args->regions_is_file)<0 ) error("Failed to read the regions: %s\n",args->regions);
    }
    if ( args->targets && bcf_sr_set_targets(args->sr, args->targets, args->targets_is_file, 0)<0 ) error("Failed to read the targets: %s\n",args->targets);
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr = bcf_sr_get_header(args->sr,0);

    if ( args->filter_str ) args->filter = filter_init(args->hdr, args->filter_str);

    int id;
    if ( (id=bcf_hdr_id2int(args->hdr, BCF_DT_ID, "PL"))<0 || !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,id) )
        error("Error: the tag FORMAT/PL is not present in %s\n", args->fname);
    if ( (args->use_model&USE_ACM) && ((id=bcf_hdr_id2int(args->hdr, BCF_DT_ID, "QS"))<0 || !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,id)) && !args->use_ppl )
        error(
            "Error:\n"
            "   The FORMAT/QS tag is not present. If you want to proceed anyway, run with the `--use ppl`\n"
            "   option at the cost of inflated false discovery rate. The QS annotation can be generated\n"
            "   at the mpileup step together with the AD annotation using the command\n"
            "       bcftools mpileup -a AD,QS -f ref.fa file.bam\n");   // Possible future todo: use AD as a proxy for QS?
    if ( (id=bcf_hdr_id2int(args->hdr, BCF_DT_ID, "AD"))<0 || !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,id) )
        fprintf(stderr, "Warning: the tag FORMAT/AD is not present in %s, the output tag FORMAT/VAF will not be added\n", args->fname);
    else
        args->has_fmt_ad = 1;

    init_priors(args,&args->priors,autosomal);
    init_priors(args,&args->priors_X,chrX);
    init_priors(args,&args->priors_XX,chrXX);

    args->hdr_out = bcf_hdr_dup(args->hdr);
    bcf_hdr_printf(args->hdr_out, "##FORMAT=<ID=%s,Number=1,Type=%s,Description=\"De-novo mutation score, bigger values = bigger confidence\">",args->dnm_score_tag,args->dnm_score_is_float?"Float":"Integer");
    bcf_hdr_printf(args->hdr_out, "##FORMAT=<ID=%s,Number=1,Type=Integer,Description=\"The de-novo allele\">",args->dnm_allele_tag);
    if ( args->has_fmt_ad )
        bcf_hdr_printf(args->hdr_out, "##FORMAT=<ID=%s,Number=1,Type=Integer,Description=\"The percentage of ALT reads\">",args->dnm_vaf_tag);

    int i, n = 0;
    char **list;
    if ( args->pfm )
    {
        args->ntrio = 1;
        args->trio  = (trio_t*) calloc(1,sizeof(trio_t));
        list = hts_readlist(args->pfm, 0, &n);
        if ( n!=3 ) error("Expected three sample names with -t\n");
        args->trio[0].idx[iCHILD]  = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, list[0]);
        args->trio[0].idx[iFATHER] = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, list[1]);
        args->trio[0].idx[iMOTHER] = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, list[2]);
        if ( args->trio[0].idx[iCHILD] < 0 )
        {
            if ( strlen(list[0])>3 && !strncasecmp(list[0],"1X:",3) )
            {
                args->trio[0].idx[iCHILD] = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, list[0]+3);
                args->trio[0].is_male = 1;
            }
            else if ( strlen(list[0])>3 && !strncasecmp(list[0],"2X:",3) )
                args->trio[0].idx[iCHILD] = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, list[0]+3);
        }
        for (i=0; i<n; i++)
        {
            if ( args->trio[0].idx[i] < 0 ) error("The sample is not present: %s\n", list[i]);
            free(list[i]);
        }
        free(list);
    }
    else
    {
        parse_ped(args,args->ped_fname);
        if ( !args->ntrio ) error("No complete trio present\n");
    }
    if ( !args->chrX_list_str || !strcasecmp("GRCh37",args->chrX_list_str) )
        args->chrX_list_str = "X:1-60000,chrX:1-60000,X:2699521-154931043,chrX:2699521-154931043";
    else if ( !strcasecmp("GRCh38",args->chrX_list_str) )
        args->chrX_list_str = "X:1-9999,chrX:1-9999,X:2781480-155701381,chrX:2781480-155701381";
    char *rmme = strdup(args->chrX_list_str), *tmp = rmme;
    while ( *tmp )
    {
        if ( *tmp==',' ) *tmp = '\n';
        tmp++;
    }
    args->chrX_idx = regidx_init_string(rmme, regidx_parse_reg, NULL, 0, NULL);
    free(rmme);

    args->out_fh = hts_open(args->output_fname,hts_bcf_wmode2(args->output_type,args->output_fname));
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( bcf_hdr_write(args->out_fh, args->hdr_out)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);

    if ( args->dnm_score_is_float )
        args->dnm_qual_float = (float*) malloc(sizeof(*args->dnm_qual_float)*bcf_hdr_nsamples(args->hdr));
    else
        args->dnm_qual_int = (int32_t*) malloc(sizeof(*args->dnm_qual_int)*bcf_hdr_nsamples(args->hdr));
    args->vaf = (int32_t*) malloc(sizeof(*args->vaf)*bcf_hdr_nsamples(args->hdr));
    args->dnm_allele = (int32_t*) malloc(sizeof(*args->dnm_allele)*bcf_hdr_nsamples(args->hdr));
}
static void destroy_data(args_t *args)
{
    if ( args->filter ) filter_destroy(args->filter);
    regidx_destroy(args->chrX_idx);
    free(args->dnm_score_tag);
    free(args->dnm_vaf_tag);
    free(args->dnm_allele_tag);
    free(args->pl3);
    free(args->aprob);
    free(args->idx);
    free(args->dnm_qual_int);
    free(args->dnm_qual_float);
    free(args->dnm_allele);
    free(args->vaf);
    free(args->trio);
    free(args->pl);
    free(args->ad);
    free(args->qs);
    free(args->qs3);
    if ( hts_close(args->out_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname);
    bcf_hdr_destroy(args->hdr_out);
    bcf_sr_destroy(args->sr);
    free(args);
}

static inline double phred2num(double phred)
{
    return pow(10,-0.1*phred);
}
static inline double log2phred(double num)
{
    return fabs(4.3429 * num);
}
static inline double phred2log(double phred)
{
    return -phred/4.3429;
}
#if 0
static inline double subtract_num_log(double a_num, double b_log)
{
    return log(a_num - exp(b_log));
}
#endif
static inline double subtract_log(double a_log, double b_log)
{
    if ( b_log==-HUGE_VAL ) return a_log;
    return log(exp(a_log - b_log) - 1) + b_log;
}
static inline double sum_log(double a, double b)    // log(exp(a)+exp(b))
{
    if ( a==-HUGE_VAL && b==-HUGE_VAL ) return -HUGE_VAL;
    if ( a>b )
        return log(1 + exp(b-a)) + a;
    else
        return log(1 + exp(a-b)) + b;
}
static double process_trio_ACM(args_t *args, priors_t *priors, int nals, double *pl[3], int npl, double *qs[3], int *al0, int *al1)
{
    assert( nals>1 && nals<=4 );

    *al0 = *al1 = 0;

    double sum = -HUGE_VAL, max = -HUGE_VAL;
    int i, ca,cb, fa,fb, ma,mb, ci=0;
    for (ca=0; ca<nals; ca++)
    {
        for (cb=0; cb<=ca; cb++)
        {
            int cals = (1<<ca)|(1<<cb);
            double cpl = pl[iCHILD][ci];
            int fi = 0;
            for (fa=0; fa<nals; fa++)
            {
                for (fb=0; fb<=fa; fb++)
                {
                    int fals = (1<<fa)|(1<<fb);
                    double fpl;
                    if ( args->use_ppl || args->use_ppl_qs ) fpl = pl[iFATHER][fi];
                    else
                    {
                        fpl = 0;
                        for (i=0; i<4; i++)
                        {
                            if ( fals&(1<<i) )
                                fpl += subtract_log(0,qs[iFATHER][i]);
                            else if ( cals&(1<<i) )
                                fpl += qs[iFATHER][i];
                            else if ( fa==fb )
                                fpl += qs[iFATHER][i];
                        }
                    } 
                    int mi = 0;
                    for (ma=0; ma<nals; ma++)
                    {
                        for (mb=0; mb<=ma; mb++)
                        {
                            int mals = (1<<ma)|(1<<mb);
                            double mpl = 0;
                            if ( args->use_ppl || args->use_ppl_qs ) mpl = pl[iMOTHER][mi];
                            else
                            {
                                mpl = 0;
                                for (i=0; i<4; i++)
                                {
                                    if ( mals&(1<<i) )
                                        mpl += subtract_log(0,qs[iMOTHER][i]);
                                    else if ( cals&(1<<i) )
                                        mpl += qs[iMOTHER][i];
                                    else if ( ma==mb )
                                        mpl += qs[iMOTHER][i];
                                }
                            }
                            double val = cpl + fpl + mpl + priors->pprob[fi][mi][ci];
                            sum = sum_log(sum,val);
#define DEBUG 0
#if DEBUG
                            if(val!=-HUGE_VAL)                            
                                fprintf(stderr,"m,f,c: %d%d+%d%d=%d%d  dn=%d (%d,%d,%d)   mpl,fpl,cpl: %+e %+e %+e \t prior:%+e \t pval=%+e  sum=%+e  %c\n",
                                    mb,ma,fb,fa,cb,ca,priors->denovo[fi][mi][ci],fi,mi,ci,mpl,fpl,cpl,priors->pprob[fi][mi][ci], val,sum,(priors->denovo[fi][mi][ci] && max < val)?'*':'-');
#endif
                            if ( priors->denovo[fi][mi][ci] && max < val )
                            {
                                max = val;
                                if ( priors->denovo_allele[fi][mi][ci] == ca )
                                    *al0 = cb, *al1 = ca;
                                else
                                    *al0 = ca, *al1 = cb;
                            }
                            mi++;
                        }
                    }
                    fi++;
                }
            }
            ci++;
        }
    }
#if DEBUG
    fprintf(stderr,"max=%e sum=%e   ret=%e\n",max,sum,max-sum);
#endif
    return log2phred(subtract_log(0,max-sum));
}
static double process_trio_DNG(args_t *args, priors_t *priors, int nals, double *pl[3], int npl, int *al0, int *al1)
{
    assert( nals>1 && nals<=4 );

    *al0 = *al1 = 0;

    double sum = -HUGE_VAL, max = -HUGE_VAL;
    int ca,cb, fa,fb, ma,mb, ci=0;
    for (ca=0; ca<nals; ca++)
    {
        for (cb=0; cb<=ca; cb++)
        {
            int fi = 0;
            for (fa=0; fa<nals; fa++)
            {
                for (fb=0; fb<=fa; fb++)
                {
                    int mi = 0;
                    for (ma=0; ma<nals; ma++)
                    {
                        for (mb=0; mb<=ma; mb++)
                        {
                            double val;
                            val = pl[iCHILD][ci] + pl[iFATHER][fi] + pl[iMOTHER][mi] + priors->pprob[fi][mi][ci];
                            sum = sum_log(val,sum);
#if DEBUG
                            if(val!=-HUGE_VAL)                            
                                fprintf(stderr,"m,f,c: %d%d+%d%d=%d%d  dn=%d (%d,%d,%d)   mpl,fpl,cpl: %+e %+e %+e \t prior:%+e \t pval=%+e  sum=%+e  %c\n",
                                    mb,ma,fb,fa,cb,ca,priors->denovo[fi][mi][ci],fi,mi,ci,pl[iMOTHER][mi],pl[iFATHER][fi],pl[iCHILD][ci],priors->pprob[fi][mi][ci], val,sum,(priors->denovo[fi][mi][ci] && max < val)?'*':'-');
#endif
                            if ( priors->denovo[fi][mi][ci] && max < val )
                            {
                                max = val;
                                if ( priors->denovo_allele[fi][mi][ci] == ca )
                                    *al0 = cb, *al1 = ca;
                                else
                                    *al0 = ca, *al1 = cb;
                            }
                            mi++;
                        }
                    }
                    fi++;
                }
            }
            ci++;
        }
    }
#if DEBUG
    fprintf(stderr,"max=%e sum=%e   ret=%e\n",max,sum,max-sum);
#endif
    return log2phred(subtract_log(0,max-sum));
}
static inline void qs_to_pl(args_t *args, double *qs, int nqs, double *pl, int npl)
{
    int i,j,k  = 0;
    double sum = 0;
    for (i=0; i<nqs; i++) sum += qs[i];
    for (i=0; i<nqs; i++)
    {
        for (j=0; j<=i; j++)
        {
            if ( i==j ) pl[k] += sum - qs[i];
            k++;
        }
    }
}
static int test_filters(args_t *args, bcf1_t *rec)
{
    uint8_t *smpl_pass;
    int i,j, pass_site = filter_test(args->filter, rec, (const uint8_t**) &smpl_pass);
    if ( args->filter_logic & FLT_EXCLUDE )
    {
        if ( pass_site )
        {
            if ( !smpl_pass ) return 0;     // no samples, -e mode, the expression failed
            pass_site = 0;
            for (i=0; i<args->ntrio; i++)
            {
                int pass_trio = 1;
                for (j=0; j<3; j++)
                {
                    int idx = args->trio[i].idx[j];
                    if ( smpl_pass[idx] ) { pass_trio = 0; break; }     // with -e one sample passes, the whole trio fails
                }
                args->trio[i].pass = pass_trio;
                if ( pass_trio ) pass_site = 1;
            }
            return pass_site;
        }
        for (i=0; i<args->ntrio; i++) args->trio[i].pass = 1;
        return 1;
    }
    if ( !pass_site ) return 0;
    if ( smpl_pass )
    {
        pass_site = 0;
        for (i=0; i<args->ntrio; i++)
        {
            int pass_trio = 1;
            for (j=0; j<3; j++)
            {
                int idx = args->trio[i].idx[j];
                if ( !smpl_pass[idx] ) { pass_trio = 0; break; }
            }
            args->trio[i].pass = pass_trio;
            if ( pass_trio ) pass_site = 1;
        }
        return pass_site;
    }
    for (i=0; i<args->ntrio; i++) args->trio[i].pass = 1;
    return 1;
}
static void process_record(args_t *args, bcf1_t *rec)
{
    int skip_site = 0;
    if ( rec->n_allele==1 || bcf_get_variant_types(rec)==VCF_REF ) skip_site = 1;
    else if ( args->filter && !test_filters(args,rec) ) skip_site = 1;
    if ( skip_site )
    {
        if ( bcf_write(args->out_fh, args->hdr_out, rec)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
        return;
    }

    static int n_ad_warned = 0;
    int nret, nsmpl = bcf_hdr_nsamples(args->hdr), n_ad = args->has_fmt_ad;
    if ( n_ad )
    {
        nret = bcf_get_format_int32(args->hdr,rec,"AD",&args->ad,&args->mad);
        if ( nret<=0 ) n_ad = 0;
        else
        {
            n_ad = nret / nsmpl;
            if ( nret != nsmpl * rec->n_allele )
            {
                if ( !n_ad_warned )
                {
                    hts_log_warning("Incorrect number of fields for FORMAT/AD at %s:%"PRId64". This warning is printed only once", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
                    n_ad_warned = 1;
                }
                if ( !args->force_ad ) n_ad = 0;
            }
        }
    }

    nret = bcf_get_format_int32(args->hdr,rec,"PL",&args->pl,&args->mpl);
    if ( nret<=0 ) error("The FORMAT/PL tag not present at %s:%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
    int npl1  = nret/nsmpl;
    if ( npl1!=rec->n_allele*(rec->n_allele+1)/2 )
        error("todo: not a diploid site at %s:%"PRId64": %d alleles, %d PLs\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1,rec->n_allele,npl1);
    hts_expand(double,3*npl1,args->mpl3,args->pl3);

    int nqs1 = 0;
    if ( args->use_model&USE_ACM && !args->use_ppl )
    {
        nret = bcf_get_format_int32(args->hdr,rec,"QS",&args->qs,&args->mqs);
        if ( nret<0 ) error("Error: the FMT/QS tag is not available at %s:%"PRId64".\n",bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        if ( nret != nsmpl * rec->n_allele ) error("Error: incorrect number of FMT/QS values at %s:%"PRId64".\n",bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        nqs1 = nret<=0 ? 0 : nret/nsmpl;
        hts_expand(double,3*nqs1,args->mqs3,args->qs3);
    }

    int is_chrX = 0;
    if ( regidx_overlap(args->chrX_idx,bcf_seqname(args->hdr,rec),rec->pos,rec->pos+rec->rlen,NULL) ) is_chrX = 1;

    int i, j, k, al0, al1, write_dnm = 0, ad_set = 0;
    if ( args->dnm_score_is_float )
        for (i=0; i<nsmpl; i++) bcf_float_set_missing(args->dnm_qual_float[i]);
    else
        for (i=0; i<nsmpl; i++) args->dnm_qual_int[i] = bcf_int32_missing;
    for (i=0; i<nsmpl; i++) args->dnm_allele[i] = bcf_int32_missing;
    for (i=0; i<args->ntrio; i++)
    {
        if ( args->filter && !args->trio[i].pass ) continue;

        // Samples can be in any order in the VCF, set PL and QS to reflect the iFATHER,iMOTHER,iCHILD indices
        double *ppl[3];
        double *pqs[3];
        for (j=0; j<3; j++) // set trio PLs
        {
            int32_t *src = args->pl + npl1 * args->trio[i].idx[j];      // j loops over iFATHER,iMOTHER,iCHILD
            double *dst = ppl[j] = args->pl3 + j*npl1;
            double sum = 0;
            for (k=0; k<npl1; k++) { dst[k] = phred2num(src[k]); sum += dst[k]; }
            for (k=0; k<npl1; k++) dst[k] = log(dst[k]/sum);
        }
        if ( args->use_model&USE_ACM )   // set trio QS
        {
            int32_t *ad_f = NULL, *ad_m = NULL;
            if ( args->pnoise_strict && args->ad )
            {
                // apply noise tolerance for alleles observed in a single parent only
                ad_f = args->ad + n_ad * args->trio[i].idx[iFATHER];
                ad_m = args->ad + n_ad * args->trio[i].idx[iMOTHER];
            }
            for (j=0; j<3; j++)
            {
                int32_t *ad = (args->pnoise_abs && args->ad ) ? args->ad + n_ad * args->trio[i].idx[j] : NULL;
                int32_t *qs = args->qs + nqs1 * args->trio[i].idx[j];
                pqs[j] = args->qs3 + j*nqs1;
                double noise_tolerance = 0;
                double sum_qs = 0, sum_ad = 0;
                if ( j!=iCHILD )
                {
                    for (k=0; k<nqs1; k++) sum_qs += qs[k];
                    noise_tolerance = sum_qs * args->pnoise_frac;
                    if ( ad )
                    {
                        for (k=0; k<n_ad; k++) sum_ad += ad[k];
                        if ( args->pnoise_abs )
                        {
                            if ( noise_tolerance < args->pnoise_abs * sum_qs / sum_ad )
                                noise_tolerance = args->pnoise_abs * sum_qs / sum_ad;
                        }
                    }
                }
                for (k=0; k<nqs1; k++)
                {
                    double val = qs[k];
                    if ( !args->pnoise_strict || !ad_f[k] || !ad_m[k] ) val -= noise_tolerance;
                    if ( val < 0 ) val = 0;
                    if ( val > 255 ) val = 255;
                    pqs[j][k] = phred2log(val);
                }
            }
            if ( args->use_ppl_qs )
            {
                qs_to_pl(args, pqs[iMOTHER], nqs1, ppl[iMOTHER], npl1);
                qs_to_pl(args, pqs[iFATHER], nqs1, ppl[iFATHER], npl1);
            }
        }
        priors_t *priors;
        if ( !is_chrX ) priors = &args->priors;
        else if ( args->trio[i].is_male ) priors = &args->priors_X;
        else priors = &args->priors_XX;

        double score;
        if ( args->use_model==USE_ACM ) score = process_trio_ACM(args, priors, rec->n_allele, ppl, npl1, pqs, &al0, &al1);
        else if ( args->use_model==USE_DNG ) score = process_trio_DNG(args, priors, rec->n_allele, ppl, npl1, &al0, &al1);
        else error("Uh, this should not happen\n");

        if ( score >= args->min_score )
        {
            write_dnm = 1;
            if ( args->dnm_score_is_float )
                args->dnm_qual_float[ args->trio[i].idx[iCHILD] ] = score==HUGE_VAL ? 0 : subtract_log(0,phred2log(score));
            else
            {
                if ( score>255 ) score = 255;
                args->dnm_qual_int[ args->trio[i].idx[iCHILD] ] = round(score);
            }
            args->dnm_allele[ args->trio[i].idx[iCHILD] ] = al1;
        }

        if ( n_ad )
        {
            if ( al0 < n_ad && al1 < n_ad )
            {
                ad_set = 1;
                for (j=0; j<3; j++)
                {
                    int32_t *src = args->ad + n_ad * args->trio[i].idx[j];
                    args->vaf[ args->trio[i].idx[j] ] = src[al0]+src[al1] ? round(src[al1]*100./(src[al0]+src[al1])) : 0;
                }
            }
            else
                for (j=0; j<3; j++) args->vaf[ args->trio[i].idx[j] ] = bcf_int32_missing;
        }
    }
    if ( write_dnm )
    {
        int ret;
        if ( args->dnm_score_is_float )
            ret = bcf_update_format_float(args->hdr_out,rec,args->dnm_score_tag,args->dnm_qual_float,nsmpl);
        else
            ret = bcf_update_format_int32(args->hdr_out,rec,args->dnm_score_tag,args->dnm_qual_int,nsmpl);
        if ( ret )
            error("Failed to write FORMAT/%s at %s:%"PRId64"\n", args->dnm_score_tag, bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        ret = bcf_update_format_int32(args->hdr_out,rec,args->dnm_allele_tag,args->dnm_allele,nsmpl);
        if ( ret )
            error("Failed to write FORMAT/%s at %s:%"PRId64"\n", args->dnm_allele_tag,bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        if ( ad_set )
        {
            if ( bcf_update_format_int32(args->hdr_out,rec,args->dnm_vaf_tag,args->vaf,nsmpl)!=0 )
                error("Failed to write FORMAT/%s at %s:%"PRId64"\n", args->dnm_vaf_tag,bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        }
    }
    if ( bcf_write(args->out_fh, args->hdr_out, rec)!=0 ) error("[%s] Error: cannot write to %s at %s:%"PRId64"\n", __func__,args->output_fname,bcf_seqname(args->hdr,rec),(int64_t)rec->pos+1);
}

static void set_option(args_t *args, char *optarg)
{
    char *tmp;
    char *opt = strdup(optarg);
    char *val = strchr(opt,'=');
    if ( val ) { *val = 0; val++; }
    if ( !strcasecmp(opt,"mrate") )
    {
        if ( !val ) error("Error: expected value with -u mrate, e.g. -u mrate=1e-8\n");
        args->mrate = strtod(val,&tmp);
        if ( *tmp ) error("Could not parse: -u %s\n", optarg);
    }
    else if ( !strcasecmp(opt,"pn") || !strcasecmp(opt,"pnoise") || !strcasecmp(opt,"pns") )
    {
        if ( !val ) error("Error: expected value with -u %s, e.g. -u %s=0.05\n",opt,opt);
        args->pnoise_frac = strtod(val,&tmp);
        if ( *tmp && *tmp==',' )
        {
            args->pnoise_abs = strtod(tmp+1,&tmp);
            if ( *tmp ) error("Could not parse: -u %s\n", optarg);
        }
        if ( args->pnoise_frac<0 || args->pnoise_frac>1 ) error("Error: expected value from the interval [0,1] for -u %s\n", optarg);
        if ( args->pnoise_abs<0 ) error("Error: expected positive value for -u %s\n", optarg);
        args->pnoise_strict = !strcasecmp(opt,"pn") ? 0 : 1;
    }
    else if ( !strcasecmp(opt,"DNG") ) { args->use_model = USE_DNG; args->use_dng_priors = 1; }
    else if ( !strcasecmp(opt,"dng-priors") ) args->use_dng_priors = 1;
    else if ( !strcasecmp(opt,"ppl") ) args->use_ppl = 1;
    else if ( !strcasecmp(opt,"tag") )
    {
        if ( !val ) error("Error: expected value with -u tag, e.g. -u tag=ANN\n");
        free(args->dnm_score_tag);
        args->dnm_score_tag = strdup(val);
    }
    else if ( !strcasecmp(opt,"vaf") )
    {
        if ( !val ) error("Error: expected value with -u vaf, e.g. -u vaf=VAF\n");
        free(args->dnm_vaf_tag);
        args->dnm_vaf_tag = strdup(val);
    }
    else if ( !strcasecmp(opt,"va") )
    {
        if ( !val ) error("Error: expected value with -u va, e.g. -u va=VA\n");
        free(args->dnm_allele_tag);
        args->dnm_allele_tag = strdup(val);
    }
    else error("Error: the option \"-u %s\" is not recognised\n",optarg);
    free(opt);
}
int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_fname = "-";
    args->dnm_score_tag  = strdup("DNM:phred");
    args->dnm_vaf_tag    = strdup("VAF");
    args->dnm_allele_tag = strdup("VA");
    args->mrate = 1e-8;
    args->pnoise_frac   = 0.045;
    args->pnoise_abs    = 0;
    args->pnoise_strict = 1;
    args->use_model = USE_ACM;
    args->dnm_score_is_float = 1;
    static struct option loptions[] =
    {
        {"chrX",required_argument,0,'X'},
        {"use",required_argument,0,'u'},
        {"force-AD",no_argument,0,1},
        {"min-score",required_argument,0,'m'},
        {"include",required_argument,0,'i'},
        {"exclude",required_argument,0,'e'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"ped",required_argument,NULL,'P'},
        {"pfm",required_argument,NULL,'p'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "p:P:o:O:s:i:e:r:R:t:T:m:au:X:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'X': args->chrX_list_str = optarg; break;
            case 'u': set_option(args,optarg);
            case  1 : args->force_ad = 1; break;
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 't': args->targets = optarg; break;
            case 'T': args->targets = optarg; args->targets_is_file = 1; break;
            case 'r': args->regions = optarg; break;
            case 'R': args->regions = optarg; args->regions_is_file = 1; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': args->output_type = FT_BCF_GZ; break;
                          case 'u': args->output_type = FT_BCF; break;
                          case 'z': args->output_type = FT_VCF_GZ; break;
                          case 'v': args->output_type = FT_VCF; break;
                          default: error("The output type \"%s\" not recognised\n", optarg);
                      };
                      break;
            case 'P': args->ped_fname = optarg; break;
            case 'p': args->pfm = optarg; break;
            case 'm': args->min_score = strtod(optarg,&tmp);
                      if ( *tmp ) error("Could not parse: -M, --min-score %s\n", optarg);
                      break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error("%s", usage_text()); }
    }
    else if ( optind+1!=argc ) error("%s", usage_text());
    else args->fname = argv[optind];

    if ( !args->ped_fname && !args->pfm ) error("Missing the -p or -P option\n");
    if ( args->ped_fname && args->pfm ) error("Expected only -p or -P option, not both\n");

    init_data(args);

    while ( bcf_sr_next_line(args->sr) )
        process_record(args, bcf_sr_get_line(args->sr,0));

    destroy_data(args);

    return 0;
}
