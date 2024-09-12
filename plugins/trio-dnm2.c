/* The MIT License

   Copyright (c) 2018-2024 Genome Research Ltd.

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
#define USE_NAIVE  3     // a naive calling model based on observed GT

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define iFATHER 0   // don't modify, QS calculations depend on this order!
#define iMOTHER 1
#define iCHILD  2

// output tag type
#define DNM_INT   (1<<0)
#define DNM_FLOAT (1<<1)
#define DNM_LOG   ((1<<2)|DNM_FLOAT)
#define DNM_PHRED ((1<<3)|DNM_INT)
#define DNM_PROB  ((1<<4)|DNM_FLOAT)
#define DNM_FLAG  ((1<<5)|DNM_INT)

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
    uint8_t denovo[10][10][10];         // is the GT combination not compatible with normal inheritance (0) or is de novo (1)
    uint8_t denovo_allele[10][10][10];  // which of the alleles is de novo for this configuration
}
priors_t;

typedef struct
{
    double abs, frac;
    double abs1, frac1;   // applied only if allele observed in a single parent, but not when observed in both
}
pnoise_t;

typedef struct
{
    int argc, filter_logic, regions_is_file, targets_is_file, output_type, record_cmd_line, clevel;
    int regions_overlap, targets_overlap;
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
    int32_t *pl,*ad,*qs,*gt, *dnm_qual_int, *dnm_allele, *vaf;    // input FMT/PL,AD,QS,GT values, output DNM and VAF
    float *dnm_qual_float;
    int mpl, mad, mqs, mgt;
    double min_score;
    double *pl3;    // normalized PLs converted to probs for iFATHER,iMOTHER,iCHILD
    double *qs3;    // QS converted to probs for iFATHER,iMOTHER,iCHILD
    int maprob, mpl3, mqs3, midx, *idx, force_ad, use_model;
    double *alt_tmp;
    int *alt_idx;
    int malt_tmp, malt_idx;
    char *dnm_score_tag,            // the argument of --dnm-tag, by default DNM:log
         *dnm_vaf_tag,
         *dnm_allele_tag;
    int dnm_score_type;             // given by e.g. --dnm-tag DNM:log
    double mrate;                   // --mrate, mutation rate
    pnoise_t pn_snv, pn_indel;      // --pn and --pns for SNVs and indels
    int with_ppl, with_pad;         // --with-pPL or --with-pAD
    int use_dng_priors;             // --dng-priors
    int need_QS;
    int strictly_novel;
    priors_t priors, priors_X, priors_XX;
    char *index_fn;
    int write_index;
}
args_t;

const char *about(void)
{
    return "Screen variants for possible de-novo mutations in trios.\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: Screen variants for possible de-novo mutations in trios\n"
        "Usage: bcftools +trio-dnm2 [OPTIONS]\n"
        "Common options:\n"
        "   -e, --exclude EXPR              Exclude trios for which the expression is true (one matching sample invalidates a trio)\n"
        "   -i, --include EXPR              Include trios for which the expression is true (one failing samples invalidates a trio)\n"
        "   -o, --output FILE               Output file name [stdout]\n"
        "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n"
        "   -r, --regions REG               Restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE         Restrict to regions listed in a file\n"
        "       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n"
        "   -t, --targets REG               Similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE         Similar to -R but streams rather than index-jumps\n"
        "       --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n"
        "       --no-version                Do not append version and command line to the header\n"
        "\n"
        "General options:\n"
        "   -m, --min-score NUM             Do not add FMT/DNM annotation if the score is smaller than NUM\n"
        "   -p, --pfm [1X:|2X:]P,F,M        Sample names of child (the proband), father, mother; \"1X:\" for male pattern of chrX inheritance [2X:]\n"
        "   -P, --ped FILE                  PED file with the columns: <ignored>,proband,father,mother,sex(1:male,2:female)\n"
        "   -X, --chrX LIST                 List of regions with chrX inheritance pattern or one of the presets: [GRCh37]\n"
        "                                      GRCh37 .. X:1-60000,chrX:1-60000,X:2699521-154931043,chrX:2699521-154931043\n"
        "                                      GRCh38 .. X:1-9999,chrX:1-9999,X:2781480-155701381,chrX:2781480-155701381\n"
        "       --dnm-tag TAG[:type]        Output tag with DNM quality score and its type [DNM:log]\n"
        "                                       log   .. log-scaled quality (-inf,0; float)\n"
        "                                       flag  .. is a DNM, implies --use-NAIVE (1; int)\n"
        "                                       phred .. phred quality (0-255; int)\n"
        "                                       prob  .. probability (0-1; float)\n"
        "       --force-AD                  Calculate VAF even if the number of FMT/AD fields is incorrect. Use at your own risk!\n"
        "       --va TAG                    Output tag name for the variant allele [VA]\n"
        "       --vaf TAG                   Output tag name for variant allele fraction [VAF]\n"
        "\n"
        "Model options:\n"
        "       --dng-priors                Use the original DeNovoGear priors (including bugs in prior assignment, but with chrX bugs fixed)\n"
        "       --mrate NUM                 Mutation rate [1e-8]\n"
        "       --pn FRAC[,NUM][:type]      Tolerance to parental noise or mosaicity, given as fraction of QS or number of reads. The type is\n"
        "                                       'snv', 'indel' or 'both' [--pn 0.005,0:snv --pn 0,0:indel]\n"
        "       --pns FRAC[,NUM][:type]     Same as --pn but applied for alleles observed in a single parent [--pns 0.045,0:snv --pns 0,0:indel]\n"
        "   -n, --strictly-novel            When Mendelian inheritance is violiated, score highly only novel alleles (e.g. in LoH regions)\n"
        "       --use-DNG                   The original DeNovoGear model, implies --dng-priors\n"
        "       --use-NAIVE                 A naive calling model which uses only FMT/GT to determine DNMs\n"
        "       --with-pAD                  Do not use FMT/QS but parental FMT/AD\n"
        "       --with-pPL                  Do not use FMT/QS but parental FMT/PL. Equals to DNG with bugs fixed (more FPs, fewer FNs)\n"
        "   -W, --write-index[=FMT]         Automatically index the output files [off]\n"
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
        "       bcftools call -mv -Ou |\n"
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
static const int8_t seq3[13]  = {-1,0,2,1,5,3,4,-1,9,6,7,-1,8};     // Lookup from (1<<ial)|(1<<jal) to iseq
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
    const double p_poly   = (1 - p_homref) * (1 - p_homref);        // p of this occurring twice for a different allele
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
    const double p_poly   = (1 - p_homref) * (1 - p_homref);        // p of this occurring twice for a different allele
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
    const double p_poly   = (1 - p_homref) * (1 - p_homref);        // p of this occurring twice for a different allele
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

    int is_novel;
    if ( args->strictly_novel )
    {
        // Consider as DNMs only variants with a novel allele. For example, if heterozygosity is lost and
        // Mendelian inheritance is violated (11+00=11), don't consider this to be a novel allele
        // (see trio-dnm.11.vcf)
        //      chr1:10000057   f,m,c 11+00=11 .. LoH region
        //      chr1:10697377   f,m,c 11+11=01 .. usually these are indel ambiguities
        is_novel = ( (ca!=fa && ca!=fb && ca!=ma && ca!=mb) || (cb!=fa && cb!=fb && cb!=ma && cb!=mb) ) ? 1 : 0;
        if ( is_novel && *denovo_allele==0 ) is_novel = 0;  // never count reference allele as novel
    }
    else
    {
        is_novel = ( ((ca==fa||ca==fb) && (cb==ma||cb==mb)) || ((ca==ma||ca==mb) && (cb==fa||cb==fb)) ) ? 0 : 1;
    }

    if ( !is_novel )
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

    if ( fa!=fb )
    {
        // this must be a genotype error, father cannot be heterozygous in X; don't flag it as a DNM unless
        // also autosomal inheritance fails
        init_tprob_mprob(args,fi,mi,ci,tprob,mprob,denovo_allele);
    }
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
    if ( !args->dnm_score_tag )
    {
        if ( args->use_model==USE_NAIVE ) args->dnm_score_tag = strdup("DNM:flag");
        else args->dnm_score_tag = strdup("DNM:log");
    }
    char *ptr = strchr(args->dnm_score_tag,':');
    if ( ptr )
    {
        if ( ptr==args->dnm_score_tag ) error("Error: could not parse --dnm-tag %s\n",ptr);
        *ptr = 0;
        if ( !strcasecmp(ptr+1,"log") ) args->dnm_score_type = DNM_LOG;
        else if ( !strcasecmp(ptr+1,"phred") ) args->dnm_score_type = DNM_PHRED;
        else if ( !strcasecmp(ptr+1,"prob") ) args->dnm_score_type = DNM_PROB;
        else if ( !strcasecmp(ptr+1,"flag") ) args->dnm_score_type = DNM_FLAG;
        else error("Error: the type \"%s\" is not supported\n",ptr+1);
    }
    else
        args->dnm_score_type = DNM_LOG;
    if ( args->dnm_score_type==DNM_FLAG )
    {
        if ( !args->use_model ) args->use_model = USE_NAIVE;
        else if ( args->use_model!=USE_NAIVE ) error("The output type FLAG can be used only with --use-NAIVE\n");
    }
    if ( args->use_model==USE_NAIVE )
    {
        if ( !args->dnm_score_type ) args->dnm_score_type = DNM_FLAG;
        else if ( args->dnm_score_type!=DNM_FLAG ) error("The output type FLAG is required with --use-NAIVE\n");
    }
    if ( !args->use_model ) args->use_model = USE_ACM;

    args->sr = bcf_sr_init();
    if ( args->regions )
    {
        args->sr->require_index = 1;
        bcf_sr_set_opt(args->sr,BCF_SR_REGIONS_OVERLAP,args->regions_overlap);
        if ( bcf_sr_set_regions(args->sr, args->regions, args->regions_is_file)<0 ) error("Failed to read the regions: %s\n",args->regions);
    }
    if ( args->targets )
    {
        bcf_sr_set_opt(args->sr,BCF_SR_TARGETS_OVERLAP,args->targets_overlap);
        if ( bcf_sr_set_targets(args->sr, args->targets, args->targets_is_file, 0)<0 ) error("Failed to read the targets: %s\n",args->targets);
    }
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr = bcf_sr_get_header(args->sr,0);

    if ( args->filter_str ) args->filter = filter_init(args->hdr, args->filter_str);

    int id;
    if ( args->use_model==USE_NAIVE  )
    {
        if ( (id=bcf_hdr_id2int(args->hdr, BCF_DT_ID, "GT"))<0 || !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,id) )
            error("Error: the tag FORMAT/GT is not present in %s\n", args->fname);
    }
    else if ( (id=bcf_hdr_id2int(args->hdr, BCF_DT_ID, "PL"))<0 || !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,id) )
        error("Error: the tag FORMAT/PL is not present in %s\n", args->fname);
    args->need_QS = ( args->use_model==USE_ACM && !args->with_ppl && !args->with_pad ) ? 1 : 0;
    if ( args->need_QS && ((id=bcf_hdr_id2int(args->hdr, BCF_DT_ID, "QS"))<0 || !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,id)) )
        error(
            "Error:\n"
            "   The FORMAT/QS tag is not present. If you want to proceed anyway, add either `--with-pAD` or\n"
            "   `--with-pPL` option, the latter at the cost of inflated false discovery rate. The QS annotation\n"
            "    can be generated at the mpileup step together with the AD annotation using the command\n"
            "       bcftools mpileup -a AD,QS -f ref.fa file.bam\n");   // Possible future todo: use AD as a proxy for QS?
    if ( args->use_model!=USE_NAIVE )
    {
        if ( (id=bcf_hdr_id2int(args->hdr, BCF_DT_ID, "AD"))<0 || !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,id) )
            fprintf(stderr, "Warning: the tag FORMAT/AD is not present in %s, the output tag FORMAT/VAF will not be added\n", args->fname);
        else
            args->has_fmt_ad = 1;
        if ( args->with_pad && !args->has_fmt_ad )
            error("Error: no FORMAT/AD is present in %s, cannot run with --with-pAD\n", args->fname);
    }

    init_priors(args,&args->priors,autosomal);
    init_priors(args,&args->priors_X,chrX);
    init_priors(args,&args->priors_XX,chrXX);

    args->hdr_out = bcf_hdr_dup(args->hdr);
    char *type = NULL;
    if ( args->dnm_score_type==DNM_LOG ) type = "log scaled value (bigger value = bigger confidence)";
    if ( args->dnm_score_type==DNM_PHRED ) type = "phred value (bigger value = bigger confidence)";
    if ( args->dnm_score_type==DNM_PROB ) type = "probability";
    if ( args->dnm_score_type==DNM_FLAG ) type = "1 for Mendelian-incompatible genotypes";
    bcf_hdr_printf(args->hdr_out, "##FORMAT=<ID=%s,Number=1,Type=%s,Description=\"De-novo mutation score given as %s\">",args->dnm_score_tag,(args->dnm_score_type&DNM_INT)?"Integer":"Float",type);
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

    if ( args->record_cmd_line )
        bcf_hdr_append_version(args->hdr_out, args->argc, args->argv, "bcftools_trio-dnm2");

    char wmode[8];
    set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
    args->out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( bcf_hdr_write(args->out_fh, args->hdr_out)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
    if ( init_index2(args->out_fh,args->hdr_out,args->output_fname,
                     &args->index_fn, args->write_index)<0 )
        error("Error: failed to initialise index for %s\n",args->output_fname);

    if ( args->dnm_score_type & DNM_FLOAT )
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
    free(args->alt_tmp);
    free(args->alt_idx);
    free(args->idx);
    free(args->dnm_qual_int);
    free(args->dnm_qual_float);
    free(args->dnm_allele);
    free(args->vaf);
    free(args->trio);
    free(args->gt);
    free(args->pl);
    free(args->ad);
    free(args->qs);
    free(args->qs3);
    if ( args->write_index )
    {
        if ( bcf_idx_save(args->out_fh)<0 )
        {
            if ( hts_close(args->out_fh)!=0 ) error("Error: close failed .. %s\n", args->output_fname?args->output_fname:"stdout");
            error("Error: cannot write to index %s\n", args->index_fn);
        }
        free(args->index_fn);
    }
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
static inline double subtract_log(double a_log, double b_log)   // log(exp(a_log)-exp(b_log))
{
    return a_log + log(1 - exp(b_log - a_log));
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
    for (ca=0; ca<nals; ca++)               // ca,cb: iterate over all possible child's genotypes
    {
        for (cb=0; cb<=ca; cb++)
        {
            int cals = (1<<ca)|(1<<cb);
            double cpl = pl[iCHILD][ci];    // genotype likelihood of this allele combination, this is the FORMAT/PL
            int fi = 0;
            for (fa=0; fa<nals; fa++)       // fa,fb: father's genotypes
            {
                for (fb=0; fb<=fa; fb++)
                {
                    int fals = (1<<fa)|(1<<fb);

                    // Father and mother genotype likelihood can be either FORMAT/PL (with --with-pPL), FORMAT/QS, or
                    // fake QS calculated from FORMAT/AD (with --with-pAD)
                    double fpl;
                    if ( args->with_ppl ) fpl = pl[iFATHER][fi];
                    else
                    {
                        fpl = 0;
                        for (i=0; i<nals; i++)
                        {
                            // Is the i-th allele expected in the tested paternal fa|fb genotype? Strong presence
                            // of the allele will add little to fpl -> strong genotype fa|fb
                            if ( fals&(1<<i) )
                                fpl += subtract_log(0,qs[iFATHER][i]);

                            // Not part of the tested paternal genotype, but of the tested proband ca|cb genotype:
                            // penalize by increasing fpl -> weaker genotype fa|fb
                            else if ( cals&(1<<i) )
                                fpl += qs[iFATHER][i];

                            // Not part of the tested paternal or proband genotypes, and the tested paternal genotype
                            // is a hom: penalize presence of different allles
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
                            if ( args->with_ppl ) mpl = pl[iMOTHER][mi];
                            else
                            {
                                mpl = 0;
                                for (i=0; i<nals; i++)
                                {
                                    if ( mals&(1<<i) )
                                        mpl += subtract_log(0,qs[iMOTHER][i]);
                                    else if ( cals&(1<<i) )
                                        mpl += qs[iMOTHER][i];
                                    else if ( ma==mb )
                                        mpl += qs[iMOTHER][i];
                                }
                            }
                            double val = cpl + fpl + mpl + priors->pprob[fi][mi][ci];   // L_{p,f,m}
                            sum = sum_log(sum,val);                                     // this is the denominator, adding to \sum L_{p,f,m}
#define DEBUG 0
#if DEBUG
                            if(val!=-HUGE_VAL)
                                fprintf(stderr,"m,f,c: %d%d+%d%d=%d%d  dn=%d (%d,%d,%d)   mpl,fpl,cpl: %+e %+e %+e \t prior:%+e \t pval=%+e  sum=%+e  %c\n",
                                    mb,ma,fb,fa,cb,ca,priors->denovo[fi][mi][ci],fi,mi,ci,mpl,fpl,cpl,priors->pprob[fi][mi][ci], val,sum,(priors->denovo[fi][mi][ci] && max < val)?'*':'-');
#endif
                            // Is this a valid de novo combination of p,f,m genotypes (ie not inherited), and is it most likely thus far?
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

    if ( args->strictly_novel )
    {
        // Downplay de novo calls with alleles present in the parents
        int ial = *al1;
        sum = sum_log(sum,qs[iMOTHER][ial] + qs[iFATHER][ial]);
        max += qs[iMOTHER][ial] + qs[iFATHER][ial];
    }

    // This is the log( 1 - (\max L_pfm) / (\sum L_pfm) ). The default output (DNM:log) prints the inverse. Note log
    // values smaller than ~40 will be output as -inf due to a loss of precision. That's OK, we are interested
    // in values very close to 0, bigger than, say, -8
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
static int process_trio_naive(args_t *args, priors_t *priors, int nals, int32_t gts[3], int *denovo_allele)
{
    int fi = seq3[gts[iFATHER]];
    int mi = seq3[gts[iMOTHER]];
    int ci = seq3[gts[iCHILD]];
    assert( fi!=-1 && mi!=-1 && ci!=-1 );
    *denovo_allele = priors->denovo_allele[fi][mi][ci];
    return priors->denovo[fi][mi][ci];
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
static void many_alts_trim(args_t *args, int *_nals, double *pl[3], int *_npl, double *qs[3])
{
    assert(*_nals > 4);

    // Find the most likely set of alleles from FORMAT/QS
    int nals = *_nals;
    hts_expand(int,nals,args->malt_idx,args->alt_idx);
    hts_expand(double,nals,args->malt_tmp,args->alt_tmp);
    hts_expand(double,10,args->malt_tmp,args->alt_tmp);
    memset(args->alt_tmp,0,sizeof(*args->alt_tmp)*nals);
    int i,j,k;
    for (i=0; i<3; i++)
        for (j=0; j<nals; j++) args->alt_tmp[j] += qs[i][j];

    // sort in ascending order, make REF allele always come first; insertion sort
    double *arr = args->alt_tmp;
    int tmp, *idx = args->alt_idx;
    for (i=0; i<nals; i++) idx[i] = i;
    for (i=2; i<nals; i++)
        for (j=i; j>1 && arr[idx[j]] < arr[idx[j-1]]; j--)
            tmp = idx[j], idx[j] = idx[j-1], idx[j-1] = tmp;

    for (i=0; i<3; i++)
    {
        for (j=0; j<4; j++) args->alt_tmp[j] = qs[i][args->alt_idx[j]];
        memcpy(qs[i],args->alt_tmp,4*sizeof(*args->alt_tmp));
    }
    for (i=0; i<3; i++)
    {
        for (j=0; j<4; j++)
            for (k=0; k<=j; k++)
            {
                int idst = bcf_alleles2gt(j,k);
                int isrc = bcf_alleles2gt(args->alt_idx[j],args->alt_idx[k]);
                args->alt_tmp[idst] = pl[i][isrc];
            }
        memcpy(pl[i],args->alt_tmp,10*sizeof(*args->alt_tmp));
    }
    *_nals = 4;
    *_npl  = 10;
}
static void many_alts_translate(args_t *args, int *al0, int *al1)
{
    *al0 = args->alt_idx[*al0];
    *al1 = args->alt_idx[*al1];
}
static void set_trio_PL(args_t *args, trio_t *trio, double *ppl[3], int npl1)
{
    int j,k;
    for (j=0; j<3; j++) // j loops over iFATHER,iMOTHER,iCHILD
    {
        int32_t *src = args->pl + npl1 * trio->idx[j];
        double *dst = ppl[j] = args->pl3 + j*npl1;
        double sum = 0;
        for (k=0; k<npl1; k++) { dst[k] = phred2num(src[k]); sum += dst[k]; }
        for (k=0; k<npl1; k++) dst[k] = log(dst[k]/sum);
    }
}
static void set_trio_QS(args_t *args, trio_t *trio, double *pqs[3], int nqs1)
{
    int j,k;
    for (j=0; j<3; j++) // j loops over iFATHER,iMOTHER,iCHILD
    {
        int32_t *qs = args->qs + nqs1 * trio->idx[j];
        double *dst = pqs[j] = args->qs3 + j*nqs1;
        for (k=0; k<nqs1; k++) dst[k] = phred2log(qs[k]);
    }
}
static void set_trio_QS_noisy(args_t *args, trio_t *trio, double *pqs[3], int nqs1, int n_ad, pnoise_t *pnoise)
{
    int32_t *ad_f = NULL, *ad_m = NULL;     // AD in father and mother

    if ( n_ad && !pnoise->abs && !pnoise->abs1 && !pnoise->frac1 ) n_ad = 0;    // AD is not required
    if ( n_ad )
    {
        // Noise tolerance will be applied only for alleles observed in a single parent (a possible mosaic
        // variant in the parent) but not for alleles observed in both parents (likely an artefact)
        ad_f = args->ad + n_ad * trio->idx[iFATHER];
        ad_m = args->ad + n_ad * trio->idx[iMOTHER];
    }

    int j,k;
    for (j=0; j<3; j++)
    {
        int32_t *ad = n_ad ? args->ad + n_ad * trio->idx[j] : NULL;
        int32_t *qs = args->qs + nqs1 * trio->idx[j];
        pqs[j] = args->qs3 + j*nqs1;
        double pn = 0, pns = 0;
        double sum_qs = 0, sum_ad = 0;
        if ( (pnoise->frac || pnoise->frac1) && j!=iCHILD )
        {
            for (k=0; k<nqs1; k++) sum_qs += qs[k];
            pn  = sum_qs * pnoise->frac;
            pns = sum_qs * pnoise->frac1;
            if ( n_ad )
            {
                // "absolute" threshold: find the average QS per read from AD and use that
                // if bigger than the relative threshold
                for (k=0; k<n_ad; k++) sum_ad += ad[k];
                if ( pn  < pnoise->abs * sum_qs / sum_ad ) pn = pnoise->abs * sum_qs / sum_ad;
                if ( pns < pnoise->abs1 * sum_qs / sum_ad ) pns = pnoise->abs1 * sum_qs / sum_ad;
            }
        }
        // Reduce QS for all alleles to account for noise
        for (k=0; k<nqs1; k++)
        {
            double val = qs[k];
            if ( n_ad && (!ad_f[k] || !ad_m[k]) ) val -= pns;
            else val -= pn;
            if ( val < 0 ) val = 0;
            pqs[j][k] = phred2log(val);
        }
#if 0
        // The original code, don't like the capping at 255
        // Reduce QS for all alleles to account for noise
        for (k=0; k<nqs1; k++)
        {
            double val = qs[k];
            if ( !args->pnoise_strict || !ad_f[k] || !ad_m[k] ) val -= noise_tolerance;
            if ( val < 0 ) val = 0;
            if ( val > 255 ) val = 255;
            pqs[j][k] = phred2log(val);
        }
#endif
#if 0
        // None of this worked better
        double max_qs = 0;
        for (k=0; k<nqs1; k++)
        {
            double tmp = qs[k];
            if ( max_qs < tmp ) max_qs = tmp;
        }
        for (k=0; k<nqs1; k++)
        {
            // This uses absolute value of parental QS, does not regard the depth. The penalty is way too
            // strict for high-coverage sites, random parental sequencing errors there are very likely.
            // Puts trio-dnm.40.vcf at -24,-17, which is too strict.
            //      pqs[j][k] = phred2log(qs[k]);


            // Proportions in log space are not well motivated and does not work well.
            // Puts trio-dnm.40.vcf at -9,-13 which is too strict
            //      pqs[j][k] = phred2log(255*(qs[k]/max_qs));

            // This seems too lenient for 1/10 alts (-4.26326e-14) and feeds -inf,0,0 to everything.
            //      double tmp = max_qs - qs[k];
            //      pqs[j][k] = log(1 - exp(phred2log(tmp)));
        }
#endif
    }
}
static int set_trio_GT(args_t *args, trio_t *trio, int32_t gts[3], int ngts, int ignore_father)
{
    int j,k;
    for (j=0; j<3; j++)     // iFATHER,iMOTHER,iCHILD
    {
        int32_t *src = args->gt + ngts * trio->idx[j];
        for (k=0; k<ngts; k++)
        {
            int ial = src[k];
            if ( ial==bcf_int32_vector_end ) break;
            if ( bcf_gt_is_missing(ial) )
            {
                if ( j!=iFATHER || !ignore_father ) return -1;
                ial = 0;     // can be anything, will not be used
            }
            else
                ial = bcf_gt_allele(ial);
            gts[j] |= 1 << ial;
            assert(gts[j]>0 && gts[j]<13);
        }
        if ( !gts[j] && (j!=iFATHER || !ignore_father) ) return -1;
    }
    return 0;
}
static int set_trio_GT_many_alts(args_t *args, trio_t *trio, int32_t gts[3], int ngts, int nals, int ignore_father)
{
    int i,j,k, nused = 0;
    hts_expand(int,nals,args->malt_idx,args->alt_idx);
    for (i=0; i<nals; i++) args->alt_idx[i] = -1;
    for (j=0; j<3; j++)     // iFATHER,iMOTHER,iCHILD
    {
        int32_t *src = args->gt + ngts * trio->idx[j];
        for (k=0; k<ngts; k++)
        {
            int ial = src[k];
            if ( ial==bcf_int32_vector_end ) break;
            if ( bcf_gt_is_missing(ial) )
            {
                if ( j!=iFATHER || !ignore_father ) return -1;
                ial = 0;    // can be anything, will not be used
            }
            else
                ial = bcf_gt_allele(ial);
            if ( ial >= nals ) error("Error: FMT/GT contains incorrect allele \"%d\" at a site with %d alleles\n",ial,nals);
            if ( args->alt_idx[ial]==-1 )
            {
                args->alt_idx[ial] = nused++;
                if ( nused > 4 ) return -1;
            }
            gts[j] |= 1<<args->alt_idx[ial];
            assert(gts[j]>0 && gts[j]<13);
        }
        if ( !gts[j] ) return -1;
    }
    return 0;
}
static void process_record_naive(args_t *args, bcf1_t *rec)
{
    int nsmpl = bcf_hdr_nsamples(args->hdr);
    int ngts = bcf_get_genotypes(args->hdr,rec,&args->gt,&args->mgt);
    if ( ngts<=0 || ngts==nsmpl )
    {
        if ( bcf_write(args->out_fh, args->hdr_out, rec)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
        return;
    }
    ngts /= nsmpl;

    if ( ngts!=2 ) error("todo: ploidy>2\n");

    int is_chrX = 0;
    if ( regidx_overlap(args->chrX_idx,bcf_seqname(args->hdr,rec),rec->pos,rec->pos+rec->rlen,NULL) ) is_chrX = 1;

    int i, write_dnm = 0;
    for (i=0; i<nsmpl; i++) args->dnm_qual_int[i] = bcf_int32_missing;
    for (i=0; i<nsmpl; i++) args->dnm_allele[i] = bcf_int32_missing;
    for (i=0; i<args->ntrio; i++)
    {
        if ( args->filter && !args->trio[i].pass ) continue;

        int ignore_father = 0;  // father is irrelevant for male proband on chrX and can have missing GT
        priors_t *priors;
        if ( !is_chrX ) priors = &args->priors;
        else if ( args->trio[i].is_male ) priors = &args->priors_X, ignore_father = 1;
        else priors = &args->priors_XX;

        int32_t gts[3] = {0,0,0};
        int ret = rec->n_allele<=4 ? set_trio_GT(args,&args->trio[i],gts,ngts,ignore_father) : set_trio_GT_many_alts(args,&args->trio[i],gts,ngts,rec->n_allele,ignore_father);
        if ( ret<0 ) continue;

        int dnm_allele;
        double is_dnm = process_trio_naive(args, priors, rec->n_allele, gts, &dnm_allele);
        args->dnm_qual_int[ args->trio[i].idx[iCHILD] ] = is_dnm;
        args->dnm_allele[ args->trio[i].idx[iCHILD] ] = dnm_allele;
        if ( is_dnm ) write_dnm = 1;
    }
    if ( write_dnm )
    {
        int ret = bcf_update_format_int32(args->hdr_out,rec,args->dnm_score_tag,args->dnm_qual_int,nsmpl);
        if ( ret )
            error("Failed to write FORMAT/%s at %s:%"PRId64"\n", args->dnm_score_tag, bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        ret = bcf_update_format_int32(args->hdr_out,rec,args->dnm_allele_tag,args->dnm_allele,nsmpl);
        if ( ret )
            error("Failed to write FORMAT/%s at %s:%"PRId64"\n", args->dnm_allele_tag,bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
    }
    if ( bcf_write(args->out_fh, args->hdr_out, rec)!=0 )
        error("[%s] Error: cannot write to %s at %s:%"PRId64"\n", __func__,args->output_fname,bcf_seqname(args->hdr,rec),(int64_t)rec->pos+1);
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
    if ( args->use_model==USE_NAIVE )
    {
        process_record_naive(args, rec);
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

    int i,j, nqs1 = 0;
    if ( args->use_model==USE_ACM || rec->n_allele > 4 )    // DNG does not use QS, but QS is needed when trimming ALTs
    {
        nret = bcf_get_format_int32(args->hdr,rec,"QS",&args->qs,&args->mqs);
        if ( nret<0 )
        {
            if ( args->need_QS )
                error("Error: the FMT/QS tag is not available at %s:%"PRId64".\n",bcf_seqname(args->hdr,rec),(int64_t)rec->pos+1);
            if ( n_ad==0 )
            {
                static int missing_AD_warned = 0;
                if ( !missing_AD_warned )
                {
                    hts_log_warning(
                        "Neither FMT/QS nor FMT/AD present at %s:%"PRId64", cannot trim the number of alleles to four, skipping.\n"
                        "This warning is printed only once", bcf_seqname(args->hdr,rec),(int64_t)rec->pos+1);
                    missing_AD_warned = 1;
                }
                return;
            }

            // fake QS from AD assuming average BQ=30, used by --with-pAD
            nret = n_ad * nsmpl;
            hts_expand(int32_t,nret,args->mqs,args->qs);
            for (i=0; i<nret; i++)
            {
                if ( args->ad[i]==bcf_int32_missing || args->ad[i]==bcf_int32_vector_end ) args->qs[i] = args->ad[i];
                else args->qs[i] = 30 * args->ad[i];
            }
        }
        else if ( nret != nsmpl * rec->n_allele )
            error("Error: incorrect number of FMT/QS values at %s:%"PRId64".\n",bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);

        nqs1 = nret<=0 ? 0 : nret/nsmpl;
        hts_expand(double,3*nqs1,args->mqs3,args->qs3);
    }

    pnoise_t *pnoise = (bcf_get_variant_types(rec) & VCF_INDEL) ? &args->pn_indel : &args->pn_snv;

    int is_chrX = 0;
    if ( regidx_overlap(args->chrX_idx,bcf_seqname(args->hdr,rec),rec->pos,rec->pos+rec->rlen,NULL) ) is_chrX = 1;

    int al0, al1, write_dnm = 0, ad_set = 0;
    if ( args->dnm_score_type & DNM_FLOAT )
        for (i=0; i<nsmpl; i++) bcf_float_set_missing(args->dnm_qual_float[i]);
    else
        for (i=0; i<nsmpl; i++) args->dnm_qual_int[i] = bcf_int32_missing;
    for (i=0; i<nsmpl; i++) args->dnm_allele[i] = bcf_int32_missing;
    for (i=0; i<nsmpl; i++) args->vaf[i] = bcf_int32_missing;
    for (i=0; i<args->ntrio; i++)
    {
        if ( args->filter && !args->trio[i].pass ) continue;

        // Samples can be in any order in the VCF, set PL and QS to reflect the iFATHER,iMOTHER,iCHILD indices
        double *ppl[3];
        set_trio_PL(args,&args->trio[i],ppl,npl1);

        double *pqs[3];
        if ( args->use_model==USE_ACM )
            set_trio_QS_noisy(args,&args->trio[i],pqs,nqs1,n_ad,pnoise);
        else if ( rec->n_allele > 4 )       // DNG does not use QS, but QS is needed when trimming ALTs
            set_trio_QS(args,&args->trio[i],pqs,nqs1);

        priors_t *priors;
        if ( !is_chrX ) priors = &args->priors;
        else if ( args->trio[i].is_male ) priors = &args->priors_X;
        else priors = &args->priors_XX;

        int nals = rec->n_allele;
        if ( rec->n_allele > 4 ) many_alts_trim(args, &nals,ppl,&npl1,pqs);

        double score;
        if ( args->use_model==USE_ACM ) score = process_trio_ACM(args, priors, nals, ppl, npl1, pqs, &al0, &al1);
        else if ( args->use_model==USE_DNG ) score = process_trio_DNG(args, priors, nals, ppl, npl1, &al0, &al1);
        else error("Uh, this should not happen\n");

        if ( rec->n_allele > 4 ) many_alts_translate(args, &al0, &al1);

        if ( score >= args->min_score )
        {
            write_dnm = 1;
            if ( args->dnm_score_type==DNM_LOG )
                args->dnm_qual_float[ args->trio[i].idx[iCHILD] ] = score==HUGE_VAL ? 0 : subtract_log(0,phred2log(score));
            else if ( args->dnm_score_type==DNM_PROB )
                args->dnm_qual_float[ args->trio[i].idx[iCHILD] ] = score==HUGE_VAL ? 1 : 1 - phred2num(score);
            else
            {
                if ( score>255 ) score = 255;
                args->dnm_qual_int[ args->trio[i].idx[iCHILD] ] = round(score);
            }
            args->dnm_allele[ args->trio[i].idx[iCHILD] ] = al1;
        }

        if ( n_ad )
        {
            if ( al1 < n_ad )
            {
                ad_set = 1;
                for (j=0; j<3; j++)
                {
                    int32_t *src = args->ad + n_ad * args->trio[i].idx[j];
                    int k, ad_sum = 0;
                    for (k=0; k<rec->n_allele; k++) ad_sum += src[k];
                    args->vaf[ args->trio[i].idx[j] ] = ad_sum ? round(src[al1]*100./ad_sum) : 0;
                }
            }
            else
                for (j=0; j<3; j++) args->vaf[ args->trio[i].idx[j] ] = bcf_int32_missing;
        }
    }
    if ( write_dnm )
    {
        int ret;
        if ( args->dnm_score_type & DNM_FLOAT )
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
    static int warn_deprecated = 1;
    if  ( warn_deprecated )
    {
        fprintf(stderr,"Note: the `-u, --use` option will be deprecated in future releases.\n");
        warn_deprecated = 0;
    }
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
    else if ( !strcasecmp(opt,"pn") || !strcasecmp(opt,"pnoise") )
    {
        if ( !val ) error("Error: expected value with -u %s, e.g. -u %s=0.05\n",opt,opt);
        double pn_frac,pn_abs = 0;
        pn_frac = strtod(val,&tmp);
        if ( *tmp && *tmp==',' )
        {
            pn_abs = strtod(tmp+1,&tmp);
            if ( *tmp ) error("Could not parse: -u %s\n", optarg);
        }
        if ( pn_frac<0 || pn_frac>1 ) error("Error: expected value from the interval [0,1] for -u %s\n", optarg);
        if ( pn_abs<0 ) error("Error: expected positive value for -u %s\n", optarg);
        args->pn_snv.frac = pn_frac;
        args->pn_snv.abs  = pn_abs;
        args->pn_indel = args->pn_snv;
    }
    else if ( !strcasecmp(opt,"pns") )
    {
        if ( !val ) error("Error: expected value with -u %s, e.g. -u %s=0.05\n",opt,opt);
        double pn_frac,pn_abs = 0;
        pn_frac = strtod(val,&tmp);
        if ( *tmp && *tmp==',' )
        {
            pn_abs = strtod(tmp+1,&tmp);
            if ( *tmp ) error("Could not parse: -u %s\n", optarg);
        }
        if ( pn_frac<0 || pn_frac>1 ) error("Error: expected value from the interval [0,1] for -u %s\n", optarg);
        if ( pn_abs<0 ) error("Error: expected positive value for -u %s\n", optarg);
        args->pn_snv.frac1 = pn_frac;
        args->pn_snv.abs1  = pn_abs;
        args->pn_indel = args->pn_snv;
    }
    else if ( !strcasecmp(opt,"DNG") ) { args->use_model = USE_DNG; args->use_dng_priors = 1; }
    else if ( !strcasecmp(opt,"dng-priors") ) args->use_dng_priors = 1;
    else if ( !strcasecmp(opt,"ppl") ) args->with_ppl = 1;
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
    args->dnm_vaf_tag    = strdup("VAF");
    args->dnm_allele_tag = strdup("VA");
    args->mrate = 1e-8;
    args->pn_snv.frac  = 0.005;
    args->pn_snv.frac1 = 0.045;
    args->pn_snv.abs   = 0;
    args->pn_snv.abs1  = 0;
    args->record_cmd_line = 1;
    args->regions_overlap = 1;
    args->targets_overlap = 0;
    args->clevel = -1;
    static struct option loptions[] =
    {
        {"use",required_argument,0,'u'},
        {"force-AD",no_argument,0,1},
        {"dnm-tag",required_argument,0,2},
        {"va",required_argument,0,3},
        {"vaf",required_argument,0,4},
        {"dng-priors",no_argument,0,5},
        {"mrate",required_argument,0,6},
        {"pn",required_argument,0,7},
        {"pns",required_argument,0,8},
        {"use-DNG",no_argument,0,9},
        {"ppl",no_argument,0,10},
        {"with-pPL",no_argument,0,10},
        {"with-ppl",no_argument,0,10},
        {"use-NAIVE",no_argument,0,11},
        {"no-version",no_argument,NULL,12},
        {"with-pAD",no_argument,0,13},
        {"with-pad",no_argument,0,13},
        {"chrX",required_argument,0,'X'},
        {"min-score",required_argument,0,'m'},
        {"strictly-novel",no_argument,0,'n'},
        {"include",required_argument,0,'i'},
        {"exclude",required_argument,0,'e'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"ped",required_argument,NULL,'P'},
        {"pfm",required_argument,NULL,'p'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"regions-overlap",required_argument,NULL,14},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"targets-overlap",required_argument,NULL,15},
        {"write-index",optional_argument,NULL,'W'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    double pn_abs, pn_frac;
    while ((c = getopt_long(argc, argv, "p:P:o:O:s:i:e:r:R:t:T:m:au:X:nW::",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case  1 : args->force_ad = 1; break;
            case  2 : free(args->dnm_score_tag); args->dnm_score_tag = strdup(optarg); break;
            case  3 : free(args->dnm_allele_tag); args->dnm_allele_tag = strdup(optarg); break;
            case  4 : free(args->dnm_vaf_tag); args->dnm_vaf_tag = strdup(optarg); break;
            case  5 : args->use_dng_priors = 1; break;
            case  6 : args->mrate = strtod(optarg,&tmp); if ( *tmp ) error("Could not parse: --mrate %s\n", optarg); break;
            case  7 :
                pn_frac = strtod(optarg,&tmp);
                pn_abs  = 0;
                if ( *tmp && *tmp==',' )
                {
                    pn_abs = strtod(tmp+1,&tmp);
                    if ( *tmp )
                    {
                        if ( *tmp!=':' ) error("Could not parse: --pn %s\n", optarg);
                        if ( !strcasecmp("snv",tmp+1) ) { args->pn_snv.frac = pn_frac; args->pn_snv.abs = pn_abs; }
                        else if ( !strcasecmp("indel",tmp+1) ) { args->pn_indel.frac = pn_frac; args->pn_indel.abs = pn_abs; }
                        else error("Could not parse: --pn %s\n", optarg);
                    }
                    else
                    {
                        args->pn_snv.frac   = pn_frac;
                        args->pn_snv.abs    = pn_abs;
                        args->pn_indel.frac = pn_frac;
                        args->pn_indel.abs  = pn_abs;
                    }
                }
                else if ( *tmp ) error("Could not parse: --pn %s\n", optarg);
                else
                {
                    args->pn_snv.frac   = pn_frac;
                    args->pn_indel.frac = pn_frac;
                }
                if ( pn_frac<0 || pn_frac>1 ) error("Error: expected value from the interval [0,1] for --pn %s\n", optarg);
                if ( pn_abs<0 ) error("Error: expected positive value for --pn %s\n", optarg);
                break;
            case 8 :
                pn_frac = strtod(optarg,&tmp);
                pn_abs  = 0;
                if ( *tmp && *tmp==',' )
                {
                    pn_abs = strtod(tmp+1,&tmp);
                    if ( *tmp )
                    {
                        if ( *tmp!=':' ) error("Could not parse: --pns %s\n", optarg);
                        if ( !strcasecmp("snv",tmp+1) ) { args->pn_snv.frac1 = pn_frac; args->pn_snv.abs1 = pn_abs; }
                        else if ( !strcasecmp("indel",tmp+1) ) { args->pn_indel.frac1 = pn_frac; args->pn_indel.abs1 = pn_abs; }
                        else error("Could not parse: --pns %s\n", optarg);
                    }
                    else
                    {
                        args->pn_snv.frac1   = pn_frac;
                        args->pn_snv.abs1    = pn_abs;
                        args->pn_indel.frac1 = pn_frac;
                        args->pn_indel.abs1  = pn_abs;
                    }
                }
                else if ( *tmp ) error("Could not parse: --pns %s\n", optarg);
                else
                {
                    args->pn_snv.frac1   = pn_frac;
                    args->pn_indel.frac1 = pn_frac;
                }
                if ( pn_frac<0 || pn_frac>1 ) error("Error: expected value from the interval [0,1] for --pns %s\n", optarg);
                if ( pn_abs<0 ) error("Error: expected positive value for --pns %s\n", optarg);
                break;
            case 9  : args->use_model = USE_DNG; args->use_dng_priors = 1; break;
            case 10 : args->with_ppl = 1; break;
            case 11 : args->use_model = USE_NAIVE; break;
            case 12 : args->record_cmd_line = 0; break;
            case 13 : args->with_pad = 1; break;
            case 14 :
                args->regions_overlap = parse_overlap_option(optarg);
                if ( args->regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case 15 :
                args->targets_overlap = parse_overlap_option(optarg);
                if ( args->targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case 'W':
                if (!(args->write_index = write_index_parse(optarg)))
                    error("Unsupported index format '%s'\n", optarg);
                break;
            case 'X': args->chrX_list_str = optarg; break;
            case 'u': set_option(args,optarg); break;
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
                          default:
                          {
                              args->clevel = strtol(optarg,&tmp,10);
                              if ( *tmp || args->clevel<0 || args->clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                          }
                      };
                      if ( optarg[1] )
                      {
                          args->clevel = strtol(optarg+1,&tmp,10);
                          if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                      }
                      break;
            case 'P': args->ped_fname = optarg; break;
            case 'p': args->pfm = optarg; break;
            case 'm': args->min_score = strtod(optarg,&tmp);
                      if ( *tmp ) error("Could not parse: -M, --min-score %s\n", optarg);
                      break;
            case 'n': args->strictly_novel = 1; break;
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
