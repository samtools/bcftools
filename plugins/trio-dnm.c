/* The MIT License

   Copyright (c) 2018-2019 Genome Research Ltd.

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
#include <getopt.h>
#include <math.h>
#include <unistd.h>     // for isatty
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <errno.h>
#include "bcftools.h"
#include "filter.h"


// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define iCHILD  0
#define iFATHER 1
#define iMOTHER 2

typedef struct
{
    int idx[3];     // VCF sample index for child, father, mother
    int pass;       // do all three pass the filters?
}
trio_t;

typedef struct
{
    int argc, filter_logic, regions_is_file, targets_is_file, output_type;
    char *filter_str;
    char **argv, *ped_fname, *pfm, *output_fname, *fname, *regions, *targets;
    htsFile *out_fh;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr, *hdr_out;
    trio_t *trio;
    int has_fmt_ad;
    int ntrio, mtrio;
    int32_t *pl, *ad, *dnm_qual, *vaf;    // input FMT/PL and AD values, output DNM and VAF
    int mpl, mad;
    double min_score;
    double *aprob;  // proband's allele probabilities
    double *pl3;    // normalized PLs converted to probs for proband,father,mother
    int maprob, mpl3, midx, *idx, force_ad;
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
        "Usage: bcftools +trio-dnm [Plugin Options]\n"
        "Plugin options:\n"
        "   -e, --exclude EXPR              exclude sites and samples for which the expression is true\n"
        "       --force-AD                  calculate VAF even if the number of FMT/AD fields is incorrect. Use at your own risk!\n"
        "   -i, --include EXPR              include sites and samples for which the expression is true\n"
        "   -m, --min-score NUM             do not add FMT/DNM annotation if the score is smaller than NUM\n"
        "   -o, --output FILE               output file name [stdout]\n"
        "   -O, --output-type <b|u|z|v>     b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n"
        "   -p, --pfm P,F,M                 sample names of proband, father, and mother\n"
        "   -P, --ped FILE                  PED file\n"
        "   -r, --regions REG               restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE         restrict to regions listed in a file\n"
        "   -t, --targets REG               similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE         similar to -R but streams rather than index-jumps\n"
        "\n"
        "Example:\n"
        "   # Annotate VCF with FORMAT/DNM, run for a single trio\n"
        "   bcftools +trio-dnm -p proband,father,mother file.bcf\n"
        "\n"
        "   # Same as above, but read the trio(s) from a PED file\n"
        "   bcftools +trio-dnm -P file.ped file.bcf\n"
        "\n"
        "   # Same as above plus extract a list of significant DNMs using the bcftools/query command\n"
        "   bcftools +trio-dnm -P file.ped file.bcf -Ou | bcftools query -i'DNM>10' -f'[%CHROM:%POS %SAMPLE %DNM\\n]'\n"
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

        args->ntrio++;
        hts_expand0(trio_t,args->ntrio,args->mtrio,args->trio);
        trio_t *trio = &args->trio[args->ntrio-1];
        trio->idx[iFATHER] = father;
        trio->idx[iMOTHER] = mother;
        trio->idx[iCHILD]  = child;
    }
    while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );

    fprintf(stderr,"Identified %d complete trio%s in the VCF file\n", args->ntrio,args->ntrio==1?"":"s");

    // sort the sample by index so that they are accessed more or less sequentially
    qsort(args->trio,args->ntrio,sizeof(trio_t),cmp_trios);

    free(str.s);
    free(off);
    if ( hts_close(fp)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,fname);
}
static void init_data(args_t *args)
{
    args->sr = bcf_sr_init();
    if ( args->regions )
    {
        args->sr->require_index = 1;
        if ( bcf_sr_set_regions(args->sr, args->regions, args->regions_is_file)<0 ) error("Failed to read the regions: %s\n",args->regions);
    }
    if ( args->targets && bcf_sr_set_targets(args->sr, args->targets, args->targets_is_file, 0)<0 ) error("Failed to read the targets: %s\n",args->targets);
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr = bcf_sr_get_header(args->sr,0);

    int id;
    if ( (id=bcf_hdr_id2int(args->hdr, BCF_DT_ID, "PL"))<0 || !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,id) )
        error("Error: the tag FORMAT/PL is not present in %s\n", args->fname);
    if ( (id=bcf_hdr_id2int(args->hdr, BCF_DT_ID, "AD"))<0 || !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,id) )
        fprintf(stderr, "Warning: the tag FORMAT/AD is not present in %s, the output tag FORMAT/VAF will not be added\n", args->fname);
    else
        args->has_fmt_ad = 1;

    args->hdr_out = bcf_hdr_dup(args->hdr);
    bcf_hdr_append(args->hdr_out, "##FORMAT=<ID=DNM,Number=1,Type=Integer,Description=\"De-novo mutation score, bigger values = bigger confidence\">");
    if ( args->has_fmt_ad )
        bcf_hdr_append(args->hdr_out, "##FORMAT=<ID=VAF,Number=1,Type=Integer,Description=\"The percentage of ALT reads\">");

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

    args->out_fh = hts_open(args->output_fname,hts_bcf_wmode(args->output_type));
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( bcf_hdr_write(args->out_fh, args->hdr_out)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);

    args->dnm_qual = (int32_t*) malloc(sizeof(*args->dnm_qual)*bcf_hdr_nsamples(args->hdr));
    args->vaf      = (int32_t*) malloc(sizeof(*args->vaf)*bcf_hdr_nsamples(args->hdr));
}
static void destroy_data(args_t *args)
{
    free(args->pl3);
    free(args->aprob);
    free(args->idx);
    free(args->dnm_qual);
    free(args->vaf);
    free(args->trio);
    free(args->pl);
    free(args->ad);
    if ( hts_close(args->out_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname);
    bcf_hdr_destroy(args->hdr_out);
    bcf_sr_destroy(args->sr);
    free(args);
}
static float process_trio(args_t *args, int nals, double *pl[3], int npl, int *al0, int *al1)
{
    assert( nals>1 );

    // determine the two most likely proband's alleles
    int i,j,k = 0,tmp;

    hts_expand(int,nals,args->midx,args->idx);
    hts_expand(double,nals,args->maprob,args->aprob);
    for (i=0; i<nals; i++) args->aprob[i] = 0;
    for (i=0; i<nals; i++)
    {
        for (j=0; j<=i; j++)
        {
            args->aprob[i] += pl[iCHILD][k];
            args->aprob[j] += pl[iCHILD][k];
            k++;
        }
    }

    // sort in descendent order
    double *arr = args->aprob;
    int *idx = args->idx;
    for (i=0; i<nals; i++) idx[i] = i;
    for (i=1; i<nals; i++)
        for (j=i; j>0 && arr[idx[j]] > arr[idx[j-1]]; j--)
            tmp = idx[j], idx[j] = idx[j-1], idx[j-1] = tmp;

    if ( idx[0] < idx[1] ) { *al0 = idx[0]; *al1 = idx[1]; }
    else { *al0 = idx[1]; *al1 = idx[0]; }

    // Calculate the probability of inheriting the 00, 01, and 11 genotype. For DNM they all will be small
    int k00 = bcf_alleles2gt(idx[0],idx[0]);
    int k01 = bcf_alleles2gt(idx[0],idx[1]);
    int k11 = bcf_alleles2gt(idx[1],idx[1]);
    double pd00 = pl[iCHILD][k00] * (pl[iFATHER][k00] + 0.5*pl[iFATHER][k01]) * (pl[iMOTHER][k00] + 0.5*pl[iMOTHER][k01]);
    double pd11 = pl[iCHILD][k11] * (pl[iFATHER][k11] + 0.5*pl[iFATHER][k01]) * (pl[iMOTHER][k11] + 0.5*pl[iMOTHER][k01]);
    double pd01 = pl[iCHILD][k01] * (pl[iFATHER][k00] * (pl[iMOTHER][k11] + 0.5*pl[iMOTHER][k01]) + pl[iFATHER][k11] * (pl[iMOTHER][k00] + 0.5*pl[iMOTHER][k01])
        + 0.5*pl[iFATHER][k01] * (pl[iMOTHER][k00] + pl[iMOTHER][k01] + pl[iMOTHER][k11]));

    double max = pd01;
    if ( max < pd00 ) max = pd00;
    if ( max < pd11 ) max = pd11;
    return fabs(4.3429 * log(max));
}
static void process_record(args_t *args, bcf1_t *rec)
{
    if ( rec->n_allele==1 )
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
        error("fixme: not a diploid site at %s:%"PRId64": %d alleles, %d PLs\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1,rec->n_allele,npl1);
    hts_expand(double,3*npl1,args->mpl3,args->pl3);
    int i, j, k, al0, al1, write_dnm = 0, ad_set = 0;
    for (i=0; i<nsmpl; i++) args->dnm_qual[i] = bcf_int32_missing;
    for (i=0; i<args->ntrio; i++)
    {
        double *ppl[3];
        for (j=0; j<3; j++)
        {
            int32_t *src = args->pl + npl1 * args->trio[i].idx[j];
            double *dst = ppl[j] = args->pl3 + j*npl1;
            double sum = 0;
            for (k=0; k<npl1; k++) { dst[k] = pow(10,-0.1*src[k]); sum += dst[k]; }
            for (k=0; k<npl1; k++) dst[k] /= sum;
        }
        int32_t score = process_trio(args, rec->n_allele, ppl, npl1, &al0, &al1);
        if ( score >= args->min_score )
        {
            write_dnm = 1;
            args->dnm_qual[ args->trio[i].idx[iCHILD] ] = score;
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
        if ( bcf_update_format_int32(args->hdr_out,rec,"DNM",args->dnm_qual,nsmpl)!=0 )
            error("Failed to write FORMAT/DNM at %s:%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        if ( ad_set )
        {
            if ( bcf_update_format_int32(args->hdr_out,rec,"VAF",args->vaf,nsmpl)!=0 )
                error("Failed to write FORMAT/VAF at %s:%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        }
    }
    if ( bcf_write(args->out_fh, args->hdr_out, rec)!=0 ) error("[%s] Error: cannot write to %s at %s:%"PRId64"\n", __func__,args->output_fname,bcf_seqname(args->hdr,rec),(int64_t)rec->pos+1);
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_fname = "-";
    static struct option loptions[] =
    {
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
    while ((c = getopt_long(argc, argv, "p:P:o:O:s:i:e:r:R:t:T:m:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case  1 : args->force_ad = 1; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
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
                      if ( *tmp ) error("Could not parse: --min-score %s\n", optarg);
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
