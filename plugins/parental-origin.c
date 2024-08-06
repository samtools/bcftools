/* The MIT License

   Copyright (c) 2019-2024 Genome Research Ltd.

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
#include <htslib/kfunc.h>
#include <errno.h>
#include "bcftools.h"
#include "filter.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define CNV_DEL 0
#define CNV_DUP 1

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
    int argc, filter_logic, cnv_type, debug, greedy;
    filter_t *filter;
    char *filter_str;
    char **argv, *pfm, *fname, *region;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr;
    trio_t trio;
    int32_t *pl, *ad, *gt;   // input FMT/PL, AD, and GT values
    int mpl, mad, mgt;
    double ppat,pmat;   // method 1: probability of paternal/maternal origin
    int ntest;          // number of informative sites
    int nmat, npat;     // method 2: number of pat/mat sites based on simple ad[0] < ad[1] comparison
    double min_pbinom;  // minimum binomial probability of paternal hets
}
args_t;

args_t args;

const char *about(void)
{
    return "Determine parental origin of a CNV region in a trio.\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: Determine parental origin of a CNV region\n"
        "Usage: bcftools +parental-origin [Plugin Options]\n"
        "Plugin options:\n"
        "   -b, --min-binom-prob FLOAT      exclude parental HETs with skewed ALT allele fraction [1e-2]\n"
        "   -d, --debug                     list informative sites\n"
        "   -e, --exclude EXPR              exclude sites and samples for which the expression is true\n"
        "   -g, --greedy                    use also ambiguous sites, e.g. het+hom parents for deletions\n"
        "   -i, --include EXPR              include sites and samples for which the expression is true\n"
        "   -p, --pfm P,F,M                 sample names of proband, father, and mother\n"
        "   -r, --region REGION             chr:beg-end\n"
        "   -t, --type <del|dup>            the CNV type\n"
        "\n"
        "Example:\n"
        "   bcftools +parental-origin -p proband,father,mother -t dup -r 14:22671179-22947951 file.bcf\n"
        "\n";
}

static void init_data(args_t *args)
{
    args->sr = bcf_sr_init();
    if ( args->region )
    {
        args->sr->require_index = 1;
        if ( bcf_sr_set_regions(args->sr, args->region, 0)<0 ) error("Failed to read the region: %s\n",args->region);
    }
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr = bcf_sr_get_header(args->sr,0);

    int id;
    if ( (id=bcf_hdr_id2int(args->hdr, BCF_DT_ID, "PL"))<0 || !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,id) )
        error("Error: the tag FORMAT/PL is not present in %s\n", args->fname);
    if ( (id=bcf_hdr_id2int(args->hdr, BCF_DT_ID, "AD"))<0 || !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,id) )
        error("Error: the tag FORMAT/AD is not present in %s\n", args->fname);
    if ( (id=bcf_hdr_id2int(args->hdr, BCF_DT_ID, "GT"))<0 || !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,id) )
        error("Error: the tag FORMAT/GT is not present in %s\n", args->fname);

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);

    int i, n = 0;
    char **list;
    list = hts_readlist(args->pfm, 0, &n);
    if ( n!=3 ) error("Expected three sample names with -t\n");
    args->trio.idx[iCHILD]  = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, list[0]);
    args->trio.idx[iFATHER] = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, list[1]);
    args->trio.idx[iMOTHER] = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, list[2]);
    for (i=0; i<n; i++)
    {
        if ( args->trio.idx[i] < 0 ) error("The sample is not present: %s\n", list[i]);
        free(list[i]);
    }
    free(list);
}
static void destroy_data(args_t *args)
{
    if ( args->filter ) filter_destroy(args->filter);
    free(args->pl);
    free(args->ad);
    free(args->gt);
    bcf_sr_destroy(args->sr);
    free(args);
}
static void process_record(args_t *args, bcf1_t *rec)
{
    if ( rec->n_allele!=2 || bcf_get_variant_types(rec)!=VCF_SNP ) return;

    int i,j;
    if ( args->filter )
    {
        uint8_t *smpl_pass = NULL;
        int pass_site = filter_test(args->filter, rec, (const uint8_t**) &smpl_pass);
        if ( args->filter_logic & FLT_EXCLUDE )
        {
            if ( pass_site )
            {
                if ( !smpl_pass ) return;
                pass_site = 0;
                for (i=0; i<3; i++)
                {
                    if ( smpl_pass[args->trio.idx[i]] ) smpl_pass[args->trio.idx[i]] = 0;
                    else { smpl_pass[args->trio.idx[i]] = 1; pass_site = 1; }
                }
                if ( !pass_site ) return;
            }
            else
                for (i=0; i<3; i++) smpl_pass[args->trio.idx[i]] = 1;
        }
        else if ( !pass_site ) return;

        if ( smpl_pass )
        {
            for (i=0; i<3; i++)
                if ( !smpl_pass[args->trio.idx[i]] ) return;
        }
    }

    int nsmpl = bcf_hdr_nsamples(args->hdr);
    int nret = bcf_get_format_int32(args->hdr,rec,"AD",&args->ad,&args->mad);
    if ( nret<=0 )
    {
        printf("The FORMAT/AD tag not present at %s:%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        return;
    }
    int nad1 = nret/nsmpl;

    nret = bcf_get_format_int32(args->hdr,rec,"PL",&args->pl,&args->mpl);
    if ( nret<=0 ) error("The FORMAT/PL tag not present at %s:%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
    int npl1 = nret/nsmpl;
    if ( npl1!=rec->n_allele*(rec->n_allele+1)/2 )
    {
        printf("todo: not a diploid site at %s:%"PRId64": %d alleles, %d PLs\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1,rec->n_allele,npl1);
        return;
    }

    nret = bcf_get_genotypes(args->hdr,rec,&args->gt,&args->mgt);
    if ( nret<=0 ) error("The FORMAT/GT tag not present at %s:%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
    int ngt1 = nret/nsmpl;
    if ( ngt1!=2 ) error("Todo: assuming diploid fields for now .. %s:%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);

    // number of ref and alt alleles in the proband
    int32_t ad[6], *adP = ad, *adF = ad+2, *adM = ad+4;
    int32_t dsg[3], *dsgP = dsg, *dsgF = dsg+1, *dsgM = dsg+2;
    double gl[9], *glP = gl, *glF = gl+3, *glM = gl+6;
    for (i=0; i<3; i++)    // trio
    {
        int isum = 0;
        int32_t *src = args->pl + npl1*args->trio.idx[i];
        double *gl_dst  = gl + 3*i;
        double sum  = 0;
        for (j=0; j<3; j++)     // iterate over PL
        {
            if ( src[j]==bcf_int32_missing || src[j]==bcf_int32_vector_end ) return;
            gl_dst[j] = pow(10,-0.1*src[j]);
            sum += gl_dst[j];
            isum += src[j];
        }
        if ( isum==0 ) return;
        for (j=0; j<3; j++) gl_dst[j] /= sum;

        int32_t *gt = args->gt + ngt1*args->trio.idx[i];
        dsg[i] = 0;
        for (j=0; j<ngt1; j++)
        {
            if ( gt[j]==bcf_int32_vector_end ) return;
            if ( bcf_gt_is_missing(gt[j]) ) return;
            if ( bcf_gt_allele(gt[j]) ) dsg[i]++;
        }

        src = args->ad + nad1*args->trio.idx[i];
        ad[2*i]   = src[0];
        ad[2*i+1] = src[1];
    }

    #define is_RR(x) (x[0]==0)
    #define is_RA(x) (x[1]==0)
    #define is_AA(x) (x[2]==0)
    if ( args->cnv_type==CNV_DEL )
    {
        if ( *dsgP!=0 && *dsgP!=2 ) return;     // proband not a hom
        if ( *dsgF == *dsgM ) return;           // cannot distinguish between parents
        if ( !args->greedy )
        {
            if ( *dsgF==1 && *dsgP==*dsgM ) return; // both parents have the proband's allele
            if ( *dsgM==1 && *dsgP==*dsgF ) return;
        }
        double pmat = glP[0]*(0.5*glM[0]*glF[0] + 2/3.*glM[0]*glF[1] + glM[0]*glF[2] + 1/3.*glM[1]*glF[0] + 0.5*glM[1]*glF[1] + glM[1]*glF[2]) +
                      glP[2]*(0.5*glM[2]*glF[2] + 2/3.*glM[2]*glF[1] + glM[2]*glF[0] + 1/3.*glM[1]*glF[2] + 0.5*glM[1]*glF[1] + glM[1]*glF[0]);
        double ppat = glP[0]*(0.5*glM[0]*glF[0] + 2/3.*glM[1]*glF[0] + glM[2]*glF[0] + 1/3.*glM[0]*glF[1] + 0.5*glM[1]*glF[1] + glM[2]*glF[1]) +
                      glP[2]*(0.5*glM[2]*glF[2] + 2/3.*glM[1]*glF[2] + glM[0]*glF[2] + 1/3.*glM[2]*glF[1] + 0.5*glM[1]*glF[1] + glM[0]*glF[1]);

        // NB: pmat/ppat is the probability of parental origin of the observed, not the deleted allele;
        //     args->pmat/ppat is the probability of parental origin of the deleted allele
        args->pmat += log(ppat);
        args->ppat += log(pmat);
        args->ntest++;

        if ( args->debug )
        {
            // output: position, paternal probability, maternal probability, PLs of child, father, mother
            printf("DBG\t%"PRId64"\t%e\t%e\t", (int64_t) rec->pos+1,ppat,pmat);
            for (i=0; i<3; i++)
            {
                for (j=0; j<3; j++)  printf(" %d",args->pl[npl1*args->trio.idx[i]+j]);
                printf("\t");
            }
            printf("\n");
        }
    }
    if ( args->cnv_type==CNV_DUP )
    {
        if ( !adP[0] || !adP[1] ) return;   // proband is homozygous or has no coverage
        if ( adP[0] == adP[1] ) return;     // proband's alleles are not informative, any or none could have been duplicated
        if ( *dsgP!=1 ) return;             // the proband's genotype is not a het
        if ( *dsgF == *dsgM ) return;       // cannot distinguish between parents

        if ( args->min_pbinom!=0 )
        {
            // exclude parental hets with skewed ALT allele proportion
            if ( *dsgF==1 && adF[0] && adF[1] && calc_binom_two_sided(adF[0],adF[1],0.5) < args->min_pbinom ) return;
            if ( *dsgM==1 && adM[0] && adM[1] && calc_binom_two_sided(adM[0],adM[1],0.5) < args->min_pbinom ) return;
        }

        double prra = glP[1] * calc_binom_one_sided(adP[1],adP[0],1/3.,1);
        double praa = glP[1] * calc_binom_one_sided(adP[1],adP[0],2/3.,0);
        double ppat = prra*(glM[1]*glF[0] + glM[2]*glF[0] + 0.5*glM[1]*glF[1] + glM[2]*glF[1]) +
                      praa*(glM[1]*glF[2] + glM[0]*glF[2] + 0.5*glM[1]*glF[1] + glM[0]*glF[1]);
        double pmat = prra*(glM[0]*glF[1] + glM[0]*glF[2] + 0.5*glM[1]*glF[1] + glM[1]*glF[2]) +
                      praa*(glM[2]*glF[1] + glM[2]*glF[0] + 0.5*glM[1]*glF[1] + glM[1]*glF[0]);
        args->pmat += log(pmat);
        args->ppat += log(ppat);
        args->ntest++;

        if ( args->debug )
        {
            // output: position; paternal probability; maternal probability; ADs of child, father,mother; PLs of child, father, mother
            printf("DBG\t%"PRId64"\t%e\t%e\t", (int64_t) rec->pos+1,ppat,pmat);
            for (i=0; i<3; i++)
            {
                printf("%d %d\t",ad[2*i],ad[2*i+1]);
            }
            for (i=0; i<3; i++)
            {
                for (j=0; j<3; j++)  printf(" %d",args->pl[npl1*args->trio.idx[i]+j]);
                printf("\t");
            }
            printf("\n");
        }
    }
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->min_pbinom = 1e-2;
    static struct option loptions[] =
    {
        {"include",required_argument,0,'i'},
        {"exclude",required_argument,0,'e'},
        {"pfm",required_argument,NULL,'p'},
        {"region",required_argument,0,'r'},
        {"type",required_argument,0,'t'},
        {"debug",no_argument,0,'d'},
        {"greedy",no_argument,0,'g'},
        {"min-binom-prob",required_argument,0,'b'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "h?e:i:p:r:t:dgb:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 't':
                if ( !strcasecmp("dup",optarg) ) args->cnv_type = CNV_DUP;
                else if ( !strcasecmp("del",optarg) ) args->cnv_type = CNV_DEL;
                break;
            case 'r': args->region = optarg; break;
            case 'p': args->pfm = optarg; break;
            case 'd': args->debug = 1; break;
            case 'g': args->greedy = 1; break;
            case 'b':
                args->min_pbinom = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -b %s\n", optarg);
                if ( args->min_pbinom<0 || args->min_pbinom>1 ) error("Expected value from the interval [0,1] with --min-binom-prob\n");
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

    if ( !args->pfm ) error("Missing the -p option\n");

    init_data(args);
    if ( args->debug )
    {
        if ( args->cnv_type==CNV_DEL ) printf("# DBG: position; paternal probability; maternal probability; PLs of child, father, mother\n");
        else printf("# DBG: position; paternal probability; maternal probability; ADs of child, father, mother; PLs of child, father, mother\n");
    }

    while ( bcf_sr_next_line(args->sr) )
        process_record(args, bcf_sr_get_line(args->sr,0));

    double qual = 4.3429*fabs(args->ppat - args->pmat);
    char *origin = "uncertain";
    if ( args->ppat > args->pmat ) origin = "paternal";
    else if ( args->ppat < args->pmat ) origin = "maternal";

    int i;
    printf("# bcftools +%s", args->argv[0]);
    for (i=1; i<args->argc; i++) printf(" %s",args->argv[i]);
    printf("\n");
    printf("# [1]type\t[2]predicted_origin\t[3]quality\t[4]nmarkers\n");
    printf("%s\t%s\t%f\t%d\n", args->cnv_type==CNV_DUP ? "dup" : "del", origin, qual, args->ntest);

    destroy_data(args);

    return 0;
}
