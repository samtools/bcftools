/* The MIT License

   Copyright (c) 2016-2020 Genome Research Ltd.

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
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/kfunc.h>
#include <htslib/kbitset.h>
#include <htslib/vcfutils.h>
#include <inttypes.h>
#include "bcftools.h"
#include "convert.h"

typedef struct
{
    int smpl,ctrl;      // VCF sample index
    const char *smpl_name, *ctrl_name;
}
pair_t;

typedef struct
{
    bcf_hdr_t *hdr;
    pair_t *pair;
    int npair, mpair, min_dp, min_alt_dp;
    int32_t *ad_arr;
    int mad_arr;
    double th;
    convert_t *convert;
    kstring_t str;
    uint64_t nsite,ncmp;
    int variant_type;
    int clean_vcf;
    kbitset_t *rm_als;
}
args_t;

args_t args;

const char *about(void)
{
    return "Find positions with wildly varying ALT allele frequency (Fisher test on FMT/AD).\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Find positions with wildly varying ALT allele frequency (Fisher test on FMT/AD).\n"
        "Usage: bcftools +ad-bias [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -a, --min-alt-dp <int>          Minimum required alternate allele depth [1]\n"
        "   -c, --clean-vcf                 Outputs VCF removing sites and ALT alleles not passing the -t threshold\n"
        "   -d, --min-dp <int>              Minimum required depth [0]\n"
        "   -f, --format <string>           Optional tags to append to output (`bcftools query` style of format)\n"
        "   -s, --samples <file>            List of sample pairs, one tab-delimited pair per line\n"
        "   -t, --threshold <float>         Output only hits with p-value smaller than <float> [1e-3]\n"
        "   -v, --variant-type <snp|indel>  Consider only variants of this type. (By default all variants are considered.)\n"
        "\n"
        "Example:\n"
        "   bcftools +ad-bias file.bcf -- -t 1e-3 -s samples.txt\n"
        "\n";
}

void parse_samples(args_t *args, char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);

    kstring_t str = {0,0,0};
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Empty file: %s\n", fname);

    int moff = 0, *off = NULL;
    do
    {
        // HPSI0513i-veqz_6    HPSI0513pf-veqz
        int ncols = ksplit_core(str.s,'\t',&moff,&off);
        if ( ncols<2 ) error("Could not parse the sample file: %s\n", str.s);

        int smpl = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[0]]);
        if ( smpl<0 ) continue;
        int ctrl = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[1]]);
        if ( ctrl<0 ) continue;

        args->npair++;
        hts_expand0(pair_t,args->npair,args->mpair,args->pair);
        pair_t *pair = &args->pair[args->npair-1];
        pair->ctrl = ctrl;
        pair->smpl = smpl;
        pair->smpl_name = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,pair->smpl);
        pair->ctrl_name = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,pair->ctrl);
    } while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );

    free(str.s);
    free(off);
    if ( hts_close(fp)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,fname);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    memset(&args,0,sizeof(args_t));
    args.hdr = in;
    args.th  = 1e-3;
    args.min_alt_dp = 1;
    char *fname = NULL, *format = NULL;
    static struct option loptions[] =
    {
        {"clean-vcf",required_argument,NULL,'c'},
        {"min-dp",required_argument,NULL,'d'},
        {"min-alt-dp",required_argument,NULL,'a'},
        {"format",required_argument,NULL,'f'},
        {"samples",required_argument,NULL,'s'},
        {"threshold",required_argument,NULL,'t'},
        {"variant-type",required_argument,NULL,'v'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "?hs:t:f:d:a:v:c",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'c': args.clean_vcf = 1; break;
            case 'a':
                args.min_alt_dp = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse: -a %s\n", optarg);
                break;
            case 'd':
                args.min_dp = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse: -d %s\n", optarg);
                break;
            case 't':
                args.th = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -t %s\n", optarg);
                break;
            case 's': fname = optarg; break;
            case 'v': 
                if ( !strcasecmp(optarg,"snp") || !strcasecmp(optarg,"snps") ) args.variant_type = VCF_SNP;
                else if ( !strcasecmp(optarg,"indel") || !strcasecmp(optarg,"indels") ) args.variant_type = VCF_INDEL;
                else error("Error: Variant type \"%s\" is not supported\n",optarg);
                break;
            case 'f': format = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( !fname ) error("Expected the -s option\n");
    parse_samples(&args, fname);
    if ( format )
    {
        if ( args.clean_vcf ) error("The option -f cannot be used together with -c\n");
        args.convert = convert_init(args.hdr, NULL, 0, format);
    }
    if ( args.clean_vcf ) return 0;

    printf("# This file was produced by: bcftools +ad-bias(%s+htslib-%s)\n", bcftools_version(),hts_version());
    printf("# The command line was:\tbcftools +ad-bias %s", argv[0]);
    for (c=1; c<argc; c++) printf(" %s",argv[c]);
    printf("\n#\n");

    int i = 1;
    printf("# FT, Fisher Test");
    printf("\t[%d]Sample", ++i);
    printf("\t[%d]Control", ++i);
    printf("\t[%d]Chrom", ++i);
    printf("\t[%d]Pos", ++i);
    printf("\t[%d]REF", ++i);
    printf("\t[%d]ALT", ++i);
    printf("\t[%d]smpl.nREF", ++i);
    printf("\t[%d]smpl.nALT", ++i);
    printf("\t[%d]ctrl.nREF", ++i);
    printf("\t[%d]ctrl.nALT", ++i);
    printf("\t[%d]P-value", ++i);
    if ( format ) printf("\t[%d-]User data: %s", ++i, format);
    printf("\n");
    return 1;
}

bcf1_t *process(bcf1_t *rec)
{
    if ( rec->n_allele < 2 ) return NULL;

    int nad = bcf_get_format_int32(args.hdr, rec, "AD", &args.ad_arr, &args.mad_arr);
    if ( nad<0 ) return NULL;
    nad /= bcf_hdr_nsamples(args.hdr);
    
    if ( args.convert ) convert_line(args.convert, rec, &args.str);
    args.nsite++;

    int keep_als = 0;
    if ( args.clean_vcf )
    {
        if ( !args.rm_als )
            args.rm_als = kbs_init(rec->n_allele);
        else if ( args.rm_als->n_max < rec->n_allele )
            kbs_resize(&args.rm_als, rec->n_allele);
        kbs_insert_all(args.rm_als);
    }

    int i,j;
    for (i=0; i<args.npair; i++)
    {
        pair_t *pair = &args.pair[i];
        int32_t *aptr = args.ad_arr + nad*pair->smpl;
        int32_t *bptr = args.ad_arr + nad*pair->ctrl;

        // Find the two most frequent alleles
        int nbig=-1,nsmall=-1,ibig=-1,ismall=-1;
        for (j=0; j<nad; j++)
        {
            if ( aptr[j]==bcf_int32_missing ) continue;
            if ( aptr[j]==bcf_int32_vector_end ) break;
            if ( ibig==-1 ) { ibig = j, nbig = aptr[j]; continue; }
            if ( nbig < aptr[j] )
            {
                if ( ismall==-1 || nsmall < nbig ) ismall = ibig, nsmall = nbig;
                ibig = j, nbig = aptr[j];
                continue;
            }
            if ( ismall==-1 || nsmall < aptr[j] ) ismall = j, nsmall = aptr[j];
        }
        for (j=0; j<nad; j++)
        {
            if ( bptr[j]==bcf_int32_missing ) continue;
            if ( bptr[j]==bcf_int32_vector_end ) break;
            if ( ibig==-1 ) { ibig = j, nbig = bptr[j]; continue; }
            if ( ibig==j )
            {
                if ( nbig < bptr[j] ) nbig = bptr[j];
                continue;
            }
            if ( nbig < bptr[j] )
            {
                if ( ismall==-1 || nsmall < nbig ) ismall = ibig, nsmall = nbig;
                ibig = j, nbig = bptr[j];
                continue;
            }
            if ( ismall==-1 || nsmall < bptr[j] ) ismall = j, nsmall = bptr[j];
        }
        if ( ibig==-1 || ismall==-1 ) continue;         // only one non-missing allele
        if ( nbig + nsmall < args.min_dp ) continue;    // low depth

        if ( aptr[ibig]==bcf_int32_missing || aptr[ibig]==bcf_int32_vector_end ) continue;
        if ( bptr[ibig]==bcf_int32_missing || bptr[ibig]==bcf_int32_vector_end ) continue;
        if ( aptr[ismall]==bcf_int32_missing || aptr[ismall]==bcf_int32_vector_end ) continue;
        if ( bptr[ismall]==bcf_int32_missing || bptr[ismall]==bcf_int32_vector_end ) continue;

        if ( args.variant_type )
        {
            if ( args.variant_type==VCF_SNP && strlen(rec->d.allele[ibig])!=strlen(rec->d.allele[ismall]) ) continue;
            if ( args.variant_type==VCF_INDEL && strlen(rec->d.allele[ibig])==strlen(rec->d.allele[ismall]) ) continue;
        }

        int iref,ialt,nalt;
        if ( ibig > ismall ) ialt = ibig, iref = ismall, nalt = nbig; 
        else ialt = ismall, iref = ibig, nalt = nsmall;

        if ( nalt < args.min_alt_dp ) continue;

        args.ncmp++;

        int n11 = aptr[iref], n12 = aptr[ialt];
        int n21 = bptr[iref], n22 = bptr[ialt];
        double left, right, fisher;
        kt_fisher_exact(n11,n12,n21,n22, &left,&right,&fisher);
        if ( fisher >= args.th ) continue;

        if ( args.clean_vcf )
        {
            keep_als = 1;
            kbs_delete(args.rm_als, ialt);
            continue;
        }

        printf("FT\t%s\t%s\t%s\t%"PRId64"\t%s\t%s\t%d\t%d\t%d\t%d\t%e",
            pair->smpl_name,pair->ctrl_name,
            bcf_hdr_id2name(args.hdr,rec->rid), (int64_t) rec->pos+1,
            rec->d.allele[iref],rec->d.allele[ialt],
            n11,n12,n21,n22, fisher
            );
        if ( args.convert ) printf("\t%s", args.str.s);
        printf("\n");
    }
    if ( keep_als )
    {
        kbs_delete(args.rm_als, 0);
        bcf_unpack(rec,BCF_UN_ALL);
        if ( bcf_remove_allele_set(args.hdr, rec, args.rm_als)!=0 )
            error("Failed to subset alleles\n");
        return rec;
    }
    return NULL;
}

void destroy(void)
{
    if ( !args.clean_vcf )
    {
        printf("# SN, Summary Numbers\t[2]Number of Pairs\t[3]Number of Sites\t[4]Number of comparisons\t[5]P-value output threshold\n");
        printf("SN\t%d\t%"PRId64"\t%"PRId64"\t%e\n",args.npair,args.nsite,args.ncmp,args.th);
    }
    if ( args.rm_als ) kbs_destroy(args.rm_als);
    if (args.convert) convert_destroy(args.convert);
    free(args.str.s);
    free(args.pair);
    free(args.ad_arr);
}
