/* The MIT License

   Copyright (c) 2015 Genome Research Ltd.

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
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include "bcftools.h"

#define SET_AN      (1<<0)
#define SET_AC      (1<<1)
#define SET_AC_Hom  (1<<2)
#define SET_AC_Het  (1<<3)
#define SET_AC_Hemi (1<<4)
#define SET_AF      (1<<5)
#define SET_NS      (1<<6)

typedef struct
{
    int nhom, nhet, nhemi, nac;
}
counts_t;

typedef struct
{
    bcf_hdr_t *in_hdr, *out_hdr;
    int tags, marr, mfarr, mcounts, gt_id;
    int32_t *arr;
    float *farr;
    counts_t *counts;
}
args_t;

static args_t args;

const char *about(void)
{
    return "Set INFO tags AF, AN, AC, NS, AC_Hom, AC_Het, AC_Hemi.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Set INFO tags AF, AN, AC, NS, AC_Hom, AC_Het, AC_Hemi.\n"
        "Usage: bcftools +fill-tags [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -t, --tags LIST         list of output tags. By default, all tags are filled.\n"
        "\n"
        "Example:\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -t AN,AC\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf\n"
        "\n";
}

int parse_tags(args_t *args, const char *str)
{
    int i, flag = 0, n_tags;
    char **tags = hts_readlist(str, 0, &n_tags);
    for(i=0; i<n_tags; i++)
    {
        if ( !strcasecmp(tags[i],"AN") ) flag |= SET_AN;
        else if ( !strcasecmp(tags[i],"AC") ) flag |= SET_AC;
        else if ( !strcasecmp(tags[i],"NS") ) flag |= SET_NS;
        else if ( !strcasecmp(tags[i],"AC_Hom") ) flag |= SET_AC_Hom;
        else if ( !strcasecmp(tags[i],"AC_Het") ) flag |= SET_AC_Het;
        else if ( !strcasecmp(tags[i],"AC_Hemi") ) flag |= SET_AC_Hemi;
        else if ( !strcasecmp(tags[i],"AF") ) flag |= SET_AF;
        else
        {
            fprintf(stderr,"Error parsing \"--tags %s\": the tag \"%s\" is not supported\n", str,tags[i]);
            exit(1);
        }
        free(tags[i]);
    }
    if (n_tags) free(tags);
    return flag;
}


int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    memset(&args,0,sizeof(args_t));
    args.in_hdr  = in;
    args.out_hdr = out;

    static struct option loptions[] =
    {
        {"tags",1,0,'t'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "?ht:T:l:cd",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 't': args.tags |= parse_tags(&args,optarg); break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( optind != argc ) error(usage());
    if ( !args.tags ) args.tags |= SET_AN|SET_AC|SET_NS|SET_AC_Hom|SET_AC_Het|SET_AC_Hemi|SET_AF;
    
    args.gt_id = bcf_hdr_id2int(args.in_hdr,BCF_DT_ID,"GT");
    if ( args.gt_id<0 ) error("Error: GT field is not present\n");

    if ( args.tags&SET_AN ) bcf_hdr_append(args.out_hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    if ( args.tags&SET_AC ) bcf_hdr_append(args.out_hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
    if ( args.tags&SET_NS ) bcf_hdr_append(args.out_hdr, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">");
    if ( args.tags&SET_AC_Hom ) bcf_hdr_append(args.out_hdr, "##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description=\"Allele counts in homozygous genotypes\">");
    if ( args.tags&SET_AC_Het ) bcf_hdr_append(args.out_hdr, "##INFO=<ID=AC_Het,Number=A,Type=Integer,Description=\"Allele counts in heterozygous genotypes\">");
    if ( args.tags&SET_AC_Hemi ) bcf_hdr_append(args.out_hdr, "##INFO=<ID=AC_Hemi,Number=A,Type=Integer,Description=\"Allele counts in hemizygous genotypes\">");
    if ( args.tags&SET_AF ) bcf_hdr_append(args.out_hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">");

    return 0;
}

bcf1_t *process(bcf1_t *rec)
{
    int i, ns = 0;

    bcf_unpack(rec, BCF_UN_FMT);
    bcf_fmt_t *fmt_gt = NULL;
    for (i=0; i<rec->n_fmt; i++)
        if ( rec->d.fmt[i].id==args.gt_id ) { fmt_gt = &rec->d.fmt[i]; break; }
    if ( !fmt_gt ) return rec;    // no GT tag

    hts_expand(int32_t,rec->n_allele,args.marr,args.arr);
    hts_expand(float,rec->n_allele,args.mfarr,args.farr);
    hts_expand(counts_t,rec->n_allele,args.mcounts,args.counts);
    memset(args.arr,0,sizeof(*args.arr)*rec->n_allele);
    memset(args.counts,0,sizeof(*args.counts)*rec->n_allele);

    #define BRANCH_INT(type_t,vector_end) { \
        for (i=0; i<rec->n_sample; i++) \
        { \
            type_t *p = (type_t*) (fmt_gt->p + i*fmt_gt->size); \
            int ial, als = 0; \
            for (ial=0; ial<fmt_gt->n; ial++) \
            { \
                if ( p[ial]==vector_end ) break; /* smaller ploidy */ \
                if ( bcf_gt_is_missing(p[ial]) ) break; /* missing allele */ \
                int idx = bcf_gt_allele(p[ial]); \
                \
                if ( idx >= rec->n_allele ) \
                    error("Incorrect allele (\"%d\") in %s at %s:%d\n",idx,args.in_hdr->samples[i],bcf_seqname(args.in_hdr,rec),rec->pos+1); \
                als |= (1<<idx);  /* this breaks with too many alleles */ \
            } \
            if ( ial==0 ) continue; /* missing alleles */ \
            ns++; \
            int is_hom  = als && !(als & (als-1)); /* only one bit is set */ \
            int is_hemi = ial==1; \
            for (ial=0; als; ial++) \
            { \
                if ( als&1 ) \
                { \
                    if ( !is_hom ) \
                        args.counts[ial].nhet++; \
                    else if ( !is_hemi ) \
                        args.counts[ial].nhom += 2; \
                    else \
                        args.counts[ial].nhemi++; \
                } \
                als >>= 1; \
            } \
        } \
    }
    switch (fmt_gt->type) {
        case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_vector_end); break;
        default: error("The GT type is not recognised: %d at %s:%d\n",fmt_gt->type, bcf_seqname(args.in_hdr,rec),rec->pos+1); break;
    }
    #undef BRANCH_INT
    if ( args.tags&SET_NS )
    {
        if ( bcf_update_info_int32(args.out_hdr,rec,"NS",&ns,1)!=0 )
            error("Error occurred while updating NS at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
    }
    if ( args.tags&SET_AN )
    {
        args.arr[0] = 0;
        for (i=0; i<rec->n_allele; i++)
            args.arr[0] += args.counts[i].nhet + args.counts[i].nhom + args.counts[i].nhemi;
        if ( bcf_update_info_int32(args.out_hdr,rec,"AN",args.arr,1)!=0 )
            error("Error occurred while updating AN at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
    }
    if ( args.tags&SET_AF )
    {
        int n = rec->n_allele-1;
        if ( n>0 )
        {
            args.arr[0] = 0;
            for (i=0; i<rec->n_allele; i++)
                args.arr[0] += args.counts[i].nhet + args.counts[i].nhom + args.counts[i].nhemi;
            for (i=1; i<rec->n_allele; i++)
                args.farr[i] = (args.counts[i].nhet + args.counts[i].nhom + args.counts[i].nhemi)*1.0/args.arr[0];
        }
        if ( args.arr[0] )
        {
            if ( bcf_update_info_float(args.out_hdr,rec,"AF",args.farr+1,n)!=0 )
                error("Error occurred while updating AF at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
        }
    }
    if ( args.tags&SET_AC )
    {
        int n = rec->n_allele-1;
        if ( n>0 )
        {
            memset(args.arr,0,sizeof(*args.arr)*rec->n_allele);
            for (i=1; i<rec->n_allele; i++)
                args.arr[i] = args.counts[i].nhet + args.counts[i].nhom + args.counts[i].nhemi;
        }
        if ( bcf_update_info_int32(args.out_hdr,rec,"AC",args.arr+1,n)!=0 )
            error("Error occurred while updating AC at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
    }
    if ( args.tags&SET_AC_Het )
    {
        int n = rec->n_allele-1;
        if ( n>0 )
        {
            memset(args.arr,0,sizeof(*args.arr)*rec->n_allele);
            for (i=1; i<rec->n_allele; i++)
                args.arr[i] += args.counts[i].nhet;
        }
        if ( bcf_update_info_int32(args.out_hdr,rec,"AC_Het",args.arr+1,n)!=0 )
            error("Error occurred while updating AC_Het at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
    }
    if ( args.tags&SET_AC_Hom )
    {
        int n = rec->n_allele-1;
        if ( n>0 )
        {
            memset(args.arr,0,sizeof(*args.arr)*rec->n_allele);
            for (i=1; i<rec->n_allele; i++)
                args.arr[i] += args.counts[i].nhom;
        }
        if ( bcf_update_info_int32(args.out_hdr,rec,"AC_Hom",args.arr+1,n)!=0 )
            error("Error occurred while updating AC_Hom at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
    }
    if ( args.tags&SET_AC_Hemi )
    {
        int n = rec->n_allele-1;
        if ( n>0 )
        {
            memset(args.arr,0,sizeof(*args.arr)*rec->n_allele);
            for (i=1; i<rec->n_allele; i++)
                args.arr[i] += args.counts[i].nhemi;
        }
        if ( bcf_update_info_int32(args.out_hdr,rec,"AC_Hemi",args.arr+1,n)!=0 )
            error("Error occurred while updating AC_Hemi at %s:%d\n", bcf_seqname(args.in_hdr,rec),rec->pos+1);
    }
    return rec;
}

void destroy(void)
{
    free(args.counts);
    free(args.arr);
    free(args.farr);
}



