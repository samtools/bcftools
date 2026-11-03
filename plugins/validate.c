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
#include <htslib/vcf.h>
#include <inttypes.h>
#include <getopt.h>
#include "bcftools.h"

typedef struct
{
    bcf_hdr_t *hdr;
    int marr, nerr_max, nerr;
    int32_t *arr;
}
args_t;

static args_t args;

const char *about(void)
{
    return "Catch some of the common errors in VCF files.\n";
}

const char *usage(void)
{
    return 
        "About: Catch some of the common errors in VCF files\n"
        "       - ranges of GT indices\n"
        "Usage: bcftools +validate [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -n, --nrecords <int>        Number of records to report before exiting [1]\n"
        "\n"
        "Example:\n"
        "   bcftools +validate file.vcf -- -n 10\n"
        "\n";
}


int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    memset(&args,0,sizeof(args_t));
    args.hdr  = in;
    args.nerr_max = 1;

    int c;
    static struct option loptions[] =
    {
        {"nrecords",1,0,'n'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "?hn:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'n': args.nerr_max = atoi(optarg);  break;
            case 'h':
            case '?':
            default: fprintf(stderr,"%s", usage()); exit(1); break;
        }
    }
    return 1;
}

void check_genotypes(args_t *args, bcf1_t *rec)
{
    if ( !rec->n_sample ) return;

    int ngts = bcf_get_genotypes(args->hdr, rec, &args->arr, &args->marr);
    ngts /= rec->n_sample;
    int i, j;
    
    for (i=0; i<rec->n_sample; i++)
    {
        int32_t *ptr = args->arr + i*ngts;
        for (j=0; j<ngts; j++)
        {
            if ( ptr[j]==bcf_int32_vector_end ) break;
            if ( bcf_gt_is_missing(ptr[j]) ) continue;
            int idx = bcf_gt_allele(ptr[j]);
            if ( idx < 0 || idx >= rec->n_allele )
            {
                args->nerr++;
                fprintf(stderr,"%s:%d\tsample %s has incorrect allele index %d\n", bcf_seqname(args->hdr,rec),rec->pos+1,args->hdr->samples[i],idx);
                break;
            }
        }
    }
}

bcf1_t *process(bcf1_t *rec)
{
    check_genotypes(&args, rec);

    if ( args.nerr >= args.nerr_max )
        error("Exiting, at least %d error%s found\n", args.nerr,args.nerr==1?"":"s");

    return NULL;
}

void destroy(void)
{
    if ( args.nerr )
        error("Found %d error%s\n", args.nerr,args.nerr==1?"":"s");

    free(args.arr);
}
