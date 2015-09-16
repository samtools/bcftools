/*  plugins/setGT.c -- set gentoypes to given values

    Copyright (C) 2015 Genome Research Ltd.

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <inttypes.h>
#include <getopt.h>
#include "bcftools.h"

bcf_hdr_t *in_hdr, *out_hdr;
int32_t *gts = NULL, mgts = 0;
int *arr = NULL, marr = 0;
uint64_t nchanged = 0;
int tgt_mask = 0, new_mask = 0, new_gt = 0;

#define GT_MISSING   1
#define GT_PARTIAL  (1<<1)
#define GT_REF      (1<<2)
#define GT_MAJOR    (1<<3)
#define GT_PHASED   (1<<4)
#define GT_UNPHASED (1<<5)
#define GT_ALL      (1<<6)

const char *about(void)
{
    return "Set genotypes: partially missing to missing, missing to ref/major allele, etc.\n";
}

const char *usage(void)
{
    return 
        "About: Sets genotypes. The target genotypes can be specified as:\n"
        "           ./.  .. completely missing (\".\" or \"./.\", depending on ploidy)\n"
        "           ./x  .. partially missing (e.g., \"./0\" or \".|1\" but not \"./.\")\n"
        "           .    .. partially or completely missing\n"
        "           a    .. all genotypes\n"
        "       and the new genotype can be one of:\n"
        "           .    .. missing (\".\" or \"./.\", keeps ploidy)\n"
        "           0    .. reference allele\n"
        "           M    .. major allele\n"
        "           p    .. phased genotype\n"
        "           u    .. unphase genotype and sort by allele (1|0 becomes 0/1)\n"
        "Usage: bcftools +setGT [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -n, --new-gt <type>         Genotypes to set, see above\n"
        "   -t, --target-gt <type>      Genotypes to change, see above\n"
        "\n"
        "Example:\n"
        "   # set missing genotypes (\"./.\") to phased ref genotypes (\"0|0\")\n"
        "   bcftools +setGT in.vcf -- -t . -n 0p\n"
        "\n"
        "   # set partially missing genotypes to completely missing\n"
        "   bcftools +setGT in.vcf -- -t ./x -n .\n"
        "\n";
}


int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    int c;
    static struct option loptions[] =
    {
        {"new-gt",1,0,'n'},
        {"target-gt",0,0,'t'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "?hn:t:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'n': new_mask = bcf_gt_phased(0); 
                if ( strchr(optarg,'.') ) new_mask |= GT_MISSING;
                if ( strchr(optarg,'0') ) new_mask |= GT_REF;
                if ( strchr(optarg,'M') ) new_mask |= GT_MAJOR;
                if ( strchr(optarg,'p') ) new_mask |= GT_PHASED;
                if ( strchr(optarg,'u') ) new_mask |= GT_UNPHASED;
                break;
            case 't':
                if ( !strcmp(optarg,".") ) tgt_mask |= GT_MISSING|GT_PARTIAL;
                if ( !strcmp(optarg,"./x") ) tgt_mask |= GT_PARTIAL;
                if ( !strcmp(optarg,"./.") ) tgt_mask |= GT_MISSING;
                if ( !strcmp(optarg,"a") ) tgt_mask |= GT_ALL;
                break;
            case 'h':
            case '?':
            default: fprintf(stderr,"%s", usage()); exit(1); break;
        }
    }
    in_hdr  = in;
    out_hdr = out;

    if ( !new_mask ) error("Expected -n option\n");
    if ( !tgt_mask ) error("Expected -t option\n");

    if ( new_mask & GT_MISSING ) new_gt = bcf_gt_missing;
    if ( new_mask & GT_REF ) new_gt = new_mask&GT_PHASED ? bcf_gt_phased(0) : bcf_gt_unphased(0);

    return 0;
}

bcf1_t *process(bcf1_t *rec)
{
    if ( !rec->n_sample ) return rec;

    int ngts = bcf_get_genotypes(in_hdr, rec, &gts, &mgts);
    ngts /= rec->n_sample;
    int i, j, changed = 0;
    
    // Calculating allele frequency for each allele and determining major allele
    // only do this if use_major is true
    int an = 0, maxAC = -1, majorAllele = -1;
    if ( new_mask & GT_MAJOR )
    {
        hts_expand(int,rec->n_allele,marr,arr);
        int ret = bcf_calc_ac(in_hdr,rec,arr,BCF_UN_FMT);
        if ( ret<= 0 )
            error("Could not calculate allele count at %s:%d\n", bcf_seqname(in_hdr,rec),rec->pos+1);

        for(i=0; i < rec->n_allele; ++i)
        {
            an += arr[i];
            if (arr[i] > maxAC)
            {
                maxAC = arr[i];
                majorAllele = i;
            }
        }

        // replacing new_gt by major allele
        new_gt = new_mask & GT_PHASED ?  bcf_gt_phased(majorAllele) : bcf_gt_unphased(majorAllele);
    }

    // replace gts
    for (i=0; i<rec->n_sample; i++)
    {
        int ploidy = 0, nmiss = 0;
        int32_t *ptr = gts + i*ngts;
        for (j=0; j<ngts; j++)
        {
            if ( ptr[j]==bcf_int32_vector_end ) break;
            ploidy++;
            if ( ptr[j]==bcf_gt_missing ) nmiss++;
        }

        int do_set = 0;
        if ( tgt_mask&GT_ALL ) do_set = 1;
        else if ( tgt_mask&GT_PARTIAL && nmiss ) do_set = 1;
        else if ( tgt_mask&GT_MISSING && ploidy==nmiss ) do_set = 1;

        if ( !do_set ) continue;

        if ( new_mask&GT_UNPHASED )
        {
            for (j=0; j<ngts; j++)
            {
                if ( ptr[j]==bcf_int32_vector_end ) break;
                if ( !bcf_gt_is_phased(ptr[j]) ) continue;
                ptr[j] = bcf_gt_unphased(bcf_gt_allele(ptr[j]));    // remove phasing
                changed++;
            }

            // insertion sort
            int k, l;
            for (k=1; k<j; k++)
            {
                int32_t x = ptr[k];
                l = k;
                while ( l>0 && ptr[l-1]>x )
                {
                    ptr[l] = ptr[l-1];
                    l--;
                }
                ptr[l] = x;
            }
            continue;
        }

        for (j=0; j<ngts; j++)
        {
            if ( ptr[j]==bcf_int32_vector_end ) break;
            ptr[j] = new_gt;
            changed++;
        }
    }
    nchanged += changed;
    if ( changed ) bcf_update_genotypes(out_hdr, rec, gts, ngts*rec->n_sample);
    return rec;
}

void destroy(void)
{
    free(arr);
    fprintf(stderr,"Filled %"PRId64" alleles\n", nchanged);
    free(gts);
}


