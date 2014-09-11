/*  missing2majorPhased.c -- sets missing genotypes to reference allele.

    Copyright (C) 2014 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>
            Warren Kretzschmar <winni@well.ox.ac.uk>

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
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <inttypes.h>


bcf_hdr_t *in_hdr, *out_hdr;
int *arr = NULL, marr = 0;
int32_t *gts = NULL, mgts = 0;
uint64_t nchanged = 0;

const char *about(void)
{
    return "Set missing genotypes (\"./.\") to phased major allele (\"0|0\" or \"1|1\") depending on whether AC/AN > 0.5.\n";
}

int init(const char *opts, bcf_hdr_t *in, bcf_hdr_t *out)
{
    in_hdr  = in;
    out_hdr = out;
    bcf_hdr_append(out_hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
    bcf_hdr_append(out_hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    return 0;
}

int process(bcf1_t *rec)
{
    hts_expand(int,rec->n_allele,marr,arr);

    // this does not handle multiallelic sites yet
    if(rec->n_allele > 2) return 1;
    
    int ret = bcf_calc_ac(in_hdr,rec,arr,BCF_UN_FMT);
    int an = 0;
    if ( ret>0 )
    {
        int i = 0;
        for (i=0; i<rec->n_allele; i++) an += arr[i];
    }
    else
        return 1;

    char majorAllele = *(arr+1) *2 > an ? 1 : 0;
    
    int ngts = bcf_get_genotypes(in_hdr, rec, &gts, &mgts);

    int i, changed = 0;
    for (i=0; i<ngts; i++)
    {
        if ( gts[i]==bcf_gt_missing ) 
        {
            gts[i] = bcf_gt_phased(majorAllele);
            changed++;
        }
    }
    nchanged += changed;
    if ( changed ) bcf_update_genotypes(out_hdr, rec, gts, ngts);

    an += changed;
    if ( majorAllele ) *(arr+1) += changed;
    else *(arr) += changed;

    // update AN and AC fields
    bcf_update_info_int32(out_hdr, rec, "AN", &an, 1);
    bcf_update_info_int32(out_hdr, rec, "AC", arr+1, rec->n_allele-1);

    return 0;
}

void destroy(void) 
{
    fprintf(stderr,"Filled %"PRId64" REF alleles\n", nchanged);
    free(gts);
}


