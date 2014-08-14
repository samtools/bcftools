/*  plugins/missing2ref.c -- sets missing genotypes to reference allele.

    Copyright (C) 2014 Genome Research Ltd.

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
#include <inttypes.h>

bcf_hdr_t *in_hdr, *out_hdr;
int32_t *gts = NULL, mgts = 0;
uint64_t nchanged = 0;

const char *about(void)
{
    return "Set missing genotypes (\"./.\") to ref allele (\"0/0\").\n";
}

int init(const char *opts, bcf_hdr_t *in, bcf_hdr_t *out)
{
    in_hdr  = in;
    out_hdr = out;
    return 0;
}

int process(bcf1_t *rec)
{
    int ngts = bcf_get_genotypes(in_hdr, rec, &gts, &mgts);

    int i, changed = 0;
    for (i=0; i<ngts; i++)
    {
        if ( gts[i]==bcf_gt_missing )
        {
            gts[i] = bcf_gt_unphased(0);
            changed++;
        }
    }
    nchanged += changed;
    if ( changed ) bcf_update_genotypes(out_hdr, rec, gts, ngts);
    return 0;
}

void destroy(void)
{
    fprintf(stderr,"Filled %"PRId64" REF alleles\n", nchanged);
    free(gts);
}


