/*  gvcf.c -- support for gVCF files.

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

#include "call.h"

void gvcf_write(htsFile *fh, gvcf_t *gvcf, bcf_hdr_t *hdr, bcf1_t *rec, int is_ref)
{
    int i, ret, nsmpl = bcf_hdr_nsamples(hdr);

    // Flush gVCF block if chr changed, non-ref call encountered, depth is too
    // low, or no more records to come
    if ( rec && is_ref )
    {
        bcf_unpack(rec, BCF_UN_ALL);

        // per-sample depth
        ret = bcf_get_format_int32(hdr, rec, "DP", &gvcf->dp, &gvcf->mdp);
        if ( ret==nsmpl )
        {
            for (i=0; i<nsmpl; i++)
                if ( gvcf->dp[i] < gvcf->min_dp ) break;
            if ( i<nsmpl )
            {
                is_ref = 0;  // the depth is too low
                rec = NULL;
            }
        }
    }

    if ( gvcf->rid!=-1 && (!rec || gvcf->rid!=rec->rid || !is_ref || rec->pos > gvcf->end+1) )
    {
        // mpileup can output two records with the same position, SNP and
        // indel. Make sure the end position does not include the non-variant
        // SNP position just before the indel.
        if ( rec && rec->rid==gvcf->rid && rec->pos==gvcf->end ) gvcf->end--;

        gvcf->end++;    // from 0-based to 1-based coordinate

        bcf_clear1(gvcf->line);
        gvcf->line->rid  = gvcf->rid;
        gvcf->line->pos  = gvcf->start;
        gvcf->line->rlen = gvcf->end - gvcf->start;
        bcf_update_alleles_str(hdr, gvcf->line, gvcf->ref);
        bcf_update_info_int32(hdr, gvcf->line, "END", &gvcf->end, 1);
        bcf_update_genotypes(hdr, gvcf->line, gvcf->gt, nsmpl*2);
        bcf_write1(fh, hdr, gvcf->line);

        gvcf->rid = -1;
    }

    if ( !rec ) return;

    if ( is_ref )
    {
        if ( gvcf->rid==-1 )
        {
            gvcf->rid    = rec->rid;
            gvcf->start  = rec->pos;
            gvcf->ref[0] = rec->d.allele[0][0];
            gvcf->ref[1] = 0;
        }
        gvcf->end = rec->pos;
        return;
    }

    bcf_write1(fh, hdr, rec);
}

