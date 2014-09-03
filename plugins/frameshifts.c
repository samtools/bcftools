/*  plugins/frameshifts.c -- annotates frameshift indels.

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
#include <htslib/synced_bcf_reader.h>
#include <wordexp.h>
#include "config.h"

bcf_hdr_t *in_hdr, *out_hdr;
bcf_sr_regions_t *exons;
int32_t *frm = NULL, nfrm = 0;

const char *about(void)
{
    return
        "Annotate frameshift indels. Run as \"-p frameshifts:exons=path/to/file.tab.gz\",\n"
        "see the \"--targets-file\" man page entry for file format description.\n";
}


int init(const char *opts, bcf_hdr_t *in, bcf_hdr_t *out)
{
    in_hdr = in;
    out_hdr = out;

    int ret = bcf_hdr_append(out_hdr,"##INFO=<ID=OOF,Number=A,Type=Integer,Description=\"Frameshift Indels: out-of-frame (1), in-frame (0), not-applicable (-1 or missing)\">");
    if ( ret!=0 )
    {
        fprintf(stderr,"Error updating the header\n");
        return -1;
    }

    char *fname = config_get_string(opts,"exons");
    if ( !fname )
    {
        fprintf(stderr,"No exons given, please run as \"bcftools annotate -p frameshifts:exons=path/to/file.tab.gz\".\n");
        return -1;
    }

    wordexp_t wexp;
    wordexp(fname, &wexp, 0);
    if ( !wexp.we_wordc )
    {
        fprintf(stderr,"No such file: %s\n", fname);
        return -1;
    }
    free(fname);
    fname = wexp.we_wordv[0];

    exons = bcf_sr_regions_init(fname,1,0,1,2);
    if ( !exons )
    {
        fprintf(stderr,"Error occurred while reading (was the file compressed with bgzip?): %s\n", fname);
        return -1;
    }

    wordfree(&wexp);
    return 0;
}


int process(bcf1_t *rec)
{
    if ( rec->n_allele<2 ) return 0;    // not a variant

    int type = bcf_get_variant_types(rec);
    if ( !(type&VCF_INDEL) ) return 0;  // not an indel

    int i, len = 0;
    for (i=1; i<rec->n_allele; i++)
        if ( len > rec->d.var[i].n ) len = rec->d.var[i].n;

    int pos_to = len!=0 ? rec->pos : rec->pos - len;    // len is negative
    if ( bcf_sr_regions_overlap(exons, bcf_seqname(in_hdr,rec),rec->pos,pos_to) ) return 0;  // no overlap

    hts_expand(int32_t,rec->n_allele-1,nfrm,frm);
    for (i=1; i<rec->n_allele; i++)
    {
        if ( rec->d.var[i].type!=VCF_INDEL ) { frm[i-1] = -1; continue; }

        int len = rec->d.var[i].n, tlen = 0;
        if ( len>0 )
        {
            // insertion
            if ( exons->start <= rec->pos && exons->end > rec->pos ) tlen = abs(len);
        }
        else if ( exons->start <= rec->pos + abs(len) )
        {
            // deletion
            tlen = abs(len);
            if ( rec->pos < exons->start )              // trim the beginning
                tlen -= exons->start - rec->pos + 1;
            if ( exons->end < rec->pos + abs(len) )     // trim the end
                tlen -= rec->pos + abs(len) - exons->end;
        }
        if ( tlen )     // there are some deleted/inserted bases in the exon
        {
            if ( tlen%3 ) frm[i-1] = 1; // out-of-frame
            else frm[i-1] = 0;          // in-frame
        }
        else frm[i-1] = -1;             // not applicable (is outside)
    }

    if ( bcf_update_info_int32(out_hdr,rec,"OOF",frm,rec->n_allele-1)<0 ) return -1;
    return 0;
}


void destroy(void)
{
    bcf_sr_regions_destroy(exons);
}


