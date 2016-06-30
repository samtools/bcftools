/*  bam_sample.c -- group data by sample.

    Copyright (C) 2010, 2011 Broad Institute.
    Copyright (C) 2013, 2016 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#include <stdlib.h>
#include <string.h>
#include <htslib/hts.h>
#include <htslib/khash_str2int.h>
#include "bam_sample.h"

bam_sample_t *bam_smpl_init(void)
{
    bam_sample_t *s;
    s = (bam_sample_t*) calloc(1, sizeof(bam_sample_t));
    s->rg2smid = khash_str2int_init();
    s->sm2id = khash_str2int_init();
    return s;
}

void bam_smpl_destroy(bam_sample_t *sm)
{
    if ( !sm ) return;
    if ( sm->rg2smid ) khash_str2int_destroy_free(sm->rg2smid);
    if ( sm->sm2id ) khash_str2int_destroy_free(sm->sm2id);
    free(sm->smpl);
    free(sm);
}

static void add_pair(bam_sample_t *sm, void *sm2id, const char *readgroup, const char *sample)
{
    if ( khash_str2int_has_key(sm->rg2smid,readgroup) ) return;   // duplicated @RG-ID
    
    int ismpl;
    if ( khash_str2int_get(sm2id, sample, &ismpl) < 0 )
    {
        // the sample encountered for the first time
        ismpl = sm->n++;
        hts_expand0(char*,sm->n,sm->m,sm->smpl);
        sm->smpl[ismpl] = strdup(sample);
        khash_str2int_set(sm2id, sm->smpl[ismpl], ismpl);
    }
    khash_str2int_set(sm->rg2smid, strdup(readgroup), ismpl);
}

int bam_smpl_add(bam_sample_t *sm, const char *fn, const char *txt, void *sample_list, int sample_logic, void *white_hash)
{
    const char *p = txt, *q, *r;
    kstring_t buf, first_sm;
    int n = 0;
    if (txt == 0) {
        add_pair(sm, sm->sm2id, fn, fn);
        return 0;
    }
    memset(&buf, 0, sizeof(kstring_t));
    memset(&first_sm, 0, sizeof(kstring_t));
    while ((q = strstr(p, "@RG")) != 0) {
        p = q + 3;
        r = q = 0;
        if ((q = strstr(p, "\tID:")) != 0) q += 4;
        if ((r = strstr(p, "\tSM:")) != 0) r += 4;
        if (r && q) {
            char *u, *v;
            int ioq, ior;
            for (u = (char*)q; *u && *u != '\t' && *u != '\n'; ++u);
            for (v = (char*)r; *v && *v != '\t' && *v != '\n'; ++v);
            ioq = *u; ior = *v; *u = *v = '\0';

            // r now points to a null terminated sample name
            int accept_rg = 1;
            if ( sample_list )
            {
                accept_rg = khash_str2int_has_key(sample_list,r);
                if ( sample_logic==0 ) accept_rg = accept_rg ? 0 : 1;
            }
            if ( accept_rg )
            {
                buf.l = 0; kputs(fn, &buf); kputc('/', &buf); kputs(q, &buf);
                add_pair(sm, sm->sm2id, buf.s, r);
                if ( !first_sm.s )
                    kputs(r,&first_sm);
                if ( sample_list ) khash_str2int_inc(white_hash,strdup(q));
            }
            *u = ioq; *v = ior;
        } else break;
        p = q > r? q : r;
        ++n;
    }
    if (n == 0) add_pair(sm, sm->sm2id, fn, fn);
    // If there is only one RG tag present in the header and reads are not annotated, don't refuse to work but
    //  use the tag instead.
    else if ( n==1 && first_sm.s )
        add_pair(sm,sm->sm2id,fn,first_sm.s);
    if ( first_sm.s )
        free(first_sm.s);

    free(buf.s);
    return 0;
}

int bam_smpl_rg2smid(const bam_sample_t *sm, const char *fn, const char *rg, kstring_t *str)
{
    int ismpl;
    if ( rg )
    {
        str->l = 0;
        kputs(fn, str); kputc('/', str); kputs(rg, str);
        if ( khash_str2int_get(sm->rg2smid, str->s, &ismpl) < 0 ) return -1;
        return ismpl;
    }
    if ( khash_str2int_get(sm->rg2smid, fn, &ismpl) < 0 ) return -1;
    return ismpl;
}
