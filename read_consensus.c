/*  read_consensus.c -- create and maintain consensus of reads

    Copyright (C) 2022 Genome Research Ltd.

    Author: pd3@sanger

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

#include <read_consensus.h>
#include <assert.h>

#define NI 10 // number of alternative insertion sequences at one position in a single sample
typedef struct
{
    char *str[NI];
    int len[NI];
    int freq[NI];
}
ins_freq_t;

typedef struct
{
    int len[NI];
    int freq[NI];
}
del_freq_t;

typedef struct
{
    int base[5];    // frequencies of A,C,G,T,N
}
base_freq_t;

struct _read_cns_t
{
    int pos, beg, end;
    int band;   // maximum absolute deviation from the diagonal, used for BAQ alignment
    base_freq_t *base_freq;
    ins_freq_t *ins_freq;
    del_freq_t *del_freq;
    char *stmp;
    int mstmp, mfreq;   // allocated size of stmp and *_freq arrays
};

void rcns_destroy(read_cns_t *rcns)
{
    int i,j;
    for (i=0; i<rcns->mfreq; i++)
    {
        ins_freq_t *ifrq = &rcns->ins_freq[i];
        for (j=0; j<NI && ifrq->str[j]; j++) free(ifrq->str[j]);
    }
    free(rcns->ins_freq);
    free(rcns->del_freq);
    free(rcns->base_freq);
    free(rcns->stmp);
    free(rcns);
}
int _rcns_init_arrays(read_cns_t *rcns)
{
    int i,j,n = rcns->end - rcns->beg + 1;
    if ( n > rcns->mfreq )
    {
        ins_freq_t *ifrq = (ins_freq_t*) realloc(rcns->ins_freq,sizeof(*rcns->ins_freq)*n);
        if ( !ifrq ) return -1;
        rcns->ins_freq = ifrq;
        memset(ifrq+rcns->mfreq,0,sizeof(*rcns->ins_freq)*(n-rcns->mfreq));

        del_freq_t *dfrq = (del_freq_t*) realloc(rcns->del_freq,sizeof(*rcns->del_freq)*n);
        if ( !dfrq ) return -1;
        rcns->del_freq = dfrq;
        memset(dfrq+rcns->mfreq,0,sizeof(*rcns->del_freq)*(n-rcns->mfreq));

        base_freq_t *bfrq = (base_freq_t*) realloc(rcns->base_freq,sizeof(*rcns->base_freq)*n);
        if ( !bfrq ) return -1;
        rcns->base_freq = bfrq;
        memset(bfrq+rcns->mfreq,0,sizeof(*rcns->base_freq)*(n-rcns->mfreq));

        rcns->mfreq = n;
    }
    memset(rcns->base_freq,0,sizeof(*rcns->base_freq)*n);
    memset(rcns->del_freq,0,sizeof(*rcns->del_freq)*n);
    for (i=0; i<n; i++)
    {
        ins_freq_t *ifrq = &rcns->ins_freq[i];
        for (j=0; j<NI && ifrq->str[j]; j++) free(ifrq->str[j]);
    }
    memset(rcns->ins_freq,0,sizeof(*rcns->ins_freq)*n);
    return 0;
}
int rcns_reset(read_cns_t *rcns, int pos, int beg, int end)
{
    rcns->band = 0;
    rcns->pos  = pos;
    rcns->beg  = beg;
    rcns->end  = end;
    return _rcns_init_arrays(rcns);
}

static inline void _rcns_add_base(read_cns_t *rcns, int ref_pos, int nt16)
{
    int i = ref_pos - rcns->beg;
    rcns->base_freq[i].base[seq_nt16_int[nt16]]++;
}
static void _rcns_add_ins(read_cns_t *rcns, int ref_pos, uint8_t *nt16_seq, int len)
{
    int i = ref_pos - rcns->beg;
    ins_freq_t *ifrq = &rcns->ins_freq[i];

    char *str;
    if ( rcns->mstmp < len )
    {
        str = realloc(rcns->stmp,len*sizeof(*str));
        if ( !str ) return;
        rcns->mstmp = len;
        rcns->stmp  = str;
    }
    else
        str = rcns->stmp;
    for (i=0; i<len; i++) str[i] = seq_nt16_int[bam_seqi(nt16_seq,i)];

    for (i=0; i<NI && ifrq->str[i]; i++)
        if ( ifrq->len[i]==len && !memcmp(ifrq->str[i],str,len) ) break;

    assert( i<NI ); // how frequent is it?
    if ( i>=NI ) return;    // too many choices; discard

    if ( !ifrq->str[i] )      // new insertion
    {
        if ( !(ifrq->str[i]=malloc(len)) ) return;
        memcpy(ifrq->str[i], str, len);
        ifrq->len[i] = len;
    }
    ifrq->freq[i]++;
}
static void _rcns_add_del(read_cns_t *rcns, int ref_pos, int len)
{
    int i = ref_pos - rcns->beg;
    del_freq_t *dfrq = &rcns->del_freq[i];

    for (i=0; i<NI && dfrq->len[i]; i++)
        if ( dfrq->len[i]==len ) break;

    assert( i<NI ); // how frequent is it?
    if ( i>=NI ) return;    // too many choices; discard

    if ( !dfrq->len[i] ) dfrq->len[i] = len;    // new deletion
    dfrq->freq[i]++;
}

read_cns_t *rcns_init(int pos, int beg, int end)
{
    read_cns_t *rcns = (read_cns_t*) calloc(1,sizeof(read_cns_t));
    rcns->pos  = pos;
    rcns->beg  = beg;
    rcns->end  = end;
    if ( _rcns_init_arrays(rcns)!=0 )
    {
        rcns_destroy(rcns);
        return NULL;
    }
    return rcns;
}

int rcns_set_reads(read_cns_t *rcns, bam_pileup1_t *plp, int nplp)
{
    // fill consensus arrays
    int i,j,k, local_band_max = 0;  // maximum absolute deviation from diagonal
    for (i=0; i<nplp; i++) // for each read...
    {
        const bam_pileup1_t *p = plp + i;
        bam1_t *b = p->b;
        int x = b->core.pos;  // ref coordinate
        int y = 0;            // seq coordinate
        uint32_t *cigar = bam_get_cigar(b);
        uint8_t *seq = bam_get_seq(b);

        int local_band = 0; // current deviation from diagonal
        for (k = 0; k < b->core.n_cigar; ++k)
        {
            int op  = cigar[k] &  BAM_CIGAR_MASK;
            int len = cigar[k] >> BAM_CIGAR_SHIFT;
            if ( op==BAM_CSOFT_CLIP ) y += len;
            else if ( op==BAM_CMATCH || op==BAM_CEQUAL || op==BAM_CDIFF )
            {
                int j_end = rcns->end < x + len - 1 ? rcns->end - x : len - 1;
                int j_beg = rcns->beg > x ? rcns->beg - x : 0;
                x += j_beg;     // ref pos
                y += j_beg;     // seq pos
                for (j=j_beg; j<=j_end; j++, x++, y++) _rcns_add_base(rcns,x,bam_seqi(seq,y));
            }
            else if ( op==BAM_CINS )
            {
                if ( x>=rcns->beg && x<rcns->end )
                {
                    local_band += p->indel;
                    _rcns_add_ins(rcns,x,seq,len);
                }
            }
            else if ( op==BAM_CDEL )
            {
                if ( x>=rcns->beg && x+len-1<=rcns->end )
                {
                    local_band += -p->indel;
                    _rcns_add_del(rcns,x,len);
                }
            }
            if ( local_band_max < local_band ) local_band_max = local_band;
        }

        // Track the biggest deviation +/- from diagonal, used in BAQ alignment step.
        if ( rcns->band < local_band_max ) rcns->band = local_band_max;
    }

    return 0;
}

// The algorithm:
//  1. Identify variant positions
//  2. Sort variants by abs(variant_allele_freq-0.5) in descending order
//  3. Take the top sorted variants (up to 8 to fit in uint8_t) and count the number of
//      corresponding reads to create frequency spectrum
//  4. Correct errors, collapse to the requested number of haplotypes (consensus sequences)
//      using majority vote for the distribution tail
const char **rcns_get_consensus(read_cns_t *rcns, const char *ref, int **npos)
{
    int i,j,n = rcns->end - rcns->beg + 1;
    for (i=0; i<n; i++)
    {
        fprintf(stderr,"%d\t",rcns->beg+i+1);
        for (j=0; j<5; j++)
            fprintf(stderr,"\t%d",rcns->base_freq[i].base[j]);
        fprintf(stderr,"\n");
    }
    return NULL;
}

