/*  bam2bcf_iaux.c -- modified indel caller

    Copyright (C) 2022 Genome Research Ltd.

    Author: pd3@sanger, jkb

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
    DEALINGS IN THE SOFTWARE
*/

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/khash_str2int.h>
#include "bam2bcf.h"

#include <htslib/ksort.h>
KSORT_INIT_STATIC_GENERIC(uint32_t)


#define MAX_TYPES 64

typedef struct
{
    int pos;    // current position
    char *chr;  // current chromosome
    int nsmpl;  // number of samples
    int *nplp;              // per-sample number of reads
    bam_pileup1_t **plp;    // per-sample reads
    bcf_callaux_t *bca;     // auxiliary bam2bcf structure
    const char *ref;        // reference genome (ASCII)
    uint32_t *uitmp;        // temporary unsigned int array
    uint32_t *inscns;       // insertions consensus [itype*max_ins_len+i]
    int muitmp, minscns;    // size of uitmp, inscns
    int iref_type, ntypes, types[MAX_TYPES];   // indel types
    int max_ins_len;
}
indel_aux_t;

static void _print_types(indel_aux_t *iaux)
{
    int i,j;
    fprintf(stderr,"types at %s:%d ... ",iaux->chr,iaux->pos+1);
    for (i=0; i<iaux->ntypes; i++)
    {
        if ( iaux->types[i]<=0 )
        {
            if ( i==iaux->iref_type ) fprintf(stderr," %d",iaux->types[i]);
            continue;
        }
        fprintf(stderr," ");
        uint32_t *cns = &iaux->inscns[i*iaux->max_ins_len];
        for (j=0; j<iaux->types[i]; j++) fprintf(stderr,"%c","ACGTN"[cns[j]]);
    }
    fprintf(stderr,"\n");
}

static int _have_indel_reads(indel_aux_t *iaux)
{
    int i,j;
    for (i=0; i<iaux->nsmpl; i++)
    {
        for (j=0; j<iaux->nplp[i]; j++)
            if ( iaux->plp[i][j].indel ) return 1;
    }
    return 0;
}

static int iaux_init_ins_types(indel_aux_t *iaux)
{
    if ( !iaux->max_ins_len ) return 0;

    uint32_t *aux;
    int naux = 5 * iaux->ntypes * iaux->max_ins_len;
    if ( iaux->muitmp < naux )
    {
        aux = (uint32_t*) realloc(iaux->uitmp,naux*sizeof(*aux));
        if ( !aux ) return -1;
        iaux->uitmp  = aux;
        iaux->muitmp = naux;
    }
    else aux = iaux->uitmp;
    memset(aux,0,naux*sizeof(*aux));

    // count the number of occurrences of each base at each position for each type of insertion
    int t,s,i,j;
    for (t=0; t<iaux->ntypes; t++)
    {
        if ( iaux->types[t] <= 0) continue;
        for (s=0; s<iaux->nsmpl; s++)
        {
            for (i=0; i<iaux->nplp[s]; i++)
            {
                bam_pileup1_t *plp = iaux->plp[s] + i;
                if ( plp->indel != iaux->types[t] ) continue;
                uint8_t *seq = bam_get_seq(plp->b);
                for (j=0; j<plp->indel; j++)
                {
                    int c = seq_nt16_int[bam_seqi(seq, plp->qpos+j+1)];
                    assert(c<5);
                    aux[5*(t*iaux->max_ins_len+j) + c]++;
                }
            }
        }
    }

    uint32_t *cns;
    int ncns = iaux->ntypes * iaux->max_ins_len;
    if ( iaux->minscns < ncns )
    {
        cns = (uint32_t*) realloc(iaux->inscns,naux*sizeof(*aux));
        if ( !cns ) return -1;
        iaux->inscns  = cns;
        iaux->minscns = ncns;
    }
    else cns = iaux->inscns;
    memset(aux,0,ncns*sizeof(*cns));

    // use the majority rule to construct the consensus
    for (t=0; t<iaux->ntypes; t++)
    {
        for (i=0; i<iaux->types[t]; i++)
        {
            uint32_t *tmp = &aux[5*(t*iaux->max_ins_len+i)], max = tmp[0], max_j = 0;
            for (j=1; j<5; j++)
                if ( max < tmp[j] ) max = tmp[j], max_j = j;
            cns[t*iaux->max_ins_len + i] = max ? max_j : 4;
            if ( max_j==4 ) { iaux->types[t] = 0; break; } // discard insertions which contain N's
        }
    }
    return 0;
}

#define MINUS_CONST 0x10000000
static int iaux_init_types(indel_aux_t *iaux)
{
    if ( !_have_indel_reads(iaux) ) return -1;

    iaux->bca->max_support = 0;

    int i,j, nreads = 0;
    for (i=0; i<iaux->nsmpl; i++) nreads += iaux->nplp[i];

    uint32_t *aux;
    if ( iaux->muitmp < nreads+1 )
    {
        aux = (uint32_t*) realloc(iaux->uitmp,(nreads+1)*sizeof(*iaux->uitmp));
        if ( !aux ) return -1;
        iaux->uitmp  = aux;
        iaux->muitmp = nreads+1;
    }
    else aux = iaux->uitmp;
    memset(aux,0,(nreads+1)*sizeof(*aux));

    int naux = 0, indel_support_ok = 0, n_alt = 0, n_tot = 0;
    int max_rd_len = 0;   // max sequence length that includes ref+del bases

    // Fill out aux[] array with all the non-zero indel sizes.
    aux[naux++] = MINUS_CONST;  // zero indel is always a type (REF)
    for (i=0; i<iaux->nsmpl; i++)
    {
        int nalt = naux, ntot = 0;  // per sample values
        for (j=0; j<iaux->nplp[i]; j++)
        {
            const bam_pileup1_t *plp = iaux->plp[i] + j;
            ntot++;
            if ( plp->indel ) aux[naux++] = MINUS_CONST + plp->indel;
            if ( !PLP_QLEN(&plp->cd) ) PLP_QLEN(&plp->cd) = bam_cigar2qlen(plp->b->core.n_cigar, bam_get_cigar(plp->b));
            if ( PLP_QLEN(&plp->cd) > max_rd_len ) max_rd_len = PLP_QLEN(&plp->cd);
        }
        nalt = naux - nalt;
        if ( iaux->bca->per_sample_flt )
        {
            double frac = (double)nalt/naux;
            if ( nalt >= iaux->bca->min_support && frac >= iaux->bca->min_frac ) indel_support_ok = 1;
            if ( nalt > iaux->bca->max_support && frac > 0 ) iaux->bca->max_support = nalt, iaux->bca->max_frac = frac;
        }
        else
        {
            n_alt += nalt;
            n_tot += ntot;
        }
    }

    // Check if the minimum required number of indel reads has been observed
    if ( !iaux->bca->per_sample_flt && n_alt >= iaux->bca->min_support && (double)n_alt/n_tot >= iaux->bca->min_frac ) indel_support_ok = 1;
    if ( naux==1 || !indel_support_ok ) return -1;

    // To prevent long stretches of N's to be mistaken for indels (sometimes thousands of bases), check the number of N's in the
    // sequence and skip places where half or more reference bases in the sequence that follows pos are Ns
    int nN = 0, i_end = iaux->pos + (iaux->bca->indel_win_size < max_rd_len ? iaux->bca->indel_win_size : max_rd_len);
    for (i=iaux->pos; i<i_end && iaux->ref[i]; i++)
        if ( iaux->ref[i] == 'N' ) nN++;
    if ( 2*nN > i - iaux->pos ) return -1;

    // Sort aux[] and dedup
    int n_types = 1;
    ks_introsort(uint32_t, naux, aux);
    for (i=1; i<naux; i++)
        if ( aux[i] != aux[i-1] ) n_types++;

    if ( n_types >= MAX_TYPES )
    {
        static int warned = 0;
        if ( !warned )
        {
            fprintf(stderr, "Warning: excessive number of INDEL alleles at %s:%d, skipping. (This warning is printed only once)\n",iaux->chr,iaux->pos+1);
            warned = 1;
        }
        return -1;
    }

    // Finally fill out the types[] array detailing the size of insertion or deletion.
    iaux->ntypes = 0;
    iaux->max_ins_len = 0;
    for (i=0; i<naux; i++)
    {
        int isize = (int32_t)(aux[i] - MINUS_CONST);
        for (j=i+1; j<naux; j++)
            if ( aux[j] != aux[i] ) break;

        // Only include the REF type and types with sufficient support. Note that the position
        // already passed, this is just to reduce the number of indel types. The check is
        // permissive, the thresholds min_support and min_frac are not enforced in per-sample mode
        int is_ok = 0;
        if ( !isize )
        {
            is_ok = 1;
            iaux->iref_type = iaux->ntypes;
        }
        else
        {
            if ( j-i >= iaux->bca->min_support ) is_ok = 1;
            if ( !iaux->bca->per_sample_flt && (double)(j-i) / n_tot < iaux->bca->min_frac ) is_ok = 0;
        }
        if ( is_ok )
        {
            iaux->types[iaux->ntypes++] = isize;
            if ( isize > 0 && isize > iaux->max_ins_len ) iaux->max_ins_len = isize;
        }
        i = j-1;
    }
    if ( iaux->ntypes <= 1 ) return -1;

    // Init insertion types
    if ( iaux_init_ins_types(iaux) < 0 ) return -1;

    _print_types(iaux);

    return iaux->ntypes;
}
#undef MINUS_CONST

static int iaux_set_consensus(indel_aux_t *iaux, int ismpl)
{
    return -1;
}
static int iaux_score_reads(indel_aux_t *iaux, int ismpl, int itype)
{
    return -1;
}
static int iaux_compute_quality(indel_aux_t *iaux)
{
    return -1;
}

/*
FIXME: with high number of samples, do we handle IMF correctly?  Is it
fraction of indels across entire data set, or just fraction for this
specific sample? Needs to check bca->per_sample_flt (--per-sample-mF) opt.
 */

/*
    notes:
    - n .. number of samples
    - the routine sets bam_pileup1_t.aux of each read as follows:
        - 6: unused
        - 6: the call; index to bcf_callaux_t.indel_types   .. (aux>>16)&0x3f
        - 8: estimated sequence quality                     .. (aux>>8)&0xff
        - 8: indel quality                                  .. aux&0xff
 */
int bcf_iaux_gap_prep(int n, int *n_plp, bam_pileup1_t **plp, int pos, bcf_callaux_t *bca, const char *ref)
{
assert(!(ref == 0 || bca == 0));    // can this ever happen? when?
    if (ref == 0 || bca == 0) return -1;

    if ( !bca->iaux ) bca->iaux = calloc(1,sizeof(indel_aux_t));
    indel_aux_t *iaux = bca->iaux;
    iaux->nsmpl = n;
    iaux->nplp  = n_plp;
    iaux->plp   = plp;
    iaux->bca   = bca;
    iaux->ref   = ref;
    iaux->pos   = pos;
    iaux->chr   = bca->chr;

    // Check if there is an indel at this position and if yes, find all indel types and determine
    // window boundaries. todo: We want this information cached so that for long reads we don't keep
    // redoing the whole analysis again and again
    int ntypes = iaux_init_types(iaux);
    if ( !ntypes ) return -1;


    // Create two template consensus sequences for each sample (assuming max diploid organism).
    // Then apply each indel type on top of the templates, realign every read and remember score
    int i,j;
    for (i=0; i<iaux->nsmpl; i++)
    {
        iaux_set_consensus(iaux, i);
        for (j=0; j<ntypes; j++) iaux_score_reads(iaux, i, j);
    }
    int nalt = iaux_compute_quality(iaux);
    return nalt > 0 ? 0 : -1;
}
void bcf_iaux_destroy(bcf_callaux_t *bca)
{
    if ( !bca->iaux ) return;
    indel_aux_t *iaux = (indel_aux_t*)bca->iaux;
    free(iaux->uitmp);
    free(iaux->inscns);
    free(iaux);
}


