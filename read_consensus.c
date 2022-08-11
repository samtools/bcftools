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

#include <assert.h>
#include <math.h>
#include "read_consensus.h"
#include "cigar_state.h"
#include "kheap.h"


// Frequency arrays for each variant type
#define NI 10 // number of alternative insertion sequences at one position in a single sample
typedef struct
{
    char *nt16_seq[NI];
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
    int base[6];    // frequencies of A,C,G,T,N,deletion
}
base_freq_t;


// Candidate variants for each interesting position to build consensus haplotypes
enum variant_type { snv, ins, del, done };
typedef struct
{
    enum variant_type vtype;
    hts_pos_t pos;      // variant position (reference sequence coordinates), indels follow VCF convention
    int idx;            // temporary 0-based index to rcns.cvar
    int which,          // base/ins/del in rcns.[base|ins|del]_freq array
        depth;          // coverage at the position
    float af, af_dev;   // variant allele frequency (just for debugging printout) and absolute af deviation from 0.5
}
candidate_var_t;
static inline int cvar_not_preferred(candidate_var_t *a, candidate_var_t *b)
{
    if ( a->af_dev == b->af_dev ) return a->depth < b->depth ? 1 : 0;
    return a->af_dev > b->af_dev ? 1 : 0;
}
KHEAP_INIT(cvh, candidate_var_t, cvar_not_preferred);
typedef khp_cvh_t cvar_heap_t;

#define MAX_NCVAR 8         // This results in alloc() of 2^MAX_NCVAR possible haplotypes
#define NHAP (1<<MAX_NCVAR) // The number of possible haplotypes
struct _read_cns_t
{
    hts_pos_t pos, beg, end;    // current position and window boundaries (0-based, inclusive, ref seq coordinates)
    int band,                   // maximum absolute deviation from the diagonal, used for BAQ alignment
        max_del;                // maximum deletion lentgth starting at the tested position
    base_freq_t *base_freq;     // frequency of each variant type: base, ins, del
    ins_freq_t *ins_freq;
    del_freq_t *del_freq;
    char *stmp;             // temporary array
    int mstmp, mfreq;       // allocated size of stmp and *_freq arrays
    cvar_heap_t *cv_heap;   // heap to maintain the top MAX_NCVAR variants
    int ncvar;              // cvar and cv_heap size
    candidate_var_t cvar[MAX_NCVAR];    // candidate variants, sorted by position and type
    int hap_freq[NHAP];     // haplotype frequencies
    bam_pileup1_t *plp;     // reads to construct consensus from
    int nplp;               // number of reads in the pileup
    int cns_hap[2], ncns;   // the top two consensus haplotypes and the number of haplotypes to use
    int mcns;               // the allocated size of cns.seq and cns.pos buffers
    cns_seq_t cns[3];       // the consensus sequences to fill
};

void rcns_destroy(read_cns_t *rcns)
{
    int i,j;
    for (i=0; i<rcns->mfreq; i++)
    {
        ins_freq_t *ifrq = &rcns->ins_freq[i];
        for (j=0; j<NI && ifrq->nt16_seq[j]; j++) free(ifrq->nt16_seq[j]);
    }
    for (i=0; i<2; i++)
    {
        free(rcns->cns[i].seq);
        free(rcns->cns[i].pos);
    }
    free(rcns->ins_freq);
    free(rcns->del_freq);
    free(rcns->base_freq);
    free(rcns->stmp);
    khp_destroy(cvh,rcns->cv_heap);
    free(rcns);
}
int init_arrays(read_cns_t *rcns)
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
        for (j=0; j<NI && ifrq->nt16_seq[j]; j++) free(ifrq->nt16_seq[j]);
    }
    memset(rcns->ins_freq,0,sizeof(*rcns->ins_freq)*n);
    return 0;
}
int rcns_reset(read_cns_t *rcns, hts_pos_t pos, hts_pos_t beg, hts_pos_t end)
{
    rcns->band = 0;
    rcns->pos  = pos;
    rcns->beg  = beg;
    rcns->end  = end;
    int i;
    for (i=0; i<2; i++) rcns->cns[i].nseq = rcns->cns[i].ipos = 0;
    // this should not be necessary if the caller did run all steps
    while (rcns->cv_heap->ndat) khp_delete(cvh, rcns->cv_heap);
    return init_arrays(rcns);
}

static inline void add_base(read_cns_t *rcns, int ref_pos, int nt16)
{
    int i = ref_pos - rcns->beg;
    rcns->base_freq[i].base[seq_nt16_int[nt16]]++;
}
static void add_ins(read_cns_t *rcns, int ref_pos, int seq_pos, uint8_t *raw_seq, int len)
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
    for (i=0; i<len; i++) str[i] = bam_seqi(raw_seq,i+seq_pos);

    for (i=0; i<NI && ifrq->nt16_seq[i]; i++)
        if ( ifrq->len[i]==len && !memcmp(ifrq->nt16_seq[i],str,len) ) break;

    assert( i<NI );         // how frequent is it to have too many insertion types?
    if ( i>=NI ) return;    // too many choices; discard

    if ( !ifrq->nt16_seq[i] )      // new insertion
    {
        if ( !(ifrq->nt16_seq[i]=malloc(len)) ) return;
        memcpy(ifrq->nt16_seq[i], str, len);
        ifrq->len[i] = len;
    }
    ifrq->freq[i]++;
}
static void add_del(read_cns_t *rcns, int ref_pos, int len)
{
    int i = ref_pos - rcns->beg;
    int j,n = rcns->end - rcns->beg + 1;
    if ( i + len + 1 < n ) n = i + len + 1;
    for (j=i+1; j<n; j++)
        rcns->base_freq[j].base[5]++;

    del_freq_t *dfrq = &rcns->del_freq[i];
    for (i=0; i<NI && dfrq->len[i]; i++)
        if ( dfrq->len[i]==len ) break;

    assert( i<NI );         // how frequent is it to have too many deletion types?
    if ( i>=NI ) return;    // too many choices; discard

    if ( !dfrq->len[i] ) dfrq->len[i] = len;    // new deletion
    dfrq->freq[i]++;
}

read_cns_t *rcns_init(hts_pos_t pos, hts_pos_t beg, hts_pos_t end)
{
    read_cns_t *rcns = (read_cns_t*) calloc(1,sizeof(read_cns_t));
    rcns->pos  = pos;
    rcns->beg  = beg;
    rcns->end  = end;
    rcns->cv_heap = khp_init(cvh);
    if ( init_arrays(rcns)!=0 )
    {
        rcns_destroy(rcns);
        return NULL;
    }
    return rcns;
}

int rcns_set_reads(read_cns_t *rcns, bam_pileup1_t *plp, int nplp)
{
    // save the reads for phasing, this can be called multiple times
    rcns->plp  = plp;
    rcns->nplp = nplp;

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
                for (j=j_beg; j<=j_end; j++, x++, y++) add_base(rcns,x,bam_seqi(seq,y));
            }
            else if ( op==BAM_CINS )
            {
                if ( x>=rcns->beg && x<rcns->end )
                {
                    local_band += p->indel;
                    add_ins(rcns,x-1,y,seq,len);    // x-1: one base before as in VCF
                    y += len;
                }
            }
            else if ( op==BAM_CDEL )
            {
                if ( x>=rcns->beg && x+len-1<=rcns->end )
                {
                    local_band += -p->indel;
                    add_del(rcns,x-1,len);          // x-1: one base before as in VCF
                    x += len;
                }
            }
            if ( local_band_max < local_band ) local_band_max = local_band;
        }

        // Track the biggest deviation +/- from diagonal, used in BAQ alignment step.
        if ( rcns->band < local_band_max ) rcns->band = local_band_max;
    }

    return 0;
}

#if DEBUG_RCNS
static void debug_print_base_freqs(read_cns_t *rcns, const char *ref)
{
    int i,j,k,n = rcns->end - rcns->beg + 1;
    base_freq_t *bfreq = rcns->base_freq;
    ins_freq_t *ifreq  = rcns->ins_freq;
    del_freq_t *dfreq  = rcns->del_freq;
    for (i=0; i<n && ref[i]; i++)
    {
        fprintf(stderr,"%"PRIhts_pos" %c\t",rcns->beg+i+1,ref[i]);
        for (j=0; j<6; j++)
            fprintf(stderr,"\t%d%s",bfreq[i].base[j],ref[i]=="ACGTNi"[j]?"*":"");
        fprintf(stderr,"\t");
        for (j=0; j<NI && dfreq[i].len[j]; j++)
            fprintf(stderr," -%d:%d",dfreq[i].len[j],dfreq[i].freq[j]);
        fprintf(stderr,"\t");
        for (j=0; j<NI && ifreq[i].len[j]; j++)
        {
            fprintf(stderr," +");
            for (k=0; k<ifreq[i].len[j]; k++) fprintf(stderr,"%c",seq_nt16_str[(int)ifreq[i].nt16_seq[j][k]]);
            fprintf(stderr,":%d",ifreq[i].freq[j]);
        }
        fprintf(stderr,"\n");
    }
}
static void debug_print_candidate_variants(read_cns_t *rcns)
{
    int i;
    fprintf(stderr,"Candidate variants:\n");
    for (i=0; i<rcns->ncvar; i++)
    {
        candidate_var_t *var = &rcns->cvar[i];
        fprintf(stderr,"\tvar%d  pos=%"PRIhts_pos" idx=%d vtype=%d which=%d depth=%d af=%f af_dev=%f\n",
            i,var->pos+1,var->idx,var->vtype,var->which,var->depth,var->af,var->af_dev);
    }
}
static void debug_print_haplotype_frequency_spectrum(read_cns_t *rcns)
{
    int i,j;
    fprintf(stderr,"Haplotype frequencies:\n");
    for (i=0; i<NHAP; i++)
    {
        if ( !rcns->hap_freq[i] ) continue;
        fprintf(stderr,"\t%d: ",i);
        for (j=0; j<rcns->ncvar; j++)
            fprintf(stderr,"%d", i&(1<<j) ? 1 : 0);
        fprintf(stderr,"\t%d\n", rcns->hap_freq[i]);
    }
}
static void debug_print_consensus(read_cns_t *rcns)
{
    int i,j;
    for (i=0; i<2; i++)
    {
        if ( !rcns->cns[i].nseq ) break;
        fprintf(stderr,"Consensus%d: ",i);
        for (j=0; j<rcns->cns[i].ipos; j++)
            fprintf(stderr,"%c","ACGTN"[(int)rcns->cns[i].seq[j]]);
        fprintf(stderr,"#");
        for (; j<rcns->cns[i].nseq; j++)
            fprintf(stderr,"%c","ACGTN"[(int)rcns->cns[i].seq[j]]);
        fprintf(stderr,"\n");
    }
}
#else
#define debug_print_base_freqs(rcns,ref)
#define debug_print_candidate_variants(rcns)
#define debug_print_haplotype_frequency_spectrum(rcns)
#define debug_print_consensus(rcns)
#endif

static int cvar_pos_cmp(const void *aptr, const void *bptr)
{
    candidate_var_t *a = (candidate_var_t*)aptr;
    candidate_var_t *b = (candidate_var_t*)bptr;
    if ( a->pos < b->pos ) return -1;
    if ( a->pos > b->pos ) return 1;
    if ( a->vtype < b->vtype ) return -1;
    if ( a->vtype > b->vtype ) return 1;
    if ( a->which < b->which ) return -1;
    if ( a->which > b->which ) return 1;
    return 0;
}
static void register_variant(read_cns_t *rcns, enum variant_type vtype, int cns_pos, int which, int depth, float freq)
{
    cvar_heap_t *cv_heap = rcns->cv_heap;
    if ( vtype==done )
    {
        rcns->ncvar = 0;
        while (cv_heap->ndat)
        {
            rcns->cvar[rcns->ncvar++] = cv_heap->dat[0];
            khp_delete(cvh,cv_heap);
        }
        // sort the variants by pos,type,which to make determination of haplotypes from reads faster
        if ( rcns->ncvar )
            qsort(rcns->cvar, rcns->ncvar, sizeof(*rcns->cvar), cvar_pos_cmp);
        return;
    }

    candidate_var_t var;
    var.pos    = cns_pos + rcns->beg;
    var.which  = which;
    var.vtype  = vtype;
    var.depth  = depth;
    var.af_dev = fabs(0.5-freq);
    var.af     = freq;

    int free_slot;

    // keep the number of variants small, maximum MAX_NCVAR
    if ( rcns->ncvar==MAX_NCVAR )
    {
        if ( cvar_not_preferred(&var,&cv_heap->dat[0]) ) return;  // no need to add, the new variant is worse than the heap's worst one
        free_slot = cv_heap->dat[0].idx;
        khp_delete(cvh,cv_heap);
    }
    else
        free_slot = rcns->ncvar++;
    var.idx =  free_slot;
    rcns->cvar[free_slot] = var;
    khp_insert(cvh,cv_heap,&var);
}

// Identify candidate variant positions. (Note that homozygous variants are not considered
// as those will be added trivially by taking the consensus base.) The detection limit is
// for now hard-wired. This has only indirect effect on sensitivity, will just not contribute
// to the consensus template when realigning.
static int select_candidate_variants(read_cns_t *rcns, const char *ref)
{
    //const float af_th = 0.1;
const float af_th = 0.05;// just for debugging
    int i,j, n = rcns->end - rcns->beg + 1;
    int max_ins_len = 0;    // maximum total length of all insertions applied to allocate big enough buffers
    base_freq_t *bfreq = rcns->base_freq;
    ins_freq_t *ifreq  = rcns->ins_freq;
    del_freq_t *dfreq  = rcns->del_freq;
    for (i=0; i<n && ref[i]; i++)
    {
        for (j=0; j<NI && ifreq[i].len[j]; j++) max_ins_len += ifreq[i].len[j];

        if ( i==rcns->pos - rcns->beg ) continue;   // creating consensus from everything but the variants at the current position

        int dp = 0;
        for (j=0; j<4; j++) dp += bfreq[i].base[j];
        for (j=0; j<NI && dfreq[i].len[j]; j++) dp += dfreq[i].freq[j];
        for (j=0; j<NI && ifreq[i].len[j]; j++) dp += ifreq[i].freq[j];
        float af = 0;   // allele frequency
        for (j=0; j<4; j++)
        {
            if ( !bfreq[i].base[j] || ref[i]=="ACGTN"[j] ) continue;   // ref base or no coverage
            af = (float)bfreq[i].base[j]/dp;
            if ( af>af_th && af<(1-af_th) ) register_variant(rcns,snv,i,j,dp,af);
        }
        for (j=0; j<NI && dfreq[i].len[j]; j++)
        {
            af = (float)dfreq[i].freq[j]/dp;
            if ( af>af_th && af<(1-af_th) ) register_variant(rcns,del,i,j,dp,af);
        }
        for (j=0; j<NI && ifreq[i].len[j]; j++)
        {
            af = (float)ifreq[i].freq[j]/dp;
            if ( af>af_th && af<(1-af_th) ) register_variant(rcns,ins,i,j,dp,af);
        }
    }
    register_variant(rcns,done,0,0,0,0);    // finalize

    // Reallocate buffers
    if ( rcns->mcns < n + max_ins_len )
    {
        n += max_ins_len;
        for (i=0; i<2; i++)
        {
            char *seq = (char*) realloc(rcns->cns[i].seq,sizeof(char)*n);
            if ( !seq ) return -1;
            rcns->cns[i].seq = seq;

            hts_pos_t *pos = (hts_pos_t*) realloc(rcns->cns[i].pos,sizeof(hts_pos_t)*n);
            if ( !pos ) return -1;
            rcns->cns[i].pos = pos;
        }
        rcns->mcns = n;
    }

    // Find the longest deletion at the query position
    i = rcns->pos - rcns->beg;
    rcns->max_del = 0;
    for (j=0; j<NI && j<dfreq[i].len[j]; j++)
    {
        if ( rcns->max_del < -dfreq[i].len[j] ) rcns->max_del = -dfreq[i].len[j];
    }

    return 0;
}
static int create_haplotype_frequency_spectrum(read_cns_t *rcns)
{
    memset(rcns->hap_freq,0,sizeof(rcns->hap_freq));

    int i;
    for (i=0; i<rcns->nplp; i++) // for each read...
    {
        const bam_pileup1_t *p = rcns->plp + i;
        cigar_state_t cigar;
        cstate_init(&cigar,p->b);

        int j,k,hap = 0;
        for (j=0; j<rcns->ncvar; j++)
        {
            candidate_var_t *cvar = &rcns->cvar[j];
            if ( cvar->vtype==snv )
            {
                int iseq = cstate_seek_fwd(&cigar, cvar->pos, BAM_CMATCH, NULL);
                if ( iseq==-2 ) break;
                if ( iseq==-1 ) continue;
                int nt16 = bam_seqi(cigar.seq, iseq);
                if ( seq_nt16_int[nt16]==cvar->which ) hap |= 1<<j;
            }
            else if ( cvar->vtype==ins )
            {
                int len;
                ins_freq_t *ifrq = &rcns->ins_freq[cvar->pos - rcns->beg];
                int iseq = cstate_seek_fwd(&cigar, cvar->pos, BAM_CINS, &len);
                if ( iseq==-2 ) break;
                if ( iseq==-1 ) continue;
                if ( len!=ifrq->len[cvar->which] ) continue;
                for (k=0; k<ifrq->len[cvar->which]; k++)
                    if ( bam_seqi(cigar.seq,iseq+k)!=ifrq->nt16_seq[cvar->which][k] ) break;
                if ( k==ifrq->len[cvar->which] ) hap |= 1<<j;
            }
            else if ( cvar->vtype==del )
            {
                int len;
                del_freq_t *dfrq = &rcns->del_freq[cvar->pos - rcns->beg];
                int ret = cstate_seek_fwd(&cigar, cvar->pos, BAM_CDEL, &len);
                if ( ret==-2 ) break;
                if ( ret==-1 ) continue;
                if ( len!=dfrq->len[cvar->which] ) continue;
                hap |= 1<<j;
            }
        }
        rcns->hap_freq[hap]++;
    }
    return 0;
}

typedef struct
{
    int idx, count;
}
ii_t;

static int ii_cmp(const void *a, const void *b)
{
    if ( ((ii_t*)a)->count > ((ii_t*)b)->count ) return -1;
    if ( ((ii_t*)a)->count < ((ii_t*)b)->count ) return 1;
    return 0;
}

static int correct_haplotype_errors(read_cns_t *rcns)
{
    int i,j, tot = 0;
    ii_t freq[NHAP];
    for (i=0; i<NHAP; i++)
    {
        freq[i].idx = i;
        freq[i].count = rcns->hap_freq[i];
        tot += rcns->hap_freq[i];
    }
    qsort(freq, NHAP, sizeof(ii_t), ii_cmp);
    for (i=NHAP-1; i>=0; i--)
    {
       //if ( freq[1].count > tot - freq[0].count ) break;   // the top2 hapotypes cannot change anymore
        if ( !freq[i].count ) continue;

        // Find a similar haplotype with the highest frequency. Assuming errors go in 0->1
        // direction only and considering one error only.
        int count = freq[i].count, max_hap = 0;
        for (j=0; j<MAX_NCVAR; j++)
        {
            if ( !(freq[i].idx & (1U<<j)) ) continue;   // j-th bit not set in this haplotype
            int hap = freq[i].idx ^ (1U<<j);
            assert( hap>=0 && hap<NHAP );
            if ( count < rcns->hap_freq[hap] ) count = rcns->hap_freq[hap], max_hap = hap;
        }
        if ( count == freq[i].count ) continue;

        // Update frequency and sort the two modified elements
        count = freq[i].count;
        freq[i].count = 0;
        rcns->hap_freq[freq[i].idx] = 0;
        rcns->hap_freq[max_hap] += count;
        for (j=i+1; j<NHAP; j++)
        {
            if ( !freq[j].count ) break;
            ii_t tmp = freq[j-1]; freq[j-1] = freq[j]; freq[j] = tmp;
        }
        for (j=i-1; j>=0; j--)
        {
            if ( freq[j].idx==max_hap ) freq[j].count += count;   // update the best matching haplotype
            if ( freq[j].count < freq[j+1].count )
            {
                ii_t tmp = freq[j]; freq[j] = freq[j+1]; freq[j+1] = tmp;
            }
        }
    }
    rcns->cns_hap[0] = freq[0].idx;
    rcns->cns_hap[1] = freq[1].idx;

    // Use only one consensus if the next best haplotype is populated by less than 10% of reads
    rcns->ncns = (freq[1].count / (freq[0].count + freq[1].count) < 0.1) ? 1 : 2;
if (freq[1].count) rcns->ncns = 2;

    return 0;
}

void create_consensus(read_cns_t *rcns, const char *ref, int ith)
{
    int n = rcns->end - rcns->beg + 1;
    cns_seq_t *cns = &rcns->cns[ith];
    base_freq_t *bfreq = rcns->base_freq;
    ins_freq_t *ifreq  = rcns->ins_freq;
    del_freq_t *dfreq  = rcns->del_freq;
    int j,k, ivar = 0;
    for (j=0; j<n; j++)
    {
        hts_pos_t ref_pos = rcns->beg + j;
        if ( rcns->pos == ref_pos ) cns->ipos = j;

        while ( ivar < rcns->ncvar && rcns->cvar[ivar].pos < ref_pos ) ivar++;

        if ( ivar >= rcns->ncvar || rcns->cvar[ivar].pos != ref_pos )
        {
            // This position is not recognised as a het variant so take the most frequent base, including
            // a deletion if that is most frequent. However, for deleted bases make sure they are not part
            // of the deletion that is being tested at this positions
            int max_freq = 0, kmax = seq_nt16_int[seq_nt16_table[(int)ref[ref_pos]]];
            int nk = ( ref_pos < rcns->pos || ref_pos > rcns->pos + rcns->max_del ) ? 6 : 5;
            for (k=0; k<nk; k++)
                if ( max_freq < bfreq[j].base[k] ) max_freq = bfreq[j].base[k], kmax = k;

            if ( kmax!=5 )  // the most frequent base can be a deletion
            {
                cns->pos[cns->nseq] = ref_pos;
                cns->seq[cns->nseq++] = kmax;
            }
            if ( rcns->pos == ref_pos ) continue;   // do not apply insertions that are being tested

            // Also check how frequent are insertions adjacent to this position. Note that reads with
            // an insertion usually increment also bfreq counts at this position, but not necessarily so,
            // therefore the counts are approximate
            int nreads = 0;
            for (k=0; k<5; k++) nreads += bfreq[j].base[k];
            max_freq = 0, kmax = 0;
            for (k=0; k<NI && ifreq[j].len[k]; k++)
                if ( max_freq < ifreq[j].freq[k] ) max_freq = ifreq[j].freq[k], kmax = k;

            if ( nreads > max_freq*2 ) continue;  // the most frequent insertion is less than half of the reads

            int len = ifreq[j].len[kmax];
            char *seq = ifreq[j].nt16_seq[kmax];
            for (k=0; k<len; k++)
            {
                cns->pos[cns->nseq] = ref_pos;
                cns->seq[cns->nseq++] = seq_nt16_int[(int)seq[k]];
            }
            continue;
        }
        if ( !(rcns->cns_hap[ith] & (1U<<ivar)) )
        {
            // This position has a heterozygous variant but not in this haplotype. Take the
            // most frequent base different from the ivar-th variant
            int max_freq = 0, kmax = seq_nt16_int[seq_nt16_table[(int)ref[ref_pos]]];
            for (k=0; k<6; k++)
            {
                if ( rcns->cvar[ivar].vtype==snv && rcns->cvar[ivar].which==k ) continue;
                if ( max_freq < bfreq[j].base[k] ) max_freq = bfreq[j].base[k], kmax = k;
            }
            if ( kmax!=5 && (!cns->nseq || cns->pos[cns->nseq-1] != ref_pos) )
            {
                cns->pos[cns->nseq] = ref_pos;
                cns->seq[cns->nseq++] = kmax;
            }
            continue;
        }
        int which = rcns->cvar[ivar].which;
        if ( rcns->cvar[ivar].vtype == snv )
        {
            cns->pos[cns->nseq] = ref_pos;
            cns->seq[cns->nseq++] = which;
            continue;
        }

        // There be multiple variants at this position, for example snv+ins. SNVs come first
        // thanks to cvar_pos_cmp(), make sure the base has not been added already.
        if ( !cns->nseq || cns->pos[cns->nseq-1] != ref_pos )
        {
            int max_freq = 0, kmax = seq_nt16_int[seq_nt16_table[(int)ref[ref_pos]]];
            for (k=0; k<6; k++)
            {
                if ( rcns->cvar[ivar].vtype==snv && rcns->cvar[ivar].which==k ) continue;
                if ( max_freq < bfreq[j].base[k] ) max_freq = bfreq[j].base[k], kmax = k;
            }
            if ( kmax!=5 )
            {
                cns->pos[cns->nseq] = ref_pos;
                cns->seq[cns->nseq++] = kmax;
            }
        }
        if ( rcns->cvar[ivar].vtype == ins )
        {
            int len = ifreq[j].len[which];
            char *seq = ifreq[j].nt16_seq[which];
            for (k=0; k<len; k++)
            {
                cns->pos[cns->nseq] = ref_pos;
                cns->seq[cns->nseq++] = seq_nt16_int[(int)seq[k]];
            }
        }
        else if ( rcns->cvar[ivar].vtype == del ) j += dfreq[j].len[which];
    }
}

// The algorithm:
//  1. Identify heterozygous variant positions
//  2. Sort variants by abs(variant_allele_freq-0.5) in descending order
//  3. Take the top sorted variants (up to 8 to fit in uint8_t) and count the number of
//      corresponding reads to create frequency spectrum
//  4. Correct errors, collapse to the requested number of haplotypes (consensus sequences)
//      using majority vote for the distribution tail
cns_seq_t *rcns_get_consensus(read_cns_t *rcns, const char *ref)
{
    debug_print_base_freqs(rcns, ref);

    select_candidate_variants(rcns, ref);
    debug_print_candidate_variants(rcns);

    if ( rcns->ncvar )
    {
        create_haplotype_frequency_spectrum(rcns);
        debug_print_haplotype_frequency_spectrum(rcns);

        correct_haplotype_errors(rcns);
        debug_print_haplotype_frequency_spectrum(rcns);
    }
    else
    {
        rcns->cns_hap[0] = 0;
        rcns->ncns = 1;
    }

    // create consensus
    int i;
    for (i=0; i<rcns->ncns; i++) create_consensus(rcns,ref,i);
    debug_print_consensus(rcns);

    return rcns->cns;
}

