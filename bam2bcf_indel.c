/*

TODO:

- Reevaluate the two STR indel-size adjusting modes.
  Maybe no longer relevant

- Explore INS_PLUS_BASE again. Prefer to disable this as it's hard to
  understand and doesn't work properly on reads ending on an
  insertion.

- Explore indelQ and the effect of STR at boundaries.  I'm not
  convined our quality calculation is correct.  Certainly QUAL appears
  to have little reality with actual indel likelihood!
  It's already there - see end of bcf_cgp_align_score.
  However try tweaking this now we've got better consensus.

- Consider limiting fract to never add more than current depth, so we
  change cons to Ns but not to another base type entirely.

- Set BAQ band width based on maximum size of ins / del observed. Do
  from *all* types, as we may realign from one type to another.

- Trim left/right down better, as we used to.  Judge this based on
  summation of various types and their consensii?

- Consider a separate rfract for lift-over of SNPs than for indels.
  SNPs is good at replacing bases with N where we're unsure on the
  data.  However ref_ins may cause issues with sizing?
  rfract*.8 is working better (so far).  Trying 0.5 too.

- Left-align indels before consensus generation.  Eg:

      /pos being studied
  AGCTGGGGGGAATCG  REF
  AGCT-GGGGGAATCG  Seq type -1
  ACGTGGGGG-AATGCG Seq type 0
      ^

  Type 0 cons shouldn't include the right hand del, but it's outside
  of "biggest_del" window.  Expand this to STR size or left-align.
*/

/*  bam2bcf_indel.c -- indel caller.

    Copyright (C) 2010, 2011 Broad Institute.
    Copyright (C) 2012-2014,2016-2017, 2021 Genome Research Ltd.

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

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/khash_str2int.h>
#include "bam2bcf.h"
#include "str_finder.h"

#include <htslib/ksort.h>
KSORT_INIT_GENERIC(uint32_t)

#define MINUS_CONST 0x10000000

#define MAX_TYPES 64

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef ABS
#  define ABS(a) ((a)<0?-(a):(a))
#endif

#ifndef MAX
#  define MAX(a,b) ((a)>(b)?(a):(b))
#endif

// Take a reference position tpos and convert to a query position (returned).
// This uses the CIGAR string plus alignment c->pos to do the mapping.
//
// *_tpos is returned as tpos if query overlaps tpos, but for deletions
// it'll be either the start (is_left) or end (!is_left) ref position.
static int tpos2qpos(const bam1_core_t *c, const uint32_t *cigar, int32_t tpos, int is_left, int32_t *_tpos)
{
    // x = pos in ref, y = pos in query seq
    int k, x = c->pos, y = 0, last_y = 0;
    *_tpos = c->pos;
    for (k = 0; k < c->n_cigar; ++k) {
        int op = cigar[k] & BAM_CIGAR_MASK;
        int l = cigar[k] >> BAM_CIGAR_SHIFT;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (c->pos > tpos) return y;
            if (x + l > tpos) {
                *_tpos = tpos;
                return y + (tpos - x);
            }
            x += l; y += l;
            last_y = y;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += l;
        else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            if (x + l > tpos) {
                *_tpos = is_left? x : x + l;
                return y;
            }
            x += l;
        }
    }
    *_tpos = x;
    return last_y;
}

// FIXME: check if the inserted sequence is consistent with the homopolymer run
// l is the relative gap length and l_run is the length of the homopolymer on the reference
static inline int est_seqQ(const bcf_callaux_t *bca, int l, int l_run)
{
    int q, qh;
    q = bca->openQ + bca->extQ * (abs(l) - 1);
    qh = l_run >= 3? (int)(bca->tandemQ * (double)abs(l) / l_run + .499) : 1000;
    return q < qh? q : qh;
}

static inline int est_indelreg(int pos, const char *ref, int l, char *ins4)
{
    int i, j, max = 0, max_i = pos, score = 0;
    l = abs(l);
    for (i = pos + 1, j = 0; ref[i]; ++i, ++j) {
        if (ins4) score += (toupper(ref[i]) != "ACGTN"[(int)ins4[j%l]])? -10 : 1;
        else score += (toupper(ref[i]) != toupper(ref[pos+1+j%l]))? -10 : 1;
        if (score < 0) break;
        if (max < score) max = score, max_i = i;
    }
    return max_i - pos;
}

// Identify spft-clip length, position in seq, and clipped seq len
static inline void get_pos(const bcf_callaux_t *bca, bam_pileup1_t *p,
                           int *sc_len_r, int *slen_r, int *epos_r, int *end) {
    bam1_t *b = p->b;
    int sc_len = 0, sc_dist = -1, at_left = 1;
    int epos = p->qpos, slen = b->core.l_qseq;
    int k;
    uint32_t *cigar = bam_get_cigar(b);
    *end = -1;
    for (k = 0; k < b->core.n_cigar; k++) {
        int op = bam_cigar_op(cigar[k]);
        if (op == BAM_CSOFT_CLIP) {
            slen -= bam_cigar_oplen(cigar[k]);
            if (at_left) {
                // left end
                sc_len += bam_cigar_oplen(cigar[k]);
                epos -= sc_len; // don't count SC in seq pos
                sc_dist = epos;
                *end = 0;
            } else {
                // right end
                int srlen = bam_cigar_oplen(cigar[k]);
                int rd = b->core.l_qseq - srlen - p->qpos;
                if (sc_dist < 0 || sc_dist > rd) {
                    // closer to right end than left
                    // FIXME: compensate for indel length too?
                    sc_dist = rd;
                    sc_len = srlen;
                    *end = 1;
                }
            }
        } else if (op != BAM_CHARD_CLIP) {
            at_left = 0;
        }
    }

    if (p->indel > 0 && slen - (epos+p->indel) < epos)
        epos += p->indel-1; // end of insertion, if near end of seq

    // slen is now length of sequence minus soft-clips and
    // epos is position of indel in seq minus left-clip.
    *epos_r = (double)epos / (slen+1) * bca->npos;

    if (sc_len) {
        // scale importance of clip by distance to closest end
        *sc_len_r = 15.0*sc_len / (sc_dist+1);
        if (*sc_len_r > 99) *sc_len_r = 99;
    } else {
        *sc_len_r = 0;
    }

    *slen_r = slen;
}

// Part of bcf_call_gap_prep.
//
// Scans the pileup to identify all the different sizes of indels
// present.
//
// Returns types and fills out n_types_r,  max_rd_len_r and ref_type_r,
//         or NULL on error.
static int *bcf_cgp_find_types(int n, int *n_plp, bam_pileup1_t **plp,
                               int pos, bcf_callaux_t *bca, const char *ref,
                               int *max_rd_len_r, int *n_types_r,
                               int *ref_type_r, int *N_r) {
    int i, j, t, s, N, m, max_rd_len, n_types;
    int n_alt = 0, n_tot = 0, indel_support_ok = 0;
    uint32_t *aux;
    int *types;

    // N is the total number of reads
    for (s = N = 0; s < n; ++s)
        N += n_plp[s];

    bca->max_support = bca->max_frac = 0;
    aux = (uint32_t*) calloc(N + 1, 4);
    if (!aux)
        return NULL;

    m = max_rd_len = 0;
    aux[m++] = MINUS_CONST; // zero indel is always a type (REF)

    // Fill out aux[] array with all the non-zero indel sizes.
    // Also tally number with indels (n_alt) and total (n_tot).
    for (s = 0; s < n; ++s) {
        int na = 0, nt = 0;
        for (i = 0; i < n_plp[s]; ++i) {
            const bam_pileup1_t *p = plp[s] + i;
            ++nt;
            if (p->indel != 0) {
                ++na;
                aux[m++] = MINUS_CONST + p->indel;
            }

            // FIXME: cache me in pileup struct.
            j = bam_cigar2qlen(p->b->core.n_cigar, bam_get_cigar(p->b));
            if (j > max_rd_len) max_rd_len = j;
        }
        double frac = (double)na/nt;
        if ( !indel_support_ok && na >= bca->min_support
             && frac >= bca->min_frac )
            indel_support_ok = 1;
        if ( na > bca->max_support && frac > 0 )
            bca->max_support = na, bca->max_frac = frac;

        n_alt += na;
        n_tot += nt;
    }

    // Sort aux[] and dedup
    ks_introsort(uint32_t, m, aux);
    for (i = 1, n_types = 1; i < m; ++i)
        if (aux[i] != aux[i-1]) ++n_types;

    // Taking totals makes it hard to call rare indels (IMF filter)
    if ( !bca->per_sample_flt )
        indel_support_ok = ( (double)n_alt / n_tot < bca->min_frac
                             || n_alt < bca->min_support )
            ? 0 : 1;
    if ( n_types == 1 || !indel_support_ok ) { // then skip
        free(aux);
        return NULL;
    }

    // Bail out if we have far too many types of indel
    if (n_types >= MAX_TYPES) {
        free(aux);
        // TODO revisit how/whether to control printing this warning
        if (hts_verbose >= 2)
            fprintf(stderr, "[%s] excessive INDEL alleles at position %d. "
                    "Skip the position.\n", __func__, pos + 1);
        return NULL;
    }

    // To prevent long stretches of N's to be mistaken for indels
    // (sometimes thousands of bases), check the number of N's in the
    // sequence and skip places where half or more reference bases are Ns.
    int nN=0, i_end = pos + (2*bca->indel_win_size < max_rd_len
                            ?2*bca->indel_win_size : max_rd_len);
    for (i=pos; i<i_end && ref[i]; i++)
        nN += ref[i] == 'N';
    if ( nN*2>(i-pos) ) {
        free(aux);
        return NULL;
    }

    // Finally fill out the types[] array detailing the size of insertion
    // or deletion.
    types = (int*)calloc(n_types, sizeof(int));
    if (!types) {
        free(aux);
        return NULL;
    }
    t = 0;
    for (i = 0; i < m; ++i) {
        int sz = (int32_t)(aux[i] - MINUS_CONST);
        int j;
        for (j = i+1; j < m; j++)
            if (aux[j] != aux[i])
                break;

        if (sz == 0
            || j-i >= bca->min_support
            // Note, doesn't handle bca->per_sample_flt yet
            || bca->per_sample_flt 
            || (double)(j-i) / n_tot >= bca->min_frac)
            types[t++] = sz;
        i = j-1;
    }
    free(aux);

    if (t <= 1)
        return NULL;
    n_types = t;

    // Find reference type; types[?] == 0)
    for (t = 0; t < n_types; ++t)
        if (types[t] == 0) break;

    *ref_type_r   = t;
    *n_types_r    = n_types;
    *max_rd_len_r = max_rd_len;
    *N_r          = N;

    return types;
}

// Increment ins["str"] and freq["str"]
#define NI 10 // number of alternative insertion sequences
// Could use a hash table too, but expectation is a tiny number of alternatives
typedef struct {
    char *str[NI];
    int len[NI];
    int freq[NI];
} str_freq;

static int bcf_cgp_append_cons(str_freq *sf, char *str, int len, int freq) {
    int j;

    for (j = 0; j < NI && sf->str[j]; j++) {
        if (sf->len[j] == len && memcmp(sf->str[j], str, len) == 0)
            break;
    }
    if (j >= NI)
        return 0; // too many choices; discard

    sf->freq[j]+=freq;
    if (!sf->str[j]) {
        // new insertion
        if (!(sf->str[j] = malloc(len+1)))
            return -1;
        memcpy(sf->str[j], str, len);
        sf->len[j] = len;
    }

    return 0;
}

/*
 * Compute the consensus for a specific indel type at pos.
 *
 * left_shift is the number of inserted(+) or deleted(-) bases added to
 * the consensus before we get to pos.  This is necessary so the alignment
 * band is correct as it's expected to start at left/right edges in
 * sync
 *
 * We accumulate into several buffers for counting base types:
 * cons_base   - consensus of data with p->indel == type, bases or gap
 * ref_base    - consensus of data with p->indel != type, bases or gap
 * cons_ins    - consensus of data with p->indel == type, insertions
 * ref_ins     - consensus of data with p->indel == type, bases or gap
 *
 * The purpose of cons_ins vs cons_base is if we have very low
 * coverage due to nearly all reads being another type, then we can
 * still get a robust consensus using the other data.  If we don't
 * have shallow data, then we'll not use as much of ref_base as we may
 * have correlated variants.
 *
 * Eg:
 * REF: AGCTATGAGGCTGATA
 * SEQ: AGGTAGGAGGGTGATA (x1)
 * SEQ: AGCTACGAGG*TGATA (x24)
 * SEQ: AGCTACTAGG*TGATA (x24)
 *
 * Cons for no-del is Cs not Gs.  Cannot trust it, so use N if shallow.
 * CON: AGCTACNAGGGTGATA
 */
static char **bcf_cgp_consensus(int n, int *n_plp, bam_pileup1_t **plp,
                                int pos, bcf_callaux_t *bca, const char *ref,
                                int left, int right,
                                int sample, int type, int biggest_del,
                                int *left_shift, int *right_shift) {
    int (*cons_base)[6] = calloc(right - left + 1, sizeof(*cons_base));// single base or del
    str_freq *cons_ins  = calloc(right - left + 1, sizeof(*cons_ins)); // multi-base insertions

    // non-indel ref for all reads on this sample, rather than those just
    // matching type.  We use this for handling the case where we have a
    // homozygous deletion being studied, but with 1 or 2 reads misaligned
    // and containing a base there.
    //
    // Eg if the type[]=0 consensus is made up of a very small sample size,
    // which is also enriched for highly error prone data.  We can use
    // the other reads from type[] != 0 to flesh out the consensus and
    // improve accuracy.
    int (*ref_base)[6]  = calloc(right - left + 1, sizeof(*ref_base));
    str_freq *ref_ins   = calloc(right - left + 1, sizeof(*ref_ins));

    int i, j, k, s = sample;

// cons_ins sequence is the insertion seq followed by the
// next match base
#define INS_PLUS_BASE


    // Accumulate sequences into cons_base and cons_ins arrays
    int last_base_ins = 0;
    for (i = 0; i < n_plp[s]; i++) {
        const bam_pileup1_t *p = plp[s] + i;
//        if (p->indel != type)
//            continue;

        //        fprintf(stderr, "p=%d\t%d/%d: Seq %3d of %3d\t", p->b->core.pos, s, type, i, n_plp[s]);

        bam1_t *b = p->b;
        int x = b->core.pos;  // ref coordinate
        int y = 0;            // seq coordinate
        uint32_t *cigar = bam_get_cigar(b);
        uint8_t *seq = bam_get_seq(b);

        last_base_ins = 0;
        for (k = 0; k < b->core.n_cigar; ++k) {
            int op  = cigar[k] &  BAM_CIGAR_MASK;
            int len = cigar[k] >> BAM_CIGAR_SHIFT;
            int base;

            switch(op) {
            case BAM_CSOFT_CLIP:
                y += len;
                break;

            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF: {
                int L[16] = {
                    // 1,2,4,8 to 0,1,2,3 plus 4 for N/ambig (and 5 for gap)
                    4,0,1,4, 2,4,4,4,  3,4,4,4, 4,4,4,4
                };

                // Can short-cut this with j_start and j_end based on x+len and left,right
                for (j = 0; j < len; j++, x++, y++) {
                    if (x < left) continue;
                    if (x >= right) break;

#ifdef INS_PLUS_BASE
                    // FIXME: need last_base_ins_type and last_base_ins_ref?
                    if (last_base_ins) {
                        last_base_ins = 0;
                        continue;
                    }
#endif
                    base = bam_seqi(seq, y);
#ifdef INS_PLUS_BASE
                    if (p->indel == type)
#else
                    if (p->indel == type || p->indel > 0) // alternative
#endif
                        cons_base[x-left][L[base]]++;
                    else if (x != pos+1) // indel being assessed question
                        ref_base[x-left][L[base]]++;
                    //                    fputc(seq_nt16_str[base], stderr);
                    // else last_base_ins=0?
                }
                break;
            }

            case BAM_CINS: {
//                if (p->indel != type) {
//                    y += len; // for when adding to ref_base
//                    break;
//                }

                char ins[1024];
                for (j = 0; j < len; j++, y++) {
                    if (x < left) continue;
                    if (x >= right) break;
                    base = bam_seqi(seq, y);
                    if (j < 1024)
                        ins[j] = seq_nt16_str[base];
                }

                // Insertions come before a ref match.
                // 5I 5M is IIIIIM M M M M events, not
                // {IIIII,M} M M M M choice.  So we need to include the
                // next match in our sequence when choosing the consensus.
#ifdef INS_PLUS_BASE
                if (y < b->core.l_qseq) {
                    base = bam_seqi(seq, y);
                    if (j < 1024)
                        ins[j++] = seq_nt16_str[base];
                }
                //last_base_ins = 1;
#endif

                //                fprintf(stderr, "<+%.*s>", j<1024?j:1024, ins);
                if (x >= left && x < right) {
                    int ilen = j<1024?j:1024;
                    if (p->indel == type) {
                        bcf_cgp_append_cons(&cons_ins[x-left], ins, ilen, 1);
                    } else  if (x != pos+1){
                        bcf_cgp_append_cons(&ref_ins[x-left],  ins, ilen, 1);
                    }
#ifdef INS_PLUS_BASE
                    last_base_ins = 1;
#endif
                }
                break;
            }

            case BAM_CDEL:
                // FIXME, not perfect for I/D combos, but likely sufficient.
                last_base_ins = 0;
                for (j = 0; j < len; j++, x++) {
                    if (x < left) continue;
                    if (x >= right) break;
                    //                    fputc('-', stderr);
                    if (p->indel == type)
                        // fixme: not p->indel==type but x==pos+1
                        cons_base[x-left][5]++;
                    else
                        ref_base[x-left][5]++;
                }
                break;
            }
        }
        //        fprintf(stderr, " %s\n", bam_get_qname(p->b));
    }

    // Expand cons_base to include depth from ref_base/ref_ins
    // Caveat: except at pos itself, where true ref is used if type != 0
    for (i = 0; i < right-left; i++) {
        // Total observed depth
        int t = cons_base[i][0] + cons_base[i][1] + cons_base[i][2] +
            cons_base[i][3] + cons_base[i][4] + cons_base[i][5];
        for (j = 0; j < NI; j++) {
            if (!cons_ins[i].str[j])
                break;
            t += cons_ins[i].freq[j];
        }
        int r = ref_base[i][0] + ref_base[i][1] + ref_base[i][2] +
            ref_base[i][3] + ref_base[i][4] + ref_base[i][5];
        for (j = 0; j < NI; j++) {
            if (!ref_ins[i].str[j])
                break;
            r += ref_ins[i].freq[j];
        }

        // When evaluating this particular indel, we don't want to
        // penalise alignments by SNP errors elsewhere.  This can
        // happen when we have low depth for a particular 'type'.
        // 
        // So add in a little data from ref_base/ref_ins.
        double rfract = (r - t*2)*.75 / (r+1);

//        // We ensure this is at least 1 fold deep, and we try to add
//        // no more than the amount of coverage in this consesnsus.
//        double rfract = (MIN(r, t*3.333+1) - t*2)*.75 / (r+1);

        //rfract*=.5; // -FN +FP/GT.  Which poison do we want?
        if (rfract < 1.01 / (r+1e-10))
            rfract = 1.01 / (r+1e-10); // low depth compensation

        // TODO: consider limiting rfract so we never drown out the
        // signal.  We want to use the remaining data only to correct
        // for sequencing errors in low depth alleles.  If we get
        // conflicts, it's better to use N than to change a base
        // incase that variant is genuine.

        if (1 || rfract > 0) { //  && !(type == 0 && i+left == pos)) {
            if (i+left >= pos+1 && i+left < pos+1-biggest_del) {
                // We're overlapping the current indel region, so
                // we don't wish to bring in evidence from the other
                // "type" data.
                continue; // 10f+
            } else {
                // Otherwise add in a portion of other data to
                // boost low population numbers.
                cons_base[i][0] += rfract * ref_base[i][0];
                cons_base[i][1] += rfract * ref_base[i][1];
                cons_base[i][2] += rfract * ref_base[i][2];
                cons_base[i][3] += rfract * ref_base[i][3];
                cons_base[i][4] += rfract * ref_base[i][4];
                cons_base[i][5] += rfract * ref_base[i][5];
            }

            for (j = 0; j < NI; j++) {
                if (!ref_ins[i].str[j])
                    break;
                bcf_cgp_append_cons(&cons_ins[i],
                                    ref_ins[i].str[j], ref_ins[i].len[j],
                                    rfract * ref_ins[i].freq[j]);
            }
        }
    }

    // Allocate consensus buffer, to worst case length
    int max_len = right-left;
    for (i = 0; i < right-left; i++) {
        if (!cons_ins[i].str[0])
            continue;

        int ins = 0;
        for (j = 0; j < NI; j++) {
            if (!cons_ins[i].str[j])
                break;
            if (cons_ins[i].str[j] && ins < cons_ins[i].len[j])
                ins = cons_ins[i].len[j];
        }
        max_len += ins;
    }
    char **cons = malloc((max_len+1)*2 + sizeof(char *)*2);
    cons[0] = (char *)&cons[2];
    cons[1] = cons[0] + max_len+1;

    // Merge insertions where they are the same length but different
    // sequences.
    // NB: we could just index by length and have accumulators for each,
    // instead of storing separately and merging later (here).
    // Ie str_freq.str is [NI][5] instead.
    for (i = 0; i < right-left; i++) {
        int ins[1024][5];
        for (j = 0; j < NI; j++) {
            if (!cons_ins[i].str[j])
                break;

            if (cons_ins[i].freq[j] == 0)
                continue; // already merged

            int l;
            for (l = 0; l < cons_ins[i].len[j]; l++) {
                // FIXME! optimise this
                ins[l][0] = ins[l][1] = ins[l][2] = ins[l][3] = ins[l][4] = 0;
                switch(cons_ins[i].str[j][l]) {
                case 'A': ins[l][0] = cons_ins[i].freq[j]; break;
                case 'C': ins[l][1] = cons_ins[i].freq[j]; break;
                case 'G': ins[l][2] = cons_ins[i].freq[j]; break;
                case 'T': ins[l][3] = cons_ins[i].freq[j]; break;
                default:  ins[l][4] = cons_ins[i].freq[j]; break;
                }
            }

            // Merge other insertions of the same length to ins[] counters
            for (k = j+1; k < NI; k++) {
                if (!cons_ins[i].str[k])
                    break;
                if (cons_ins[i].len[k] != cons_ins[i].len[j])
                    continue;
                if (cons_ins[i].freq[k] == 0)
                    continue; // redundant?

                // Merge str[j] and str[k]
                for (l = 0; l < cons_ins[i].len[k]; l++) {
                    // FIXME! optimise this
                    switch(cons_ins[i].str[k][l]) {
                    case 'A': ins[l][0]+=cons_ins[i].freq[k]; break;
                    case 'C': ins[l][1]+=cons_ins[i].freq[k]; break;
                    case 'G': ins[l][2]+=cons_ins[i].freq[k]; break;
                    case 'T': ins[l][3]+=cons_ins[i].freq[k]; break;
                    default:  ins[l][4]+=cons_ins[i].freq[k]; break;
                    }
                }
                cons_ins[i].freq[j] += cons_ins[i].freq[k];
                cons_ins[i].freq[k] = 0;
            }

            // Now replace ins[j] with the consensus insertion of this len.
            for (l = 0; l < cons_ins[i].len[j]; l++) {
                int max_v = 0, base = 0;
                int tot = ins[l][0] + ins[l][1] + ins[l][2]
                        + ins[l][3] + ins[l][4];
                if (max_v < ins[l][0]) max_v = ins[l][0], base = 0;
                if (max_v < ins[l][1]) max_v = ins[l][1], base = 1;
                if (max_v < ins[l][2]) max_v = ins[l][2], base = 2;
                if (max_v < ins[l][3]) max_v = ins[l][3], base = 3;
                if (max_v < ins[l][4]) max_v = ins[l][4], base = 4;

                cons_ins[i].str[j][l] = (max_v > 0.6*tot) ?"ACGTN"[base] :'N';
            }
        }
    }

// TODO: try CONS_CUTOFF higher, eg .6, to force more Ns?
#define CONS_CUTOFF      .40 // 40% needed for base vs N
#define CONS_CUTOFF2     .80 // 80% needed for gap in cons[1]
#define CONS_CUTOFF_INC  .40 // 40% to include any insertion cons[0]
#define CONS_CUTOFF_INC2 .80 // 80% to include any insertion cons[1] HOM
#define CONS_CUTOFF_INS  .60 // and then 60% needed for it to be bases vs N
    // Walk through the frequency arrays to call the consensus.
    // We produce cons[0] and cons[1].  Both include strongly
    // homozygous indels.  Both also include the indel at 'pos'.
    // However for heterozygous indels we call the most likely event
    // for cons[0] and the less-likely alternative in cons[1].
    // TODO: a proper phase analysis so multiple events end up
    // combining together into the correct consensus.
    *left_shift = 0;
    *right_shift = 0;
    int cnum;

    // Het call filled out in cnum==0 (+ve or -ve)
    // Used in cnum==1 to do the opposite of whichever way we did before.
    int heti[1024] = {0}, hetd[1024] = {0};

    for (cnum = 0; cnum < 2; cnum++) {
        for (i = k = 0; i < right-left; i++) {
            //        fprintf(stderr, "%d\t", i);
            int max_v = 0, max_v2 = 0, max_j = 4, max_j2 = 4, tot = 0;
            for (j = 0; j < 6; j++) {
                // Top 2 consensus calls
                if (max_v < cons_base[i][j]) {
                    max_v2 = max_v, max_j2 = max_j;
                    max_v = cons_base[i][j], max_j = j;
                } else if (max_v2 < cons_base[i][j]) {
                    max_v2 = cons_base[i][j], max_j2 = j;
                }
                tot += cons_base[i][j];
                //            if (cons_base[i][j])
                //                fprintf(stderr, "%c%d ", "ACGTN*"[j], cons_base[i][j]);
            }

            // +INS
            int max_v_ins = 0, max_j_ins = 0;
            int tot_ins = 0;
            for (j = 0; j < NI; j++) {
                if (!cons_ins[i].str[j])
                    break;
                if (cons_ins[i].freq[j] == 0)
                    continue; // previously merged

                if (max_v_ins < cons_ins[i].freq[j])
                    //if (i != pos-left+1 || cons_ins[i].len[j] == type)
                    max_v_ins = cons_ins[i].freq[j], max_j_ins = j;
                tot_ins += cons_ins[i].freq[j];

                //            fprintf(stderr, "%.*s%d ", cons_ins[i].len[j], cons_ins[i].str[j],
                //                    cons_ins[i].freq[j]);
            }
            // NB: tot is based on next matching base, so it includes
            // everything with or without the insertion.
#ifdef INS_PLUS_BASE
//            if (max_v_ins > CONS_CUTOFF_INC *(tot+tot_ins) && (cnum==0 ||
//                max_v_ins > CONS_CUTOFF_INC2*(tot+tot_ins) ||
//                i == pos-left+1)) {
            int always_ins =
                (i == pos-left+1 && type>0) ||             // current eval
                max_v_ins > CONS_CUTOFF_INC2*(tot+tot_ins);// HOM
            int het_ins = 0;
            if (!always_ins && max_v_ins >= bca->min_support) {
                // Candidate HET ins.
                if (cnum == 0) {
                    het_ins = max_v_ins > CONS_CUTOFF_INC *(tot+tot_ins);
                    if (i < 1024) heti[i] = het_ins
                                      ? 1
                                      : (max_v_ins > .2*(tot+tot_ins) ? -1:0);
                } else {
                    het_ins = (heti[i] == -1); // HET but uncalled before
                }
            }
//            if (max_v_ins)
//                fprintf(stderr, "Cons @ %d: type %d cnum %d always %d het_ins %d // max_v %d vs %d+%d\n", i, type, cnum, always_ins, het_ins, max_v_ins, tot, tot_ins);
            if (always_ins || het_ins) {
//            if ((i == pos-left+1 && type) || // current 'type' at pos
//                max_v_ins > CONS_CUTOFF_INC2*(tot+tot_ins) ||  // HOM
//                (max_v_ins > bca->min_support &&
//                 (cnum != 0) ^ (max_v_ins > CONS_CUTOFF_INC *(tot+tot_ins)))) { // HET
#else
            if ((i == pos-left+1 && type) || // current 'type' at pos
                max_v_ins > CONS_CUTOFF_INC2*tot ||  // HOM
                (max_v_ins > bca->min_support &&
                 (cnum != 0) ^ max_v_ins > CONS_CUTOFF_INC*tot)) { // HET
#endif
                if (max_v_ins > CONS_CUTOFF_INS*tot_ins) {
                    // Insert bases
                    for (j = 0; j < cons_ins[i].len[max_j_ins]; j++) {
                        // FIXME: commented out to deliberate get consensus shift.
                        // Need to know how to get aligner working properly in that
                        // scenario, as it'll happen sometimes!
                        if (cnum == 0) {
                            if (k < pos-left+*left_shift)
                                (*left_shift)++;
                            else
                                (*right_shift)++;
                        }
                        cons[cnum][k++] = cons_ins[i].str[max_j_ins][j];
                    }
                } else {
                    for (j = 0; j < cons_ins[i].len[max_j_ins]; j++)
                        cons[cnum][k++] = 'N';
                }
                // don't call next base as included in insertion
#ifdef INS_PLUS_BASE
                continue;
#endif
            }

            // Call deletions
            int always_del = (type < 0 && i > pos-left && i <= pos-left-type)
                || cons_base[i][5] > CONS_CUTOFF2 * tot; // HOM del
            int het_del = 0;
            if (!always_del && cons_base[i][5] >= bca->min_support) {
                // Candidate HET del.
                if (cnum == 0) {
                    het_del = cons_base[i][5] >= CONS_CUTOFF * tot;
                    if (i < 1024) {
                        if (i >= pos-left && i <= pos-left-biggest_del)
                            hetd[i] = 0;
                        else
                            hetd[i] = het_del
                                ? 1
                                : (cons_base[i][5] >= .2 * tot ? -1 : 0);
                    }
                } else {
                    het_del = (hetd[i] == -1); // HET del uncalled on cnum 0
                    if (max_j == 5 && het_del == 0) {
                        max_v = max_v2;
                        max_j = max_j2;
                    }
                }
            }
            if (always_del || het_del) {
                // Deletion
                if (k < pos-left+*left_shift)
                    (*left_shift)--;
                else
                    (*right_shift)++;
            } else {
                // Finally the easy case - a non-indel base or an N
                if (max_v > CONS_CUTOFF*tot)
                    cons[cnum][k++] = "ACGTN*"[max_j];
                else
                    cons[cnum][k++] = 'N';
            }
        }
        cons[cnum][k++] = '\0';
    }
    //    fprintf(stderr, "Cons:           %s\n", cons);
    free(cons_base);
    free(ref_base);

    for (i = 0; i < right-left; i++) {
        for (j = 0; j < NI; j++) {
            // FIXME: replace by string pool
            if (cons_ins[i].str[j])
                free(cons_ins[i].str[j]);
            if (ref_ins[i].str[j])
                free(ref_ins[i].str[j]);
        }
    }
    free(cons_ins);
    free(ref_ins);

    return cons;
}

// The length of the homopolymer run around the current position
static int bcf_cgp_l_run(const char *ref, int pos) {
    int i, l_run;

    int c = seq_nt16_table[(int)ref[pos + 1]];
    if (c == 15) {
        l_run = 1;
    } else {
        for (i = pos + 2; ref[i]; ++i)
            if (seq_nt16_table[(int)ref[i]] != c) break;
        l_run = i;
        for (i = pos; i >= 0; --i)
            if (seq_nt16_table[(int)ref[i]] != c) break;
        l_run -= i + 1;
    }

    return l_run;
}


// Compute the consensus for this sample 's', minus indels which
// get added later.
static char *bcf_cgp_calc_cons(int n, int *n_plp, bam_pileup1_t **plp,
                               int pos, int *types, int n_types,
                               int max_ins, int s) {
    int i, j, t, k;
    int *inscns_aux = (int*)calloc(5 * n_types * max_ins, sizeof(int));
    if (!inscns_aux)
        return NULL;

    // Count the number of occurrences of each base at each position for
    // each type of insertion.
    for (t = 0; t < n_types; ++t) {
        if (types[t] > 0) {
            for (s = 0; s < n; ++s) {
                for (i = 0; i < n_plp[s]; ++i) {
                    bam_pileup1_t *p = plp[s] + i;
                    if (p->indel == types[t]) {
                        uint8_t *seq = bam_get_seq(p->b);
                        for (k = 1; k <= p->indel; ++k) {
                            int c = seq_nt16_int[bam_seqi(seq, p->qpos + k)];
                            assert(c<5);
                            ++inscns_aux[(t*max_ins+(k-1))*5 + c];
                        }
                    }
                }
            }
        }
    }

    // Use the majority rule to construct the consensus
    char *inscns = (char *)calloc(n_types * max_ins, 1);
    for (t = 0; t < n_types; ++t) {
        for (j = 0; j < types[t]; ++j) {
            int max = 0, max_k = -1, *ia = &inscns_aux[(t*max_ins+j)*5];
            for (k = 0; k < 5; ++k)
                if (ia[k] > max)
                    max = ia[k], max_k = k;
            inscns[t*max_ins + j] = max ? max_k : 4;
            if (max_k == 4) {
                // discard insertions which contain N's
                types[t] = 0;
                break;
            }
        }
    }
    free(inscns_aux);

    return inscns;
}

// Part of bcf_call_gap_prep.
//
// Realign using BAQ to get an alignment score of a single read vs
// a haplotype consensus.  TODO: replace BAQ with something more robust.
//
// There are many coordinates, so let's explain them.
// - left, right, tbeg, tend, r_start and r_end are in aligned reference
//   coordinates.
//   left/right start from pos +/- indel_win_size.
//   r_start/r_end are the BAM first and last mapped coord on the reference.
//   tbeg and tend are the intersection of the two.
// - qbeg and qend are in BAM sequence coordinates
// - qpos is in sequence coordinates, relative to qbeg.
//
// To see what this means, we have illustrations with coordinates
// above the seqs in reference space and below the seqs in BAM seq space.
//
// Overlap left:
//                     tbeg                        tend
//      r_start        left                 pos    r_end          right
// REF  :..............|--------------------#------:--------------|...
// SEQ  :..............|--------------------#------|
//      0              qbeg                 qpos   qend
//
// Overlap right:
//                        r_start                     tend
//         left           tbeg  pos                   right       r_end
// REF  ...|--------------:-----#---------------------|...........:
// SEQ                    |-----#---------------------|...........:
//                        qbeg  qpos                  qend
//                        0
//
// The "-" sequence is the bit passed in.
// Ie ref2 spans left..right and query spans qbeg..qend.
// We need to adjust ref2 therefore to tbeg..tend.
//
// Fills out score
// Returns 0 on success,
//        <0 on error
static int bcf_cgp_align_score(bam_pileup1_t *p, bcf_callaux_t *bca, int type,
                               uint8_t *ref1, uint8_t *ref2, uint8_t *query,
                               int r_start, int r_end, int long_read,
                               int tbeg, int tend1, int tend2,
                               int left, int right,
                               int qbeg, int qend,
                               int qpos, int max_deletion,
                               int *score) {
    // Illumina
    probaln_par_t apf = { 1e-4, 1e-2, 10 };

    // Parameters that work better on PacBio CCS 15k.
    // We should consider querying the header and RG PU field.
    // See also htslib/realn.c:sam_prob_realn()
    if (long_read) {
        apf.d = 1e-3;
        apf.e = 1e-1;
    }

    type = abs(type);
    apf.bw = type + 3; // or abs(l_ref - l_query), so we want to keep similar
    int l, sc1, sc2;
    const uint8_t *qual = bam_get_qual(p->b), *bq;
    uint8_t *qq;

    // Trim poly_Ns at ends of ref.
    // This helps to keep len(ref) and len(query) similar, to reduce
    // band size and reduce the chance of -ve BAQ scores.

    // FIXME  Maybe instead of l>ABS(type) it should be l>query_len/2 ?
    // TODO: no difference to result, but what difference is there to
    // speed? Is this worth it?
#if 1
    for (l = 0; l < tend1-tbeg && l < tend2-tbeg; l++)
        if (ref1[l + tbeg-left] != 4 || ref2[l + tbeg-left] != 4)
            break;
    if (l > ABS(type)) {
        fprintf(stderr, "Prune %d N to left\n", l-ABS(type));
        tbeg += l-ABS(type);
    }

    for (l = tend1-tbeg-1; l >= 0; l--)
        if (ref1[l + tbeg-left] != 4)
            break;
    l = tend1-tbeg-1 - l;
    if (l > ABS(type)) {
        fprintf(stderr, "Prune %d N to right 1\n", l-ABS(type));
        tend1 -= l-ABS(type);
    }

    for (l = tend2-tbeg-1; l >= 0; l--)
        if (ref2[l + tbeg-left] != 4)
            break;
    l = tend2-tbeg-1 - l;
    if (l > ABS(type)) {
        fprintf(stderr, "Prune %d N to right 2\n", l-ABS(type));
        tend2 -= l-ABS(type);
    }
#endif

    // Get segment of quality, either ZQ tag or if absent QUAL.
    if (!(qq = (uint8_t*) calloc(qend - qbeg, 1)))
        return -1;
    bq = (uint8_t*)bam_aux_get(p->b, "ZQ");
    if (bq) ++bq; // skip type
    for (l = qbeg; l < qend; ++l) {
        int qval = bq? qual[l] + (bq[l] - 64) : qual[l];
        if (qval > 30)
            qval = 30;
        if (qval < 7)
            qval = 7;
        qq[l - qbeg] = qval;

        // Skew qq at qpos to be higher than background and qq at
        // other regions to be lower.  This means the alignment of
        // indel we are currently assessing takes precedence over
        // alignment of flanking regions.
        //
        // Ins; type = +ve
        // Ref  AGCTAG---CTGA
        // Qry  AGCTAGGGGCTGA  (qpos..qpos+type)
        //
        // Del; type = -ve
        // Ref  AGCTAGGGGCTGA
        // Qry  AGCTAG---CTGA  (qpos..qpos)

//        // Tests over 1-47MB
//        // shift8b            FP/GT/FN = 290/296/2310
//        // develop                     = 264/326/2282
//        if (l >= qpos-2 && l <= qpos+2+(type>0?type:0))
//            //qq[l-qbeg] += 15;  //qq2 = 282/312/2334
//            qq[l-qbeg] *= 1.5;   //qq3 = 284/305/2326
//            //qq[l-qbeg] *= 0.75;//qq4 = 287/333/2347
////        else
////          qq[l-qbeg] *= 0.67;  // qq = 269/343/2413 (qq3 with else clause)
    }

    // The bottom 8 bits are length-normalised score while
    // the top bits are unnormalised.
    //
    // Try original cons and new cons and pick best.
    // This doesn't removed FN much (infact maybe adds very slightly),
    // but it does reduce GT errors and some slight reduction to FP.
    sc2 = probaln_glocal(ref2 + tbeg - left, tend2 - tbeg,
                         query, qend - qbeg, qq, &apf, 0, 0);

    if (tend1 != tend2 ||
        memcmp((char *)ref1 + tbeg - left, (char *)ref2 + tbeg - left,
               tend1 - tbeg + type) != 0)
        sc1 = probaln_glocal(ref1 + tbeg - left, tend1 - tbeg,
                             query, qend - qbeg, qq, &apf, 0, 0);
    else
        sc1 = INT_MAX; // skip

#if 1
#define CONS_DEBUG
    fprintf(stderr, "\nref1: ");
    for (int j = 0; j < tend1-tbeg; j++)
        putc("ACGTN"[(uint8_t)ref1[j+tbeg-left]], stderr);
    putc('\n', stderr);
    fprintf(stderr, "ref2: ");
    for (int j = 0; j < tend2-tbeg; j++)
        putc("ACGTN"[(uint8_t)ref2[j+tbeg-left]], stderr);
    putc('\n', stderr);
    fprintf(stderr, "qury: ");
    for (int j = 0; j < qend-qbeg; j++)
        putc("ACGTN"[(uint8_t)query[j]], stderr);
    putc('\n', stderr);
    fprintf(stderr, "sc1 %-9d sc2 %-9d ", sc1, sc2);
#endif

    if (sc1 < 0 && sc2 < 0) {
        *score = 0xffffff;
        free(qq);
        return 0;
    }
    if (sc1 < 0) {
        // sc2 is already correct
    } else if (sc2 < 0) {
        sc2 = sc1;
    } else {
        // sc1 and sc2 both pass, so use best
        if (sc2 > sc1)
            sc2 = sc1;
    }

    // used for adjusting indelQ below
    l = (int)((100. * sc2 / (qend - qbeg) + .499) * bca->indel_bias);
    *score = sc2<<8 | MIN(255, l);
    //fprintf(stderr, "score = %d, qend-qbeg = %d, => adj score %d\n", sc, qend-qbeg, l);

    rep_ele *reps, *elt, *tmp;
    uint8_t *seg = ref2 + tbeg - left;
    int seg_len = tend2 - tbeg + type;

    // Note: although seg moves (tbeg varies), ref2 is reused many times
    // so we could factor out some find_STR calls.  However it's not the
    // bottleneck for now.

    // FIXME: need to make this work on IUPAC.
    reps = find_STR((char *)seg, seg_len, 0);
    int iscore = 0;

    // Identify STRs in ref covering the indel up to
    // (or close to) the end of the sequence.
    // Those having an indel and right at the sequence
    // end do not confirm the total length of indel
    // size.  Specifically a *lack* of indel at the
    // end, where we know indels occur in other
    // sequences, is a possible reference bias.
    //
    // This is emphasised further if the sequence ends with
    // soft clipping.
    DL_FOREACH_SAFE(reps, elt, tmp) {
        if (elt->start <= qpos && elt->end >= qpos) {
            iscore += (elt->end-elt->start) / elt->rep_len;  // c
            if (elt->start+tbeg <= r_start ||
                elt->end+tbeg   >= r_end) {
                //iscore += 2*(elt->end-elt->start); //h5  (STR2)
                //iscore += 4*(elt->end-elt->start);  //h5STR4
                iscore += (elt->end-elt->start);  //h5STR1
            }
       }

        DL_DELETE(reps, elt);
        free(elt);
    }

    // Apply STR score to existing indelQ
    l  =  (*score&0xff)*.8 + iscore*2;
    *score = (*score & ~0xff) | MIN(255, l);
#ifdef CONS_DEBUG
    fprintf(stderr, " iscore %-4d  l %d  %d..%d\n\n", iscore, l, r_start, r_end);
#endif

    free(qq);

    return 0;
}

// Part of bcf_call_gap_prep.
//
// Returns n_alt on success
//         -1 on failure
static int bcf_cgp_compute_indelQ(int n, int *n_plp, bam_pileup1_t **plp,
                                  bcf_callaux_t *bca, char *inscns,
                                  int l_run, int max_ins,
                                  int ref_type, int *types, int n_types,
                                  int *score) {
    // FIXME: n_types has a maximum; no need to alloc - use a #define?
    int sc[MAX_TYPES], sumq[MAX_TYPES], s, i, j, t, K, n_alt, tmp;
    memset(sumq, 0, n_types * sizeof(int));
    for (s = K = 0; s < n; ++s) {
        for (i = 0; i < n_plp[s]; ++i, ++K) {
            bam_pileup1_t *p = plp[s] + i;
            // Labelling is confusing here.
            //    sct is short for score.
            //    sc is score + t(type)
            // Why aren't these variable names reversed?
            int *sct = &score[K*n_types], seqQ, indelQ;
            for (t = 0; t < n_types; ++t) sc[t] = sct[t]<<6 | t;
            for (t = 1; t < n_types; ++t) // insertion sort
                for (j = t; j > 0 && sc[j] < sc[j-1]; --j)
                    tmp = sc[j], sc[j] = sc[j-1], sc[j-1] = tmp;

            /* errmod_cal() assumes that if the call is wrong, the
             * likelihoods of other events are equal. This is about
             * right for substitutions, but is not desired for
             * indels. To reuse errmod_cal(), I have to make
             * compromise for multi-allelic indels.
             */
            if ((sc[0]&0x3f) == ref_type) {
                // sc >> 14 is the total score.  It's been shifted by 8
                // from normalised score and 6 from type.
                indelQ = (sc[1]>>14) - (sc[0]>>14);
                seqQ = est_seqQ(bca, types[sc[1]&0x3f], l_run);
            } else {
                for (t = 0; t < n_types; ++t) // look for the reference type
                    if ((sc[t]&0x3f) == ref_type) break;
                indelQ = (sc[t]>>14) - (sc[0]>>14);
                seqQ = est_seqQ(bca, types[sc[0]&0x3f], l_run);
            }

            tmp = sc[0]>>6 & 0xff; // normalised score

            // reduce indelQ
            // high score = bad, low score = good.
            // low normalised scores leave indelQ unmodified
            // high normalised scores set indelQ to 0
            // inbetween scores have a linear scale from indelQ to 0
            indelQ = tmp > 111? 0 : (int)((1. - tmp/111.) * indelQ + .499);

            // Doesn't really help accuracy, but permits -h to take
            // affect still.
            if (indelQ > seqQ) indelQ = seqQ;
            if (indelQ > 255) indelQ = 255;
            if (seqQ > 255) seqQ = 255;
            p->aux = (sc[0]&0x3f)<<16 | seqQ<<8 | indelQ; // use 22 bits in total
            sumq[sc[0]&0x3f] += indelQ < seqQ? indelQ : seqQ; // FIXME: redunctant; always indelQ
//            fprintf(stderr, "  read=%d:%d name=%s call=%d indelQ=%d seqQ=%d\n", s, i, bam_get_qname(p->b), types[sc[0]&0x3f], indelQ, seqQ);
        }
    }
    // determine bca->indel_types[] and bca->inscns
    bca->maxins = max_ins;
    bca->inscns = (char*) realloc(bca->inscns, bca->maxins * 4);
    if (bca->maxins && !bca->inscns)
        return -1;
    for (t = 0; t < n_types; ++t)
        sumq[t] = sumq[t]<<6 | t;
    for (t = 1; t < n_types; ++t) // insertion sort
        for (j = t; j > 0 && sumq[j] > sumq[j-1]; --j)
            tmp = sumq[j], sumq[j] = sumq[j-1], sumq[j-1] = tmp;
    for (t = 0; t < n_types; ++t) // look for the reference type
        if ((sumq[t]&0x3f) == ref_type) break;
    if (t) { // then move the reference type to the first
        tmp = sumq[t];
        for (; t > 0; --t) sumq[t] = sumq[t-1];
        sumq[0] = tmp;
    }
    for (t = 0; t < 4; ++t) bca->indel_types[t] = B2B_INDEL_NULL;
    for (t = 0; t < 4 && t < n_types; ++t) {
        bca->indel_types[t] = types[sumq[t]&0x3f];
        if (bca->maxins)
            memcpy(&bca->inscns[t * bca->maxins],
                   &inscns[(sumq[t]&0x3f) * max_ins], bca->maxins);
    }
    // update p->aux
    for (s = n_alt = 0; s < n; ++s) {
        for (i = 0; i < n_plp[s]; ++i) {
            bam_pileup1_t *p = plp[s] + i;
            int x = types[p->aux>>16&0x3f];
            for (j = 0; j < 4; ++j)
                if (x == bca->indel_types[j]) break;
            p->aux = j<<16 | (j == 4? 0 : (p->aux&0xffff));
            if ((p->aux>>16&0x3f) > 0) ++n_alt;
            //fprintf(stderr, "X pos=%d read=%d:%d name=%s call=%d type=%d seqQ=%d indelQ=%d\n", pos, s, i, bam_get_qname(p->b), (p->aux>>16)&0x3f, bca->indel_types[(p->aux>>16)&0x3f], (p->aux>>8)&0xff, p->aux&0xff);
        }
    }

    return n_alt;
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
int bcf_call_gap_prep(int n, int *n_plp, bam_pileup1_t **plp, int pos,
                      bcf_callaux_t *bca, const char *ref)
{
    if (ref == 0 || bca == 0) return -1;

    int i, s, t, n_types, *types, max_rd_len, left, right, max_ins;
    int *score, max_ref2;
    int N, K, l_run, ref_type, n_alt;
    char *inscns = 0, *ref1, *ref2, *query;

    // determine if there is a gap
    for (s = N = 0; s < n; ++s) {
        for (i = 0; i < n_plp[s]; ++i)
            if (plp[s][i].indel != 0) break;
        if (i < n_plp[s]) break;
    }
    if (s == n)
        // there is no indel at this position.
        return -1;

    // find out how many types of indels are present
    types = bcf_cgp_find_types(n, n_plp, plp, pos, bca, ref,
                               &max_rd_len, &n_types, &ref_type, &N);
    if (!types)
        return -1;


    // calculate left and right boundary
    left = pos > bca->indel_win_size ? pos - bca->indel_win_size : 0;
    right = pos + bca->indel_win_size;
    int del_size = types[0]<0 ? -types[0] : 0;
    right += del_size;

    // in case the alignments stand out the reference
    for (i = pos; i < right; ++i)
        if (ref[i] == 0) break;
    right = i;

    // FIXME: move to own function: STR_adj_left_right?
    if (0) {
        rep_ele *reps, *elt, *tmp;

        // Convert ASCII to 0,1,2,3 seq for find_STR usage
        int j;
        char ref4[1024]; // FIXME, check!
        if (right > left+1024)
            right = left+1024;
        for (j = 0, i = left; i < right; i++, j++) {
            switch(ref[i]) {
            case 'A': ref4[j] = 0; break;
            case 'C': ref4[j] = 1; break;
            case 'G': ref4[j] = 2; break;
            case 'T': ref4[j] = 3; break;
            default:  ref4[j] = j%4; break; // mix N across all 4 types
            }
        }
        reps = find_STR(ref4, right-left, 0);

        //fprintf(stderr, "\nRef at %d: %.*s\n", left, right-left, ref+left);

#if 0
        int adjusted = 1;
        int over_l = pos-1;
        int over_r = pos+del_size+1;
        int ins_size = types[0]>0 ? types[0] : 0;
        DL_FOREACH_SAFE(reps, elt, tmp) {
            //fprintf(stderr, "rep %d..%d: %.*s\n", elt->start, elt->end,
            //        elt->end-elt->start+1, ref+left+elt->start);
            if (elt->start + left < over_l && elt->end + left >= pos-1) {
                over_l = elt->start + left;
                //fprintf(stderr, "Adj left\n");
                adjusted=1;
            }
            if (elt->end + left > over_r && elt->start + left <= pos+1) {
                over_r = elt->end + left;
                //fprintf(stderr, "Adj right\n");
                adjusted=1;
            }
            //DL_DELETE(reps, elt);
            //free(elt);
        }

        // 2nd pass, adjusting to next STR so require 2 STRs out
        if (adjusted) {
            int pos_l = over_l;
            int pos_r = over_r;
            DL_FOREACH_SAFE(reps, elt, tmp) {
                if (elt->start + left < over_l && elt->end + left >= pos_l-1)
                    over_l = elt->start + left;
                if (elt->end + left > over_r && elt->start + left <= pos_r+1)
                    over_r = elt->end + left;
                DL_DELETE(reps, elt);
                free(elt);
            }
        }
        //fprintf(stderr, "STR overlap = %d..(%d)..%d\n", over_l, pos, over_r);

        // FIXME adjustable param
        over_l = pos - (pos-over_l)*2;
        over_r = pos + (over_r-pos)*2;
        //over_l -= 5+del_size+ins_size;
        //over_r += 5+del_size+ins_size;

        over_l -= 5+3*(del_size+ins_size);
        over_r += 5+3*(del_size+ins_size);
        //fprintf(stderr, "=>  overlap = %d..(%d)..%d\n", over_l, pos, over_r);
        if (left < over_l)
            left = over_l;
        if (right > over_r)
            right = over_r;
#else
        // Too many FNs, but OK otherwise.
        char str[1024] = {0};
        const int n = 3;
        DL_FOREACH_SAFE(reps, elt, tmp) {
            int i, i_start = MAX(elt->start-n, 0), i_end = MIN(elt->end+n, 1024);
//            fprintf(stderr, "rep %d..%d: %.*s\n", elt->start, elt->end,
//                    elt->end-elt->start+1, ref+left+elt->start);
            for (i = i_start; i < i_end; i++)
                str[i]=1;
            DL_DELETE(reps, elt);
            free(elt);
        }
        int score;
        for (score = 3, i = pos; i > left && score; i--)
            score -= str[i-left]==0;
        int left_new = i;

        for (score = 3, i = pos; i < right && score; i++)
            score -= str[i-left]==0;
        int right_new = i;

        fprintf(stderr, "left %d, %d, pos %d, %d, right %d\n",
                left, left_new, pos, right_new, right);

        left = left_new;
        right = right_new;
#endif
    }

//    fprintf(stderr, "=== POS %d, left/right = len %d\n", pos, right-left);

    // compute the likelihood given each type of indel for each read
    max_ins = types[n_types - 1];   // max_ins is at least 0
    max_ref2 = right - left + 2 + 2 * (max_ins > -types[0]? max_ins : -types[0]);
    // FIXME: add fudge to permit some extra neighbouring indels
    max_ref2 += 50;

    // The length of the homopolymer run around the current position
    l_run = bcf_cgp_l_run(ref, pos);

    // construct the consensus sequence (minus indels, which are added later)
    if (max_ins > 0) {
        inscns = bcf_cgp_calc_cons(n, n_plp, plp, pos,
                                   types, n_types, max_ins, s);
        if (!inscns)
            return -1;
    }

    ref1  = (char*) calloc(max_ref2, 1);
    ref2  = (char*) calloc(max_ref2, 1);
    query = (char*) calloc(right - left + max_rd_len + max_ins + 2, 1);
    score = (int*) calloc(N * n_types, sizeof(int));
    bca->indelreg = 0;
    double nqual_over_60 = bca->nqual / 60.0;

    // FIXME: need additional types, or rather to amend the type 0 case?
    //
    // We have types matching indel, plus type 0 which is ref.
    // What about type 0 which matches consensus?
    // Eg we have a small (wrong) 1bp insertion at current location,
    // and a larger (correct) homozygous insertion say 10 bp away.
    //
    // We don't want the alignment of seqs vs wrong indel-hypothesis to be
    // scoring higher than against ref.  So need a consensus with the large
    // insertion and no small hypothesised one.

    int biggest_del = 0;
    for (t = 0; t < n_types; t++)
        if (biggest_del > types[t])
            biggest_del = types[t];

    for (t = 0; t < n_types; ++t) {
        int l, ir;

        // compute indelreg
        if (types[t] == 0)
            ir = 0;
        else if (types[t] > 0)
            ir = est_indelreg(pos, ref, types[t], &inscns[t*max_ins]);
        else
            ir = est_indelreg(pos, ref, -types[t], 0);

        if (ir > bca->indelreg)
            bca->indelreg = ir;

        // Identify max deletion length
        int max_deletion = 0;
        for (s = 0; s < n; ++s) {
            for (i = 0; i < n_plp[s]; ++i, ++K) {
                bam_pileup1_t *p = plp[s] + i;
                if (max_deletion < -p->indel)
                    max_deletion = -p->indel;
            }
        }

        // Realignment score, computed via BAQ
        for (s = K = 0; s < n; ++s) {
            char **tcons, *cp;
            int left_shift, right_shift;
            tcons = bcf_cgp_consensus(n, n_plp, plp, pos, bca, ref,
                                      left, right, s, types[t], biggest_del, 
                                      &left_shift, &right_shift);
            fprintf(stderr, "Cons0 @ %d %4d/%3d %s\n", pos, types[t], left_shift, tcons[0]);
            fprintf(stderr, "Cons1 @ %d %4d/%3d %s\n", pos, types[t], left_shift, tcons[1]);

            // FIXME: map from ascii to 0,1,2,3,4.
            // This is only needed because bcf_cgp_consensus is reporting in ASCII
            // currently, for ease of debugging.
            int tcon_len[2], cnum;
            for (cnum = 0; cnum < 2; cnum++) {
                for (cp = tcons[cnum]; *cp; cp++) {
                    switch(*cp) {
                    case 'A': *cp = 0; break;
                    case 'C': *cp = 1; break;
                    case 'G': *cp = 2; break;
                    case 'T': *cp = 3; break;
                    default : *cp = 4; break;
                    }
                }
                tcon_len[cnum] = cp-tcons[cnum];
            }

            // original consensus method
            //memcpy(ref1, ref2, right-left+(types[t]>0?types[t]:0));
            memcpy(ref1, tcons[1], MIN(tcon_len[1], max_ref2));
            if (tcon_len[1] < right-left+(types[t]>0?types[t]:0)) {
                memset(ref1+tcon_len[1], 4,
                       right-left+(types[t]>0?types[t]:0) - tcon_len[1]);
            }
//            fprintf(stderr, "Type %d = %2d\t", t, types[t]);
//            for (j = 0; j < right-left+(types[t]>0?types[t]:0); j++)
//                putc("ACGTN"[(uint8_t)ref2[j]], stderr);
//            putc('\n', stderr);

            // Our computed consensus may start/end in slightly different
            // positions due to indels.
            // We pad it out with Ns so sequences overlapping don't
            // carry penalties.  (Ideally we'd pad with the reference, but
            // this suffices and it's tricky to track.)
            int ref2_pos = 0;
            int rright = left + tcon_len[0]; // ref left/right
            if (left_shift > 0) {
                memset(ref2, 4/*N*/, MIN(left_shift, max_ref2));
                ref2_pos += MIN(left_shift, max_ref2);
            }
            memcpy(ref2 + ref2_pos, tcons[0], MIN(tcon_len[0], max_ref2-ref2_pos));
            ref2_pos += MIN(tcon_len[0], max_ref2-ref2_pos);
            if (right_shift > 0) {
                memset(ref2 + ref2_pos, 4/*N*/,
                       MIN(right_shift, max_ref2-ref2_pos));
            }

//            fprintf(stderr, "TYPE %d = %2d\t", t, types[t]);
//            for (j = 0; j < rright-left && j < max_ref2; j++)
//                putc("ACGTN"[(uint8_t)ref2[j]], stderr);
//            putc('\n', stderr);

            // align each read to ref2
            for (i = 0; i < n_plp[s]; ++i, ++K) {
                bam_pileup1_t *p = plp[s] + i;

                // Some basic ref vs alt stats.
                int imq = p->b->core.qual > 59 ? 59 : p->b->core.qual;
                imq *= nqual_over_60;

                int sc_len, slen, epos, sc_end;

                // Only need to gather stats on one type, as it's
                // identical calculation for all the subsequent ones
                // and we're sharing the same stats array
                if (t == 0) {
                    // Gather stats for INFO field to aid filtering.
                    // mq and sc_len not very helpful for filtering, but could
                    // help in assigning a better QUAL value.
                    //
                    // Pos is slightly useful.
                    // Base qual can be useful, but need qual prior to BAQ?
                    // May need to cache orig quals in aux tag so we can fetch
                    // them even after mpileup step.
                    get_pos(bca, p, &sc_len, &slen, &epos, &sc_end);

                    assert(imq >= 0 && imq < bca->nqual);
                    assert(epos >= 0 && epos < bca->npos);
                    assert(sc_len >= 0 && sc_len < 100);
                    if (p->indel) {
                        bca->ialt_mq[imq]++;
                        bca->ialt_scl[sc_len]++;
                        bca->ialt_pos[epos]++;
                    } else {
                        bca->iref_mq[imq]++;
                        bca->iref_scl[sc_len]++;
                        bca->iref_pos[epos]++;
                    }
                }

                int qbeg, qpos, qend, tbeg, tend, kk;
                uint8_t *seq = bam_get_seq(p->b);
                uint32_t *cigar = bam_get_cigar(p->b);
                if (p->b->core.flag & BAM_FUNMAP) continue;

                // FIXME: the following loop should be better moved outside;
                // nonetheless, realignment should be much slower anyway.
                for (kk = 0; kk < p->b->core.n_cigar; ++kk)
                    if ((cigar[kk]&BAM_CIGAR_MASK) == BAM_CREF_SKIP)
                        break;
                if (kk < p->b->core.n_cigar)
                    continue;

                // determine the start and end of sequences for alignment
                // FIXME: loops over CIGAR multiple times
                int left2 = left, right2 = right;
                if (p->b->core.l_qseq > 1000) {
                    // long read data needs less context.  It also tends to
                    // have many more candidate indels to investigate so
                    // speed here matters more.
                    if (pos - left >= bca->indel_win_size)
                        left2 += bca->indel_win_size/2;
                    if (right-pos >= bca->indel_win_size)
                        right2 -= bca->indel_win_size/2;
                }

                int r_start = p->b->core.pos;
                int r_end = bam_cigar2rlen(p->b->core.n_cigar,
                                           bam_get_cigar(p->b))
                            -1 + r_start;

                qbeg = tpos2qpos(&p->b->core, bam_get_cigar(p->b), left2,
                                 0, &tbeg);
                qpos = tpos2qpos(&p->b->core, bam_get_cigar(p->b), pos,
                                     0, &tend) - qbeg;
                qend = tpos2qpos(&p->b->core, bam_get_cigar(p->b), right2,
                                 1, &tend);

                if (types[t] < 0) {
                    int l = -types[t];
                    tbeg = tbeg - l > left?  tbeg - l : left;
                }
                if (left_shift < 0)
                    tbeg = tbeg + left_shift > left ? tbeg + left_shift : left;

                // FIXME: Why +20?  tbeg-left_shift to tend+right_shift
                // is still insufficient.  Why?  Check tpos2qpos maybe?
                if (left_shift+20 > 0)
                    tbeg = tbeg - (left_shift+20) > left
                         ? tbeg - (left_shift+20)
                         : left;
                if (right_shift+20 > 0)
                    tend = tend + right_shift+20 < rright
                         ? tend + right_shift+20
                         : rright;

                // write the query sequence
                for (l = qbeg; l < qend; ++l)
                    query[l - qbeg] = seq_nt16_int[bam_seqi(seq, l)];

                // A fudge for now.  Consider checking SAM header for
                // RG platform field.
                int long_read = p->b->core.l_qseq > 1000;

                // FIXME: we can improve (see above).
                // Maybe use tbeg/tend as before, but with adjustment for
                // difference between right-left and tcon_len.
                // For now we just brute force it and do full ref range.
                // It doesn't seem to impact on band at all.  *Why?*
                int tend1 = left + tcon_len[0] - (left2-left);
                int tend2 = left + tcon_len[1] - (left2-left);


                // do realignment; this is the bottleneck.
                //
                // Note low score = good, high score = bad.
                if (tend > tbeg) {
                    if (bcf_cgp_align_score(p, bca, types[t],
                                            //(uint8_t *)ref1 + left2-left,
                                            //(uint8_t *)ref2 + left2-left,
                                            (uint8_t *)tcons[0] + left2-left,
                                            (uint8_t *)tcons[1] + left2-left,
                                            (uint8_t *)query,
                                            r_start, r_end, long_read,
                                            tbeg, tend1, tend2, left2, rright,
                                            qbeg, qend, qpos, max_deletion,
                                            &score[K*n_types + t]) < 0) {
                        score[K*n_types + t] = 0xffffff;
                        return -1;
                    }
                } else {
                    // place holder large cost for reads that cover the
                    // region entirely within a deletion (thus tend < tbeg).
                    score[K*n_types + t] = 0xffffff;
                }
#if 0
                for (l = 0; l < tend - tbeg + abs(types[t]); ++l) {
                    if (tbeg-left+l >= max_ref2)
                        break;
                    fputc("ACGTN"[(int)ref2[tbeg-left+l]], stderr);
                }
                fputc('\n', stderr);
                for (l = 0; l < qend - qbeg; ++l)
                    fputc("ACGTN"[(int)query[l]], stderr);
                fputc('\n', stderr);
                fprintf(stderr, "pos=%d type=%d read=%d:%d name=%s "
                        "qbeg=%d tbeg=%d score=%d,%d\n",
                        pos, types[t], s, i, bam_get_qname(p->b),
                        qbeg, tbeg, score[K*n_types + t]>>8, score[K*n_types + t]&0xff);
#endif
            }
            free(tcons);
        }
    }

    // compute indelQ
    n_alt = bcf_cgp_compute_indelQ(n, n_plp, plp, bca, inscns, l_run, max_ins,
                                   ref_type, types, n_types, score);

    // free
    free(ref1);
    free(ref2);
    free(query);
    free(score);
    free(types);
    free(inscns);

    return n_alt > 0? 0 : -1;
}
