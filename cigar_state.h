/*  cigar_state.h -- API for efficient parsing of CIGAR strings

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

#ifndef CIGAR_STATE_H
#define CIGAR_STATE_H

#include <stdint.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

typedef struct
{
    bam1_t *bam;
    uint32_t *cigar;
    uint8_t *seq;
    int ncig;
    int icig;           // position in the cigar string
    int iseq;           // the cigar[icigar] operation refers to seq[iseq+1]
    hts_pos_t ref_pos;  // reference coordinate, corresponds to iseq
}
cigar_state_t;

inline void cstate_init(cigar_state_t *cs, bam1_t *bam)
{
    cs->bam   = bam;
    cs->cigar = bam_get_cigar(bam);
    cs->seq   = bam_get_seq(bam);
    cs->ncig  = bam->core.n_cigar;
    cs->icig  = 0;
    cs->iseq  = 0;
    cs->ref_pos = bam->core.pos;
    int op  = cs->cigar[0] &  BAM_CIGAR_MASK;
    int len = cs->cigar[0] >> BAM_CIGAR_SHIFT;
    if ( op==BAM_CINS || op==BAM_CSOFT_CLIP ) cs->ref_pos -= len;
}

// Move in the cigar forward to find query index that matches the
// seek operator and the reference position.
//
// Returns the index to the query sequence cs->seq
// on success; -1 when there is no such matching position but the cigar
// is still not entirely consumed (e.g. a deletion or a soft-clip); -2
// when there is no overlap (i.e. the read ends before the position).
inline int cstate_seek_fwd(cigar_state_t *cs, hts_pos_t pos, int seek_op, int *oplen)
{
    while ( cs->ref_pos <= pos )
    {
        if ( cs->icig >= cs->ncig ) return -2;

        int op  = cs->cigar[cs->icig] &  BAM_CIGAR_MASK;
        int len = cs->cigar[cs->icig] >> BAM_CIGAR_SHIFT;
        if ( op==BAM_CMATCH || op==BAM_CEQUAL || op==BAM_CDIFF )
        {
            if ( cs->ref_pos + len <= pos )
            {
                cs->ref_pos += len;
                cs->iseq += len;
                cs->icig++;
                continue;
            }
            if ( seek_op==BAM_CMATCH ) return pos - cs->ref_pos + cs->iseq;
            return -1;
        }
        if ( op==BAM_CINS || op==BAM_CSOFT_CLIP )
        {
            if ( cs->ref_pos == pos + 1 && seek_op==op )
            {
                if ( oplen ) *oplen = len;
                return cs->iseq;
            }
            if ( cs->ref_pos >= pos ) return -1;
            cs->iseq += len;
            cs->icig++;
            continue;
        }
        if ( op==BAM_CDEL || op==BAM_CREF_SKIP )
        {
            if ( cs->ref_pos == pos && seek_op==op )
            {
                if ( oplen ) *oplen = len;
                return cs->iseq;
            }
            if ( cs->ref_pos >= pos ) return -1;
            cs->ref_pos += len;
            cs->icig++;
            continue;
        }
    }
    return cs->icig < cs->ncig ? -1 : -2;
}

#endif
