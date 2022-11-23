/*  mileup2.h -- mpileup v2.0 API

   Copyright (C) 2022 Genome Research Ltd.

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
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

#ifndef __MPILEUP20_H__
#define __MPILEUP20_H__

typedef struct _mpileup_t mpileup_t;

// Various options. Key marked as 'can fail' shoud be set with mpileup_set() and the
// return status checked as these can fail.
typedef enum
{
    // required
    FASTA_REF,          // char*, fasta reference
    BAM,                // char*, add another BAM/CRAM. Can fail.
    BAM_FNAME,          // char*, read list of BAMs/CRAMs from a file. Can fail.

    // optional
    REGIONS,            // char*, comma-separated regions
    REGIONS_FNAME,      // char*, file with a list of regions
    TARGETS,
    TARGETS_FNAME,
    SAMPLES,            // char*, comma-separated samples to include
    SAMPLES_FNAME,      // char*, file of samples to include
    READ_GROUPS_FNAME,  // char*, read groups to include or exclude, see `bcftools mpileup -G` documentation

    // bit flags
    SMART_OVERLAPS,     // int {0,1}, 1:disable read-pair overlap detection [0]
    SMPL_IGNORE_RG,     // int {0,1}, 1:ignore read groups, one bam = one sample [0]

    // numeric parameters
    MAX_DP_PER_FILE,    // int, maximum depth per file, regardless of the number of samples [250]
    ADJUST_MQ,          // int, --adjust-MQ value for downgrading reads with excessive mismatches, 0 to disable [0]
    MIN_MQ,             // int, skip alignments with mapQ smaller than MIN_MQ [0]
    MIN_BQ,             // int, skip bases with baseQ smaller than MIN_BQ [1]
    MAX_BQ,             // int, cap BQ for overly optimistics platforms [60]
    DELTA_BQ,           // int, --delta-BQ value for downgrading high BQ if neighbor_BQ is < BQ - 30  [30]
    MIN_REALN_FRAC,     // float, minimum fraction of reads with evidence for an indel to consider the site needing realignment [0.05]
    MIN_REALN_DP,       // int, minimum number of reads with evidence for an indel to consider the site needing realignment [2]
    MAX_REALN_DP,       // int, subsample to this many reads for realignment [250]
    MAX_REALN_LEN,      // int, realign reads this long max [500]
}
mpileup_opt_t;

#define mpileup_set_opt(mplp,type,key,value) { type tmp = value; mpileup_set(mplp, key, (void*)&tmp); }
#define mpileup_get_opt(mplp,type,key) (*(type*)mpileup_get(mplp, key))

mpileup_t *mpileup_init(void);
int mpileup_set(mpileup_t *mplp, mpileup_opt_t key, void *value);   // returns 0 on success
void *mpileup_get(mpileup_t *mplp, mpileup_opt_t key);
void mpileup_destroy(mpileup_t *mplp);

#endif
