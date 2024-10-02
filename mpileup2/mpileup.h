/*  mileup2.h -- mpileup v2.0 API

   Copyright (C) 2022-2023 Genome Research Ltd.

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

#include <htslib/sam.h>

typedef struct mpileup_t_ mpileup_t;

typedef enum
{
    // required
    FASTA_REF,          // char*; fasta reference
    BAM,                // char*; add another BAM/CRAM
    BAM_FNAME,          // char*; read list of BAMs/CRAMs from a file

    // optional
    REGIONS,            // char*; comma-separated regions
    REGIONS_FNAME,      // char*; file with a list of regions
    TARGETS,
    TARGETS_FNAME,
    SAMPLES,            // char*; comma-separated samples to include
    SAMPLES_FNAME,      // char*; file of samples to include
    READ_GROUPS_FNAME,  // char*; read groups to include or exclude, see `bcftools mpileup -G` documentation

    // bit flags
    SMART_OVERLAPS,     // int {0,1}; 1:disable read-pair overlap detection [0]
    SMPL_IGNORE_RG,     // int {0,1}; 1:ignore read groups, one bam = one sample [0]
    SKIP_ANY_UNSET,     // skip a read if any of the bits is not set [0]
    SKIP_ALL_UNSET,     // skip a read if all these bits are not set [0]
    SKIP_ANY_SET,       // skip a read if any of the bits is set [BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP]
    SKIP_ALL_SET,       // skip a read if all these bits are set [0]

    // numeric parameters
    MAX_DP_PER_SAMPLE,  // int; maximum depth per sample [250]
    ADJUST_MQ,          // int; --adjust-MQ value for downgrading reads with excessive mismatches, 0 to disable [0]
    MIN_MQ,             // int; skip alignments with mapQ smaller than MIN_MQ [0]
    MIN_BQ,             // int; skip bases with baseQ smaller than MIN_BQ [1]
    MAX_BQ,             // int; cap BQ for overly optimistics platforms [60]
    DELTA_BQ,           // int; --delta-BQ value for downgrading high BQ if neighbor_BQ is < BQ - 30  [30]
    MIN_REALN_FRAC,     // float; minimum fraction of reads with evidence for an indel to consider the site needing realignment [0.05]
    MIN_REALN_DP,       // int; minimum number of reads with evidence for an indel to consider the site needing realignment [2]
    MAX_REALN_DP,       // int; subsample to this many reads for realignment [250]
    MAX_REALN_LEN,      // int; realign reads this long max [500]

    // pileup related values
    N_SAMPLES,          // int; number of samples in the pileup
    N_READS,            // int*; number of reads per sample at the current position; can be called only once
    PILEUP,             // (todo)
    CHROM,              // char*; current chromosome
    POS,                // hts_pos_t; current positition in the pileup
    REF,                // char*,int*; get reference at the current position

    LEGACY_MODE,        // int; 0:new mode (unfinished, todo); 1:original mpileup code and API [0]
    LEGACY_PILEUP,      // bam_pileup1_t**; reads grouped by sample; fills N_READS array
}
mpileup_opt_t;

/**
 *  mpileup_alloc() - allocate the mpileup structure
 *
 *  Note that any default values set by the function can change without a warning and
 *  therefore should be set explicitly with `mpileup_set()` by the caller.
 */
mpileup_t *mpileup_alloc(void);

/**
 *  mpileup_destroy - destroy the mpileup structure
 */
void mpileup_destroy(mpileup_t *mplp);

/**
 *  mpileup_set() - set various options, see the mpileup_opt_t keys for the complete list
 *
 *  Returns 0 if the call succeeded, or negative number on error.
 */
int mpileup_set(mpileup_t *mplp, mpileup_opt_t key, ...);   // returns 0 on success

/**
 *  mpileup_get()     - get various options, see the mpileup_opt_t keys
 *  mpileup_get_val() - wrapper for `mpileup_get()` to return typed value
 *
 *  The former returns pointer to the memory area populated by the requested setting,
 *  its type can be inferred from the mpileup_opt_t documentation.
 */
void *mpileup_get(mpileup_t *mplp, mpileup_opt_t key, ...);
#define mpileup_get_val(mplp,type,key) (*(type*)mpileup_get(mplp, key))

/**
 *  mpileup_init() - inits regions, bams, iterators, etc.
 *
 *  Should be called after fully set up, i.e. after all `mpileup_set()` calls.
 *
 *  Returns 0 on success or negative number on error.
 */
int mpileup_init(mpileup_t *mplp);

/**
 *  mpileup_next() - positions mpileup to the next genomic position
 *
 *  Returns 1 when next position is available, 0 when all regions are done, negative value on error
 */
int mpileup_next(mpileup_t *mplp);

#endif
