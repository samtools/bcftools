/* The MIT License

   Copyright (c) 2023 Genome Research Ltd.

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
/*
    For estimating confidence intervals we need to sample efficiently from a
    distribution (a histogram of values). This we do by keeping an array with
    a representative sample of the values and use this to sample from. This
    looses granularity but is a reasonable approximation when done well.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <htslib/hts_os.h>
#include <htslib/kroundup.h>
#include "estimate_ci.h"

extern void error(const char *format, ...);

struct eci_t_
{
    // the distribution of values we sample from
    uint32_t nsmpl;     // number of data points we can afford to sample from
    double *smpl;       // representative sample of values to sample from

    // the simulations
    uint32_t imin,imax;         // indices to exp for the desired size of CI
    size_t kexp;                // number of parallel simulation experiments to run
    double *exp, *tmp;          // the cumulative outcome of nCI experiments, and temporary array
    double *ci_min, *ci_max;    // confidence intervals: min and max
    size_t mci, nci;            // number of allocated and precalculated CIs

    // short lived, only for initialization by init_vals()
    double *dist_val;
    uint64_t *dist_cnt;
    uint32_t dist_n;
};

// initialize the array with a representative sample of values
static void init_vals(eci_t *eci)
{
    int i, j;
    uint64_t nvals = 0, nflush = 0;
    for (i=0; i<eci->dist_n; i++) nvals += eci->dist_cnt[i];
    eci->nsmpl = nvals < 1e4 ? nvals : 1e4;
    eci->smpl  = (double*) malloc(eci->nsmpl*sizeof(*eci->smpl));
    double val, nsmpl_per_nvals = (double)eci->nsmpl / nvals;
    i = 0, j = 0;
    while (i<eci->nsmpl && j<eci->dist_n)
    {
        while (j<eci->dist_n)
        {
            val = eci->dist_val[j];
            nflush += eci->dist_cnt[j++];
            if ( nflush >= nvals/eci->nsmpl || j>=eci->dist_n ) break;
        }
        int nset = nflush * nsmpl_per_nvals;
        nflush = 0;
        while (nset-- && i<eci->nsmpl) eci->smpl[i++] = val;
    }
    eci->nsmpl = i;
}

eci_t *eci_init(double *val, uint64_t *cnt, uint32_t n)
{
    eci_t *eci = (eci_t*) calloc(1,sizeof(eci_t));
    eci->dist_val = val;
    eci->dist_cnt = cnt;
    eci->dist_n   = n;

    // Number of experiments (100) and CI percentiles (5-95)
    // todo: make these customizable
    eci->kexp = 100;
    eci->imin = 4;
    eci->imax = 94;
    eci->exp = (double*) calloc(eci->kexp,sizeof(*eci->exp));
    eci->tmp = (double*) malloc(eci->kexp*sizeof(*eci->tmp));

    init_vals(eci);
    return eci;
}

void eci_destroy(eci_t *eci)
{
    if ( !eci ) return;
    free(eci->ci_min);
    free(eci->ci_max);
    free(eci->exp);
    free(eci->tmp);
    free(eci->smpl);
    free(eci);
}

static int cmp_double(const void *a, const void *b)
{
    if ( *((double*)a) < *((double*)b) ) return -1;
    if ( *((double*)a) > *((double*)b) ) return 1;
    return 0;
}

/*
    eci_calc() - estimate confidence intervals given N missing values and
                 desired confidence level

    This is done by running K parallel "experiments" N times, i.e. we sample
    K values from the eci.smpl distribution N times. We save the intermediate
    results so that repeated calls are faster.

    todo: if nmissing is too big and either too much memory or too much time
    is required, add another layer of approximation and use mean and variance
    to estimate the estimators
 */
int eci_calc(eci_t *eci, uint32_t nmissing, double *ci_min, double *ci_max)
{
    int i,j;

    if ( !nmissing )
    {
        *ci_min = *ci_max = 0;
        return 0;
    }

    // expand arrays when necessary
    if ( eci->mci < nmissing )
    {
        size_t new_size = nmissing;
        kroundup_size_t(new_size);
        eci->ci_min = realloc(eci->ci_min,new_size*sizeof(*eci->ci_min));
        if ( !eci->ci_min ) error("Failed to allocate %zu bytes\n",new_size*sizeof(*eci->ci_min));
        eci->ci_max = realloc(eci->ci_max,new_size*sizeof(*eci->ci_min));
        if ( !eci->ci_max ) error("Failed to allocate %zu bytes\n",new_size*sizeof(*eci->ci_max));
        eci->mci = new_size;
    }

    // run the sampling (nmissing-nci)x times, 0 times if already have it
    for (i=eci->nci; i<nmissing; i++)
    {
        for (j=0; j<eci->kexp; j++)     // increment all K sampling experiments
        {
            eci->exp[j] += eci->smpl[(int)(eci->nsmpl*hts_drand48())];
            eci->tmp[j] = eci->exp[j];
        }
        // sort to determine the confidence interval
        qsort(eci->tmp, eci->kexp, sizeof(*eci->tmp), cmp_double);
        eci->ci_min[i] = eci->tmp[eci->imin];
        eci->ci_max[i] = eci->tmp[eci->imax];
    }
    if ( eci->nci < nmissing ) eci->nci = nmissing;
    *ci_min = eci->ci_min[nmissing-1];
    *ci_max = eci->ci_max[nmissing-1];
    return 0;
}

