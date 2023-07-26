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
#include "estimate_ci.h"

extern void error(const char *format, ...);

struct eci_t_
{
    double *dist_val;
    uint64_t *dist_cnt;
    uint32_t dist_n;
    uint32_t nsmpl;     // number of data points we can afford to sample from
    uint32_t *smpl;     // representative sample of values to sample from
};

// initialize the array with a representative sample of values
static void init_vals(eci_t *eci)
{
    int i, j;
    uint64_t nvals = 0, nflush = 0;
    for (i=0; i<eci->dist_n; i++) nvals += eci->dist_cnt[i];
    eci->nsmpl = nvals < 1e4 ? nvals : 1e4;
    eci->smpl  = (uint32_t*) malloc(eci->nsmpl*sizeof(*eci->smpl));
    double val;
    i = 0, j = 0;
    while (i<eci->nsmpl && j<eci->dist_n)
    {
        while (j<eci->dist_n)
        {
            val = eci->dist_val[j];
            nflush += eci->dist_cnt[j++];
            if ( nflush >= nvals/eci->nsmpl || j>=eci->dist_n ) break;
        }
        int nset = nflush*eci->nsmpl/nvals;
        nflush = 0;
        fprintf(stderr,"ECI\t%d\t%d\t%f\n",i,nset,val);
        while (nset-- && i<eci->nsmpl) eci->smpl[i++] = val;
    }
}

eci_t *eci_init(double *val, uint64_t *cnt, uint32_t n)
{
    eci_t *eci = (eci_t*) calloc(1,sizeof(eci_t));
    eci->dist_val = val;
    eci->dist_cnt = cnt;
    eci->dist_n   = n;
    init_vals(eci);
    return eci;
}

void eci_destroy(eci_t *eci)
{
    if ( !eci ) return;
    free(eci->smpl);
    free(eci);
}

