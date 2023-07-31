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
    Estimate confidence intervals given a distribution and the number of missing
    data points. Used by gtcheck.
 */

#ifndef ESTIMATE_CI_H__
#define ESTIMATE_CI_H__

#include <stdio.h>
#include <inttypes.h>

typedef struct eci_t_ eci_t;

/*
    eci_init() - initialize from the provided distribution
    @val:   array of values (x-values of a histogram)
    @cnt:   counts (y-values of a histogram)
    @n:     array size
 */
eci_t *eci_init(double *val, uint64_t *cnt, uint32_t n);
int eci_calc(eci_t *eci, uint32_t nmissing, double *ci_min, double *ci_max);
void eci_destroy(eci_t *eci);

#endif

