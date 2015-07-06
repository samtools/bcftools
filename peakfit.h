/* The MIT License

   Copyright (c) 2013-2015 Genome Research Ltd.

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

#ifndef PEAKFIT_H
#define PEAKFIT_H

typedef struct _peakfit_t peakfit_t;

peakfit_t *peakfit_init(void);
void peakfit_reset(peakfit_t *pkf);
void peakfit_destroy(peakfit_t *pkf);

void peakfit_add_gaussian(peakfit_t *pkf, double a, double b, double c, int fit_mask);
void peakfit_add_bounded_gaussian(peakfit_t *pkf, double a, double b, double c, double d, double e, int fit_mask);
void peakfit_add_exp(peakfit_t *pkf, double a, double b, double c, int fit_mask);
void peakfit_set_mc(peakfit_t *pkf, double xmin, double xmax, int iparam, int niter);
double peakfit_run(peakfit_t *pkf, int nvals, double *xvals, double *yvals);
double peakfit_evaluate(peakfit_t *pkf);

void peakfit_verbose(peakfit_t *pkf, int level);
void peakfit_get_params(peakfit_t *pkf, int ipk, double *params, int nparams);
void peakfit_set_params(peakfit_t *pkf, int ipk, double *params, int nparams);
const char *peakfit_sprint_func(peakfit_t *pkf);

#endif

