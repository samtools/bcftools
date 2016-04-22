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

#include "peakfit.h"
#include <stdio.h>
#include <gsl/gsl_version.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <assert.h>
#include <math.h>

#define NPARAMS 5

// gauss params: sqrt(scale), center, sigma
typedef struct _peak_t
{
    int fit_mask;
    double params[NPARAMS], ori_params[NPARAMS];    // current and input parameters
    struct { int scan; double min, max, best; } mc[NPARAMS];    // monte-carlo settings and best parameter
    void (*calc_f)  (int nvals, double *xvals, double *yvals, void *args);
    void (*calc_df) (int nvals, double *xvals, double *yvals, double *dfvals, int idf, void *args);
    void (*print_func) (struct _peak_t *pk, kstring_t *str);
    void (*convert_get) (struct _peak_t *pk, double *params);
    double (*convert_set) (struct _peak_t *pk, int iparam, double value);
}
peak_t;

struct _peakfit_t
{
    int npeaks, mpeaks, nparams, mparams;
    peak_t *peaks;
    double *params;
    int nvals, mvals;
    double *xvals, *yvals, *vals;
    kstring_t str;
    int verbose, nmc_iter;
};


/*
    Gaussian peak with the center bound in the interval <d,e>:
        yi = scale^2 * exp(-(xi-z)^2/sigma^2)

        dy/dscale  = 2*scale * EXP
        dy/dcenter = -scale^2 * sin(center) * (e-d) * (xi - z) * EXP / sigma^2
        dy/dsigma  = 2*scale^2 * (xi - z)^2 * EXP / sigma^3

    where
        z   = 0.5*(cos(center)+1)*(e-d) + d
        EXP = exp(-(xi-z)^2/sigma^2)
*/
void bounded_gaussian_calc_f(int nvals, double *xvals, double *yvals, void *args)
{
    peak_t *pk = (peak_t*) args;

    double scale2 = pk->params[0] * pk->params[0];
    double center = pk->params[1];
    double sigma  = pk->params[2];
    double d = pk->params[3];
    double e = pk->params[4];
    double z = 0.5*(cos(center)+1)*(e-d) + d;

    int i;
    for (i=0; i<nvals; i++)
    {
        double tmp = (xvals[i] - z)/sigma;
        yvals[i] += scale2 * exp(-tmp*tmp);
    }
}
void bounded_gaussian_calc_df(int nvals, double *xvals, double *yvals, double *dfvals, int idf, void *args)
{
    peak_t *pk = (peak_t*) args;

    double scale  = pk->params[0];
    double center = pk->params[1];
    double sigma  = pk->params[2];
    double d = pk->params[3];
    double e = pk->params[4];
    double z = 0.5*(cos(center)+1)*(e-d) + d;

    int i;
    for (i=0; i<nvals; i++)
    {
        double EXP = exp(-(xvals[i]-z)*(xvals[i]-z)/sigma/sigma);
        double zi  = xvals[i] - z;
        if ( idf==0 )       // dscale
            dfvals[i] += 2*scale*EXP;
        else if ( idf==1 )  // dcenter
            dfvals[i] -= scale*scale*sin(center)*(e-d)*zi*EXP/sigma/sigma;
        else if ( idf==2 )  // dsigma
            dfvals[i] += 2*scale*scale*zi*zi*EXP/sigma/sigma/sigma;
    }
}
void bounded_gaussian_sprint_func(peak_t *pk, kstring_t *str)
{
    double center = pk->params[1];
    double d = pk->params[3];
    double e = pk->params[4];
    double z = 0.5*(cos(center)+1)*(e-d) + d;
    ksprintf(str,"%f**2 * exp(-(x-%f)**2/%f**2)",fabs(pk->params[0]),z,fabs(pk->params[2]));
}
double bounded_gaussian_convert_set(peak_t *pk, int iparam, double value)
{
    if ( iparam!=1 ) return value;
    double d = pk->ori_params[3];
    double e = pk->ori_params[4];
    if ( value<d ) value = d;
    else if ( value>e ) value = e;
    return acos(2*(value-d)/(e-d) - 1);
}
void bounded_gaussian_convert_get(peak_t *pk, double *params)
{
    params[0] = fabs(pk->params[0]);
    params[2] = fabs(pk->params[2]);
    double center = pk->params[1];
    double d = pk->params[3];
    double e = pk->params[4];
    params[1] = 0.5*(cos(center)+1)*(e-d) + d;
}

void peakfit_add_bounded_gaussian(peakfit_t *pkf, double a, double b, double c, double d, double e, int fit_mask)
{
    pkf->npeaks++;
    hts_expand0(peak_t,pkf->npeaks,pkf->mpeaks,pkf->peaks);

    int i, nfit = 0;
    for (i=0; i<NPARAMS; i++) 
        if ( fit_mask & (1<<i) ) nfit++;

    assert( d<e );

    pkf->nparams += nfit;
    hts_expand0(double,pkf->nparams,pkf->mparams,pkf->params);

    peak_t *pk = &pkf->peaks[pkf->npeaks-1];
    memset(pk, 0, sizeof(peak_t));

    pk->calc_f      = bounded_gaussian_calc_f;
    pk->calc_df     = bounded_gaussian_calc_df;
    pk->print_func  = bounded_gaussian_sprint_func;
    pk->convert_set = bounded_gaussian_convert_set;
    pk->convert_get = bounded_gaussian_convert_get;
    pk->fit_mask    = fit_mask;
    pk->ori_params[0] = a;
    pk->ori_params[2] = c;
    pk->ori_params[3] = d;
    pk->ori_params[4] = e;
    pk->ori_params[1] = pk->convert_set(pk, 1, b);
}


/*
    Gaussian peak:
        yi = scale^2 * exp(-(x-center)^2/sigma^2)

        dy/dscale  = 2 * scale * EXP
        dy/dcenter = 2 * scale^2 * (x-center) * EXP / sigma^2
        dy/dsigma  = 2 * scale^2 * (x-center)^2 * EXP / sigma^3

    where 
        EXP = exp(-(x-center)^2/sigma^2)
*/
void gaussian_calc_f(int nvals, double *xvals, double *yvals, void *args)
{
    peak_t *pk = (peak_t*) args;

    double scale2 = pk->params[0] * pk->params[0];
    double center = pk->params[1];
    double sigma  = pk->params[2];

    int i;
    for (i=0; i<nvals; i++)
    {
        double tmp = (xvals[i] - center)/sigma;
        yvals[i] += scale2 * exp(-tmp*tmp);
    }
}
void gaussian_calc_df(int nvals, double *xvals, double *yvals, double *dfvals, int idf, void *args)
{
    peak_t *pk = (peak_t*) args;

    double scale  = pk->params[0];
    double center = pk->params[1];
    double sigma  = pk->params[2];

    int i;
    for (i=0; i<nvals; i++)
    {
        double zi  = xvals[i] - center;
        double EXP = exp(-zi*zi/(sigma*sigma));
        if ( idf==0 )       // dscale
            dfvals[i] += 2*scale*EXP;
        else if ( idf==1 )  // dcenter
            dfvals[i] += 2*scale*scale*zi*EXP/(sigma*sigma);
        else if ( idf==2 )  // dsigma
            dfvals[i] += 2*scale*scale*zi*zi*EXP/(sigma*sigma*sigma);
    }
}
void gaussian_sprint_func(struct _peak_t *pk, kstring_t *str)
{
    ksprintf(str,"%f**2 * exp(-(x-%f)**2/%f**2)",fabs(pk->params[0]),pk->params[1],fabs(pk->params[2]));
}
void gaussian_convert_get(peak_t *pk, double *params)
{
    params[0] = fabs(pk->params[0]);
    params[1] = fabs(pk->params[1]);
    params[2] = fabs(pk->params[2]);
}


void peakfit_add_gaussian(peakfit_t *pkf, double a, double b, double c, int fit_mask)
{
    pkf->npeaks++;
    hts_expand0(peak_t,pkf->npeaks,pkf->mpeaks,pkf->peaks);

    int i, nfit = 0;
    for (i=0; i<NPARAMS; i++) 
        if ( fit_mask & (1<<i) ) nfit++;

    pkf->nparams += nfit;
    hts_expand0(double,pkf->nparams,pkf->mparams,pkf->params);

    peak_t *pk = &pkf->peaks[pkf->npeaks-1];
    memset(pk, 0, sizeof(peak_t));

    pk->calc_f      = gaussian_calc_f;
    pk->calc_df     = gaussian_calc_df;
    pk->print_func  = gaussian_sprint_func;
    pk->convert_get = gaussian_convert_get;
    pk->fit_mask    = fit_mask;
    pk->ori_params[0] = a;
    pk->ori_params[1] = b;
    pk->ori_params[2] = c;
}


/*
    exp peak:
        yi = scale^2 * exp((x-center)/sigma^2)

        dy/dscale  = 2 * scale * EXP
        dy/dcenter = -scale^2  * EXP / sigma^2
        dy/dsigma  = -2 * scale^2 * (x-center) * EXP / sigma^3

    where 
        EXP = exp((x-center)/sigma^2)
*/
void exp_calc_f(int nvals, double *xvals, double *yvals, void *args)
{
    peak_t *pk = (peak_t*) args;

    double scale2 = pk->params[0] * pk->params[0];
    double center = pk->params[1];
    double sigma  = pk->params[2];

    int i;
    for (i=0; i<nvals; i++)
    {
        yvals[i] += scale2 * exp((xvals[i]-center)/sigma/sigma);
    }
}
void exp_calc_df(int nvals, double *xvals, double *yvals, double *dfvals, int idf, void *args)
{
    peak_t *pk = (peak_t*) args;

    double scale  = pk->params[0];
    double center = pk->params[1];
    double sigma  = pk->params[2];

    int i;
    for (i=0; i<nvals; i++)
    {
        double EXP = exp((xvals[i]-center)/sigma/sigma);
        if ( idf==0 )       // dscale
            dfvals[i] += 2*scale*EXP;
        else if ( idf==2 )  // dsigma
            dfvals[i] -= 2*scale*scale*(xvals[i]-center)*EXP/sigma/sigma/sigma;
    }
}
void exp_sprint_func(struct _peak_t *pk, kstring_t *str)
{
    ksprintf(str,"%f**2 * exp((x-%f)/%f**2)",fabs(pk->params[0]),pk->params[1],fabs(pk->params[2]));
}
void exp_convert_get(peak_t *pk, double *params)
{
    params[0] = fabs(pk->params[0]);
    params[1] = fabs(pk->params[1]);
    params[2] = fabs(pk->params[2]);
}

void peakfit_add_exp(peakfit_t *pkf, double a, double b, double c, int fit_mask)
{
    pkf->npeaks++;
    hts_expand0(peak_t,pkf->npeaks,pkf->mpeaks,pkf->peaks);

    int i, nfit = 0;
    for (i=0; i<NPARAMS; i++) 
        if ( fit_mask & (1<<i) ) nfit++;

    assert( !(fit_mask&2) );

    pkf->nparams += nfit;
    hts_expand0(double,pkf->nparams,pkf->mparams,pkf->params);

    peak_t *pk = &pkf->peaks[pkf->npeaks-1];
    memset(pk, 0, sizeof(peak_t));

    pk->calc_f      = exp_calc_f;
    pk->calc_df     = exp_calc_df;
    pk->print_func  = exp_sprint_func;
    pk->convert_get = exp_convert_get;
    pk->fit_mask    = fit_mask;
    pk->ori_params[0] = a;
    pk->ori_params[1] = b;
    pk->ori_params[2] = c;
}


void peakfit_set_params(peakfit_t *pkf, int ipk, double *params, int nparams)
{
    peak_t *pk = &pkf->peaks[ipk];
    int i;
    if ( pk->convert_set )
        for (i=0; i<nparams; i++) pk->params[i] = pk->convert_set(pk, i, params[i]);
    else
        for (i=0; i<nparams; i++) pk->params[i] = params[i];
}

void peakfit_get_params(peakfit_t *pkf, int ipk, double *params, int nparams)
{
    peak_t *pk = &pkf->peaks[ipk];
    if ( pk->convert_get ) pk->convert_get(pk, params);
    else
    {
        int i;
        for (i=0; i<nparams; i++) params[i] = pk->params[i];
    }
}

peakfit_t *peakfit_init(void)
{
    return (peakfit_t*)calloc(1,sizeof(peakfit_t));
}

void peakfit_reset(peakfit_t *pkf)
{
    pkf->npeaks = pkf->nparams = 0;
    memset(pkf->peaks,0,sizeof(peak_t)*pkf->mpeaks);
}

void peakfit_destroy(peakfit_t *pkf)
{
    free(pkf->str.s);
    free(pkf->vals);
    free(pkf->params);
    free(pkf->peaks);
    free(pkf);
}

int peakfit_calc_f(const gsl_vector *params, void *data, gsl_vector *yvals)
{
    peakfit_t *pkf = (peakfit_t *) data;

    int i,j;
    for (i=0; i<pkf->nvals; i++)
        pkf->vals[i] = 0;

    int iparam = 0;
    for (i=0; i<pkf->npeaks; i++)
    {
        peak_t *pk = &pkf->peaks[i];
        for (j=0; j<NPARAMS; j++)
        {
            if ( !(pk->fit_mask & (1<<j)) ) continue;
            pk->params[j] = gsl_vector_get(params,iparam);
            iparam++;
        }
        pk->calc_f(pkf->nvals, pkf->xvals, pkf->vals, pk);
    }

    for (i=0; i<pkf->nvals; i++)
        gsl_vector_set(yvals, i, (pkf->vals[i] - pkf->yvals[i])/0.01);

    return GSL_SUCCESS;
}
int peakfit_calc_df(const gsl_vector *params, void *data, gsl_matrix *jacobian)
{
    peakfit_t *pkf = (peakfit_t *) data;

    int i,j,k,iparam = 0;
    for (i=0; i<pkf->npeaks; i++)
    {
        peak_t *pk = &pkf->peaks[i];
        int iparam_prev = iparam;
        for (j=0; j<NPARAMS; j++)
        {
            if ( !(pk->fit_mask & (1<<j)) ) continue;
            pk->params[j] = gsl_vector_get(params,iparam);
            iparam++;
        }
        iparam = iparam_prev;
        for (j=0; j<NPARAMS; j++)
        {
            if ( !(pk->fit_mask & (1<<j)) ) continue;
            for (k=0; k<pkf->nvals; k++) pkf->vals[k] = 0;
            pk->calc_df(pkf->nvals, pkf->xvals, pkf->yvals, pkf->vals, j, pk);
            for (k=0; k<pkf->nvals; k++) gsl_matrix_set(jacobian, k, iparam, pkf->vals[k]);
            iparam++;
        }
    }
    return GSL_SUCCESS;
}
int peakfit_calc_fdf(const gsl_vector *params, void *data, gsl_vector *yvals, gsl_matrix *jacobian)
{
    peakfit_calc_f(params, data, yvals);
    peakfit_calc_df(params, data, jacobian);
    return GSL_SUCCESS;
}

double peakfit_evaluate(peakfit_t *pkf)
{
    int i;
    for (i=0; i<pkf->nvals; i++)
        pkf->vals[i] = 0;

    for (i=0; i<pkf->npeaks; i++)
        pkf->peaks[i].calc_f(pkf->nvals, pkf->xvals, pkf->vals, &pkf->peaks[i]);

    double sum = 0;
    for (i=0; i<pkf->nvals; i++)
        sum += fabs(pkf->vals[i] - pkf->yvals[i]);

    return sum;
}

const char *peakfit_sprint_func(peakfit_t *pkf)
{
    pkf->str.l = 0;
    int i;
    for (i=0; i<pkf->npeaks; i++)
    {
        if ( i>0 ) kputs(" + ", &pkf->str);
        pkf->peaks[i].print_func(&pkf->peaks[i], &pkf->str);
    }
    return (const char*)pkf->str.s;
}

void peakfit_verbose(peakfit_t *pkf, int level)
{
    pkf->verbose = level;
}

void peakfit_set_mc(peakfit_t *pkf, double xmin, double xmax, int iparam, int niter)
{
    peak_t *pk = &pkf->peaks[ pkf->npeaks-1 ];
    pk->mc[iparam].scan = 1;
    pk->mc[iparam].min  = xmin;
    pk->mc[iparam].max  = xmax;
    pkf->nmc_iter = niter;
}

double peakfit_run(peakfit_t *pkf, int nvals, double *xvals, double *yvals)
{
    srand(0);   // for reproducibility

    pkf->nvals = nvals;
    pkf->xvals = xvals;
    pkf->yvals = yvals;
    hts_expand0(double,pkf->nvals,pkf->mvals,pkf->vals);
    if ( !pkf->nparams ) return peakfit_evaluate(pkf);

    gsl_multifit_function_fdf mfunc;
    mfunc.f   = &peakfit_calc_f;
    mfunc.df  = &peakfit_calc_df;
    mfunc.fdf = &peakfit_calc_fdf;
    mfunc.n   = nvals;
    mfunc.p   = pkf->nparams;
    mfunc.params = pkf;

    const gsl_multifit_fdfsolver_type *solver_type;
    gsl_multifit_fdfsolver *solver;
    solver_type = gsl_multifit_fdfsolver_lmsder;
    solver = gsl_multifit_fdfsolver_alloc(solver_type, nvals, mfunc.p);
    gsl_vector *grad = gsl_vector_alloc(pkf->nparams);

    int imc_iter, i,j, iparam;
    double best_fit = HUGE_VAL;
    for (imc_iter=0; imc_iter<=pkf->nmc_iter; imc_iter++)   // possibly multiple monte-carlo iterations
    {
        // set GSL parameters
        iparam = 0;
        for (i=0; i<pkf->npeaks; i++)
        {
            peak_t *pk = &pkf->peaks[i];
            for (j=0; j<NPARAMS; j++)
            {
                pk->params[j] = pk->ori_params[j];
                if ( pk->mc[j].scan )
                {
                    pk->params[j] = rand()*(pk->mc[j].max - pk->mc[j].min)/RAND_MAX + pk->mc[j].min;
                    if ( pk->convert_set ) pk->params[j] = pk->convert_set(pk, j, pk->params[j]);
                }
                if ( !(pk->fit_mask & (1<<j)) ) continue;
                pkf->params[iparam] = pk->params[j];
                iparam++;
            }
        }

        gsl_vector_view vview = gsl_vector_view_array(pkf->params, mfunc.p);
        gsl_multifit_fdfsolver_set(solver, &mfunc, &vview.vector);

        // iterate until convergence (or lack of it)
        int ret, test1 = 0, test2 = 0, niter = 0, niter_max = 500;
        do
        {
            ret = gsl_multifit_fdfsolver_iterate(solver);
            if ( pkf->verbose >1 )
            {
                fprintf(stderr, "%d: ", niter);
                for (i=0; i<pkf->npeaks; i++)
                {
                    peak_t *pk = &pkf->peaks[i];
                    fprintf(stderr,"\t%f %f %f", pk->params[0],pk->params[1],pk->params[2]);
                }
                fprintf(stderr, "\t.. %s\n", gsl_strerror(ret));
            }
            if ( ret ) break;

#if GSL_MAJOR_VERSION >= 2
            int info;
            test1 = gsl_multifit_fdfsolver_test(solver, 1e-8,1e-8, 0.0, &info);
#else
            gsl_multifit_gradient(solver->J, solver->f, grad);
            test1 = gsl_multifit_test_gradient(grad, 1e-8);
            test2 = gsl_multifit_test_delta(solver->dx, solver->x, 1e-8, 1e-8);
#endif
        }
        while ((test1==GSL_CONTINUE || test2==GSL_CONTINUE) && ++niter<niter_max);
        if ( pkf->verbose >1 )
        {
            fprintf(stderr,"test1=%s\n", gsl_strerror(test1));
            fprintf(stderr,"test2=%s\n", gsl_strerror(test2));
        }

        // recover parameters
        iparam = 0;
        for (i=0; i<pkf->npeaks; i++)
        {
            peak_t *pk = &pkf->peaks[i];
            for (j=0; j<NPARAMS; j++)
            {
                if ( !(pk->fit_mask & (1<<j)) ) continue;
                pk->params[j] = gsl_vector_get(solver->x, iparam++);
            }
        }

        // evaluate fit, update best parameters
        double fit = peakfit_evaluate(pkf);
        if ( fit<best_fit )
        {
            for (i=0; i<pkf->npeaks; i++)
            {
                peak_t *pk = &pkf->peaks[i];
                for (j=0; j<NPARAMS; j++) pk->mc[j].best = pk->params[j];
            }
        }
        if ( fit<best_fit ) best_fit = fit;
    }
    gsl_multifit_fdfsolver_free(solver);
    gsl_vector_free(grad);

    for (i=0; i<pkf->npeaks; i++)
    {
        peak_t *pk = &pkf->peaks[i];
        for (j=0; j<NPARAMS; j++) pk->params[j] = pk->mc[j].best;
    }
    return best_fit;
}


