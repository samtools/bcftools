/* The MIT License

   Copyright (c) 2013-2014 Genome Research Ltd.
   Authors:  see http://github.com/samtools/bcftools/blob/master/AUTHORS

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <htslib/hts.h>
#include "HMM.h"

static inline void multiply_matrix(int n, double *a, double *b, double *dst)
{
    double *out = dst;
    if ( a==dst || b==dst )
        out = (double*) malloc(sizeof(double)*n*n);

    int i,j,k;
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            double val = 0;
            for (k=0; k<n; k++) val += MAT(a,n,i,k)*MAT(b,n,k,j);
            MAT(out,n,i,j) = val;
        }
    }
    if ( out!=dst )
    {
        memcpy(dst,out,sizeof(double)*n*n);
        free(out);
    }
}

hmm_t *hmm_init(int nstates, double *tprob, int ntprob)
{
    hmm_t *hmm = (hmm_t*) calloc(1,sizeof(hmm_t));
    hmm->nstates = nstates;
    hmm->tprob = (double*) malloc(sizeof(double)*nstates*nstates);

    hmm_set_tprob(hmm, tprob, ntprob);

    return hmm;
}

void hmm_set_tprob(hmm_t *hmm, double *tprob, int ntprob)
{
    hmm->ntprob_arr = ntprob;
    if ( ntprob<=0 ) ntprob = 1;

    if ( !hmm->tprob_arr )
        hmm->tprob_arr  = (double*) malloc(sizeof(double)*hmm->nstates*hmm->nstates*ntprob);

    memcpy(hmm->tprob_arr,tprob,sizeof(double)*hmm->nstates*hmm->nstates);

    int i;
    for (i=1; i<ntprob; i++)
        multiply_matrix(hmm->nstates, hmm->tprob_arr, hmm->tprob_arr+(i-1)*hmm->nstates*hmm->nstates, hmm->tprob_arr+i*hmm->nstates*hmm->nstates);
}

void hmm_set_tprob_func(hmm_t *hmm, set_tprob_f set_tprob, void *data)
{
    hmm->set_tprob = set_tprob;
    hmm->set_tprob_data = data;
}

static void _set_tprob(hmm_t *hmm, int pos_diff)
{
    assert( pos_diff>=0 );

    int i, n;

    n = hmm->ntprob_arr ? pos_diff % hmm->ntprob_arr : 0;  // n-th precalculated matrix
    memcpy(hmm->tprob, hmm->tprob_arr+n*hmm->nstates*hmm->nstates, sizeof(*hmm->tprob)*hmm->nstates*hmm->nstates);

    if ( hmm->ntprob_arr > 0  )
    {
        n = pos_diff / hmm->ntprob_arr;  // number of full blocks to jump
        for (i=0; i<n; i++)
            multiply_matrix(hmm->nstates, hmm->tprob_arr+(hmm->ntprob_arr-1)*hmm->nstates*hmm->nstates, hmm->tprob, hmm->tprob);
    }
}

void hmm_run_viterbi(hmm_t *hmm, int n, double *eprobs, uint32_t *sites)
{
    // Init arrays when run for the first time
    if ( hmm->nsites < n )
    {
        hmm->nsites = n;
        hmm->vpath = (uint8_t*) realloc(hmm->vpath, sizeof(uint8_t)*hmm->nsites*hmm->nstates);
    }
    if ( !hmm->vprob )
    {
        hmm->vprob     = (double*) malloc(sizeof(double)*hmm->nstates);
        hmm->vprob_tmp = (double*) malloc(sizeof(double)*hmm->nstates);
    }


    // Init all states with equal likelihood
    int i,j, nstates = hmm->nstates;
    for (i=0; i<nstates; i++) hmm->vprob[i] = 1./nstates;


    // Run Viterbi
    uint32_t prev_pos = sites[0];
    for (i=0; i<n; i++)
    {
        uint8_t *vpath = &hmm->vpath[i*nstates];
        double *eprob  = &eprobs[i*nstates];

        int pos_diff = sites[i] == prev_pos ? 0 : sites[i] - prev_pos - 1;

        _set_tprob(hmm, pos_diff);
        if ( hmm->set_tprob ) hmm->set_tprob(hmm, prev_pos, sites[i], hmm->set_tprob_data);

        double vnorm = 0;
        for (j=0; j<nstates; j++)
        {
            double vmax = 0;
            int k, k_vmax = 0;
            for (k=0; k<nstates; k++)
            {
                double pval = hmm->vprob[k] * MAT(hmm->tprob,hmm->nstates,j,k);
                if ( vmax < pval ) { vmax = pval; k_vmax = k; }
            }
            vpath[j] = k_vmax;
            hmm->vprob_tmp[j] = vmax * eprob[j];
            vnorm += hmm->vprob_tmp[j];
        }
        for (j=0; j<nstates; j++) hmm->vprob_tmp[j] /= vnorm;
        double *tmp = hmm->vprob; hmm->vprob = hmm->vprob_tmp; hmm->vprob_tmp = tmp;

        #if 0
            fprintf(stderr,"%d: vprob=", sites[i]+1);
            for (j=0; j<nstates; j++)  fprintf(stderr," %f", hmm->vprob[j]);
            fprintf(stderr,"\teprob=");
            for (j=0; j<nstates; j++)  fprintf(stderr," %f", eprob[j]);
            fprintf(stderr,"\tvpath=");
            for (j=0; j<nstates; j++)  fprintf(stderr," %d", vpath[j]);
            fprintf(stderr,"\ttprob=");
            for (j=0; j<nstates; j++) 
            {
                int k;
                for (k=0; k<nstates; k++) fprintf(stderr," %f", MAT(hmm->tprob,hmm->nstates,j,k));
            }
            fprintf(stderr,"\t %d\n", pos_diff);
        #endif

        prev_pos = sites[i];
    }

    // Find the most likely state
    int iptr = 0;
    for (i=1; i<nstates; i++) 
        if ( hmm->vprob[iptr] < hmm->vprob[i] ) iptr = i;

    // Trace back the Viterbi path, we are reusing vpath for storing the states (vpath[i*nstates])
    for (i=n-1; i>=0; i--)
    {
        assert( iptr<nstates && hmm->vpath[i*nstates + iptr]<nstates );
        iptr = hmm->vpath[i*nstates + iptr];
        hmm->vpath[i*nstates] = iptr;     // reusing the array for different purpose here
    }
}

void hmm_destroy(hmm_t *hmm)
{
    free(hmm->vprob);
    free(hmm->vpath);
    free(hmm->tprob);
    free(hmm->tprob_arr);
    free(hmm);
}

