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

#ifndef __HMM_H__
#define __HMM_H__

#define MAT(matrix,ndim,i,j) (matrix)[(ndim)*(i)+(j)]       // P(i|j), that is, transition j->i

typedef struct _hmm_t hmm_t;

typedef void (*set_tprob_f) (hmm_t *hmm, uint32_t prev_pos, uint32_t pos, void *data);

struct _hmm_t
{
    int nstates;    // number of states

    double *vprob, *vprob_tmp;  // viterbi probs [nstates]
    uint8_t *vpath;             // viterbi path [nstates*nsites]
    double *bwd, *bwd_tmp;      // bwd probs [nstates]
    double *fwd;                // fwd probs [nsites]
    int nsites;

    int ntprob_arr;             // number of pre-calculated tprob matrices
    double *tprob, *tprob_arr;  // array of transition matrices, precalculated to ntprob_arr positions
    set_tprob_f set_tprob;      // optional user function to set / modify transition probabilities
    void *set_tprob_data;
};

hmm_t *hmm_init(int nstates, double *tprob, int ntprob);
void hmm_set_tprob(hmm_t *hmm, double *tprob, int ntprob);
void hmm_set_tprob_func(hmm_t *hmm, set_tprob_f set_tprob, void *data);
void hmm_run_viterbi(hmm_t *hmm, int nsites, double *eprob, uint32_t *sites);
void hmm_destroy(hmm_t *hmm);

#endif

