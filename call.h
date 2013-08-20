#ifndef __CALL_H__
#define __CALL_H__

#include <htslib/vcf.h>

#define CALL_KEEPALT 1
#define CALL_VARONLY (1<<1)

typedef struct
{
    // mcall only
    double min_ma_lrt;  // variant accepted if P(chi^2)>=FLOAT [0.99]
    int *PLs, nPLs, nqsum;
    double *pl2p, *qsum;

    // ccall only
    double indel_frac, theta, min_lrt, min_perm_p; 
    double prior_type, pref;
    int ngrp1_samples, n_perm;
    char *prior_file;

    // shared
    bcf1_t *rec;
    bcf_hdr_t *hdr_in, *hdr_out;

    uint32_t flag;
    uint8_t *ploidy;
    int nsamples;
    uint32_t trio;
}
call_t;

void error(const char *format, ...);

int mcall(call_t *call, bcf1_t *rec);    // multiallic and rare-variant calling model
int ccall(call_t *call, bcf1_t *rec);    // the default consensus calling model
int qcall(call_t *call, bcf1_t *rec);    // QCall output

void mcall_init(call_t *call);
void ccall_init(call_t *call);
void qcall_init(call_t *call);

void mcall_destroy(call_t *call);
void ccall_destroy(call_t *call);
void qcall_destroy(call_t *call);


#endif
