#include <math.h>
#include "call.h"

int ccall(call_t *call, bcf1_t *rec) { return 0; }
int qcall(call_t *call, bcf1_t *rec) 
{ 
    // QCall format: 
    //  chromosome, position, reference allele, depth, mapping quality, 0, ..
    error("TODO: qcall output\n");
    return 0; 
}

void mcall_init(call_t *call) 
{ 
    call->pl2p = (double*) malloc(sizeof(double)*256);
    int i;
    for (i=0; i<256; i++)
        call->pl2p[i] = pow(10., -i/10.);
    return; 
}
void ccall_init(call_t *call) { return; }
void qcall_init(call_t *call) { return; }

void mcall_destroy(call_t *call) 
{ 
    free(call->pl2p);
    free(call->PLs);
    return; 
}
void ccall_destroy(call_t *call) { return; }
void qcall_destroy(call_t *call) { return; }


// Inits P(D|G): convert from LOG space and normalize
void set_pdg(call_t *call, int *PLs, double *pdg, int n_smpl, int n_gt)
{
    int i, j;
    for (i=0; i<n_smpl; i++)
    {
        double sum = 0;
        for (j=0; j<n_gt; j++)
        {
            assert( PLs[j]!=bcf_int32_vector_end && PLs[j]!=bcf_int32_missing && PLs[j]<256 );

            pdg[j] = call->pl2p[ PLs[j] ];
            sum += pdg[j];
        }
        assert( sum!=0 );
        for (j=0; j<n_gt; j++) pdg[j] /= sum;

        PLs += n_gt;
        pdg += n_gt;
    }
}

int mcall(call_t *call, bcf1_t *rec)
{
    int nals = rec->n_allele;
    if ( nals>4 )
        error("FIXME: Not ready for more than 4 alleles at %s:%d (%d)\n", call->hdr_in->id[BCF_DT_CTG][rec->rid].key,rec->pos+1, nals);

    // Get the genotype likelihoods
    bcf_fmt_t *pl_fmt = bcf_get_fmt(call->hdr_in, rec, "PL");
    if ( !pl_fmt ) error("PL tag not present at %s:%d\n", call->hdr_in->id[BCF_DT_CTG][rec->rid].key,rec->pos+1);
    assert( pl_fmt->n==nals*(nals+1)/2 );
//    bcf_get_format_int(pl_fmt, call->hdr_out->n[BCF_DT_SAMPLE], call->PLs, &call->nPLs);
    
    // Get sum of qualities
//    int nqsum = bcf_
//    bcf_info_t *qs_info = bcf_get_info(call->hdr_in, rec, "QS");

    // log Lb/La vs 0/La vs La/La

    return 0;
}

