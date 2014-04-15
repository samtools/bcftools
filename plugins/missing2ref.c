#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <inttypes.h>

bcf_hdr_t *in_hdr, *out_hdr;
int32_t *gts = NULL, mgts = 0;
uint64_t nchanged = 0;

const char *about(void)
{
    return "Set missing genotypes (\"./.\") to ref allele (\"0/0\").\n";
}

int init(const char *opts, bcf_hdr_t *in, bcf_hdr_t *out)
{
    in_hdr  = in;
    out_hdr = out;
    return 0;
}

int process(bcf1_t *rec)
{
    int ngts = bcf_get_genotypes(in_hdr, rec, &gts, &mgts);

    int i, changed = 0;
    for (i=0; i<ngts; i++)
    {
        if ( gts[i]==bcf_gt_missing ) 
        {
            gts[i] = bcf_gt_unphased(0);
            changed++;
        }
    }
    nchanged += changed;
    if ( changed ) bcf_update_genotypes(out_hdr, rec, gts, ngts);
    return 0;
}

void destroy(void) 
{
    fprintf(stderr,"Filled %"PRId64" REF alleles\n", nchanged);
    free(gts);
}


