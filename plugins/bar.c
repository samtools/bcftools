#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>

bcf_hdr_t *hdr;

void init(bcf_hdr_t *h)
{
    hdr = h;
    fprintf(stderr,"Plugin \"bar\" initialized\n");
}

int process(bcf1_t *rec)
{
    fprintf(stderr,"bar: %d\n", rec->pos+1);
    return 1;
}

void destroy(void)
{
}


