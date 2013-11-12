#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>

bcf_hdr_t *hdr;

void init(bcf_hdr_t *h)
{
    hdr = h;
    fprintf(stderr,"Plugin \"foo\" initialized\n");
}

int process(bcf1_t *rec)
{
    return 1;
}

