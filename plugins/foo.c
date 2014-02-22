#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>

bcf_hdr_t *hdr;
int nsnps, nindels;

/* 
    This short description is used to generate the output of `bcftools annotate -l`.
*/
const char *about(void)
{
    return "A minimal plugin which counts number of SNPs and Indels and prints the result to STDERR.\n";
}


/* 
    Called once at startup allows to initialize any local variables.
*/
void init(bcf_hdr_t *h)
{
    hdr = h;
    nsnps = nindels = 0;
}


/*
    Called for each VCF record after all standard annotation things are finished.
    Return 0 on success, 1 to suppress the line from printing, -1 on critical errors.
*/
int process(bcf1_t *rec)
{
    int type = bcf_get_variant_types(rec);
    if ( type & VCF_SNP ) nsnps++;
    if ( type & VCF_INDEL ) nindels++;
    return 1;
}


/*
    Clean up.
*/
void destroy(void)
{
    fprintf(stderr,"Number of SNPs:    %d\n", nsnps);
    fprintf(stderr,"Number of Indels:  %d\n", nindels);
}


