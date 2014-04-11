#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>

bcf_hdr_t *hdr;
int nsnps, nindels, nsites;

/* 
    This short description is used to generate the output of `bcftools annotate -l`.
*/
const char *about(void)
{
    return 
        "A minimal plugin which counts number of SNPs, Indels, and\n"
        "total number of sites.\n";
}


/* 
    Called once at startup, allows to initialize local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/
int init(const char *opts, bcf_hdr_t *h)
{
    hdr  = h;
    nsnps = nindels = nsites = 0;
    return 1;
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
    nsites++;
    return 1;
}


/*
    Clean up.
*/
void destroy(void)
{
    printf("Number of SNPs:    %d\n", nsnps);
    printf("Number of Indels:  %d\n", nindels);
    printf("Number of sites:   %d\n", nsites);
}


