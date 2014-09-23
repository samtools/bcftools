/*  plugins/counts.c -- counts SNPs, Indels, and total number of sites.

    Copyright (C) 2013, 2014 Genome Research Ltd.

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>

int nsnps, nindels, nsites;

/*
    This short description is used to generate the output of `bcftools plugin -l`.
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
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    nsnps = nindels = nsites = 0;
    return 1;
}


/*
    Called for each VCF record. Return rec to output the line or NULL
    to suppress output.
*/
bcf1_t *process(bcf1_t *rec)
{
    int type = bcf_get_variant_types(rec);
    if ( type & VCF_SNP ) nsnps++;
    if ( type & VCF_INDEL ) nindels++;
    nsites++;
    return NULL;
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


