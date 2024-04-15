/*  plugins/allele-length.c -- Calculate stats about the length of alleles

    Copyright (C) 2017-2018 GENOMICS plc.
    Copyright (C) 2023 Genome Research Ltd.

    Author: Nicola Asuni <nicola.asuni@genomicsplc.com>

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
#include <inttypes.h>

#define MAXLEN 512

static uint64_t numvar;
static uint64_t numxvar;
static uint64_t reflen[MAXLEN];
static uint64_t altlen[MAXLEN];
static uint64_t refaltlen[MAXLEN];
static uint64_t xrefaltlen[MAXLEN];

const char *about(void)
{
    return "Count the frequency of the length of REF, ALT and REF+ALT\n";
}

const char *usage(void)
{
    return
        "\n"
        "About: Count the frequency of the length of alleles.\n"
        "Usage: bcftools +allele-length [General Options] \n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Example:\n"
        "   bcftools +allele-length in.vcf\n"
        "\n";
}

// return 0 if the string contains characters other than standard ACGT base letters
int contain_non_base(const char *str)
{
    int c;
    while ((c = *str++))
    {
        if ((c != 'A') && (c != 'a') && (c != 'C') && (c != 'c') && (c != 'G') && (c != 'g') && (c != 'T') && (c != 't'))
        {
            return 1;
        }
    }
    return 0;
}

// Called once at startup, it initializes local variables.
// Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    numvar = 0;
    int i = 0;
    for(i = 0; i < MAXLEN; i++) {
        reflen[i] = 0;
        altlen[i] = 0;
        refaltlen[i] = 0;
        xrefaltlen[i] = 0;
    }
    return 1;
}

// Called for each VCF record. Return rec to output the line or NULL to suppress output.
bcf1_t *process(bcf1_t *rec)
{
    int rl = strlen(rec->d.allele[0]);
    int al = strlen(rec->d.allele[1]);
    int ral = rl + al;
    if ( rl >= MAXLEN ) rl = MAXLEN - 1;
    if ( al >= MAXLEN ) al = MAXLEN - 1;
    if ( ral >= MAXLEN ) ral = MAXLEN - 1;
    reflen[rl] += 1;
    altlen[al] += 1;
    refaltlen[ral] += 1;
    if ((contain_non_base(rec->d.allele[0])) || (contain_non_base(rec->d.allele[1])))
    {
        xrefaltlen[ral] += 1;
        numxvar++;
    }
    numvar++;
    return NULL;
}

// Print final output
void destroy(void)
{
    int i = 0;
    printf("LENGTH\tREF\tALT\tREF+ALT\tREF+ALT WITH NON-BASE NUCLEOTIDES\n");
    for(i = 0; i < MAXLEN; i++) {
        printf("%d\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n", i, reflen[i], altlen[i], refaltlen[i], xrefaltlen[i]);
    }
    printf("\t\t\t%"PRIu64"\t%"PRIu64"\n", numvar, numxvar);
}
