/*  plugins/variantkey-hex.c -- Generate unsorted VariantKey lookup tables files in hexadecimal format.

    Copyright (C) 2017-2018 GENOMICS plc.

    Author: Nicola Asuni <nicola.asuni@tecnick.com>

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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <inttypes.h>
#include "../variantkey.h"

const char *FILE_VKRS = "vkrs.unsorted.hex";
const char *FILE_RSVK = "rsvk.unsorted.hex";
const char *FILE_NRVK = "nrvk.unsorted.tsv";

FILE *fp_vkrs; // VariantKey -> rsID
FILE *fp_rsvk; // rsID -> VariantKey
FILE *fp_nrvk; // VariantKey non-reversible map (maps VariantKey to REF and ALT)

static uint64_t numvar; // number of variants
static uint64_t nrv; // number of non-reversible variants

bcf_hdr_t *in_hdr;

const char *about(void)
{
    return "Generate VariantKey index files\n";
}

const char *usage(void)
{
    return
        "\n"
        "About: Generate unsorted VariantKey lookup tables files in hexadecimal format.\n"
        "Usage: bcftools +variantkey-hex [General Options] \n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Example:\n"
        "   bcftools +variantkey-hex in.vcf\n"
        "\n";
}

// Called once at startup, it initializes local variables.
// Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    in_hdr  = in;
    numvar = 0;
    char path[1024];
    char dir[1024] = "./";
    if (argc > 1)
    {
        strcpy(dir, argv[1]);
    }
    strcpy(path, dir);
    strcat(path, FILE_VKRS);
    fp_vkrs = fopen(path, "w");
    if (!fp_vkrs)
    {
        fprintf(stderr, "%s: %s\n", path, strerror(errno));
    }
    strcpy(path, dir);
    strcat(path, FILE_RSVK);
    fp_rsvk = fopen(path, "w");
    if (!fp_rsvk)
    {
        fprintf(stderr, "%s: %s\n", path, strerror(errno));
    }
    strcpy(path, dir);
    strcat(path, FILE_NRVK);
    fp_nrvk = fopen(path, "w");
    if (!fp_nrvk)
    {
        fprintf(stderr, "%s: %s\n", path, strerror(errno));
    }
    return 1;
}

// Called for each VCF record. Return rec to output the line or NULL to suppress output.
bcf1_t *process(bcf1_t *rec)
{
    int len_ref = strlen(rec->d.allele[0]);
    int len_alt = strlen(rec->d.allele[1]);
    uint64_t vk = variantkey(
                      in_hdr->id[BCF_DT_CTG][rec->rid].key,
                      strlen(in_hdr->id[BCF_DT_CTG][rec->rid].key),
                      rec->pos,
                      rec->d.allele[0],
                      len_ref,
                      rec->d.allele[1],
                      len_alt);
    char *ptr = rec->d.id;
    ptr += 2; // remove 'rs'
    uint32_t rs = (uint32_t)strtoul(ptr, NULL, 10);
    fprintf(fp_vkrs, "%016" PRIx64 "\t%08" PRIx32 "\n", vk, rs); // map VariantKey to rsID
    fprintf(fp_rsvk, "%08" PRIx32 "\t%016" PRIx64 "\n", rs, vk); // map rsID to VariantKey
    if (vk & 1)
    {
        // map VariantKey to REF and ALT
        fprintf(fp_nrvk, "%016" PRIx64 "\t%s\t%s\n", vk, rec->d.allele[0], rec->d.allele[1]);
        nrv++;
    }
    numvar++;
    return NULL;
}

void destroy(void)
{
    fclose(fp_vkrs);
    fclose(fp_rsvk);
    printf("VariantKeys: %" PRIu64 "\n", numvar);
    printf("Non-reversible VariantKeys: %" PRIu64 "\n", nrv);
}
