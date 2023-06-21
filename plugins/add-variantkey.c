/*  plugins/add-variantkey.c -- add VariantKey INFO field.

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

#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <inttypes.h>
#include "../variantkey.h"

bcf_hdr_t *in_hdr, *out_hdr;

const char *about(void)
{
    return "Add VariantKey INFO fields VKX and RSX.\n";
}

const char *usage(void)
{
    return
        "\n"
        "About: Add VKX and RSX columns.\n"
        "Usage: bcftools +add-variantkey [General Options] \n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Example:\n"
        "   bcftools +add-variantkey in.vcf\n"
        "\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    in_hdr  = in;
    out_hdr = out;
    bcf_hdr_append(out_hdr, "##INFO=<ID=VKX,Number=1,Type=String,Description=\"Hexadecimal representation of 64 bit VariantKey\">");
    bcf_hdr_append(out_hdr, "##INFO=<ID=RSX,Number=1,Type=String,Description=\"Hexadecimal representation of ID minus the 'rs' prefix (32bit)\">");
    return 0;
}

bcf1_t *process(bcf1_t *rec)
{
    uint64_t vk = variantkey(
                      in_hdr->id[BCF_DT_CTG][rec->rid].key,
                      strlen(in_hdr->id[BCF_DT_CTG][rec->rid].key),
                      rec->pos,
                      rec->d.allele[0],
                      strlen(rec->d.allele[0]),
                      rec->d.allele[1],
                      strlen(rec->d.allele[1]));
    char vs[17];
    variantkey_hex(vk, vs);
    bcf_update_info_string(out_hdr, rec, "VKX", vs);
    char rsid[9];
    char *ptr = rec->d.id;
    ptr += 2; // remove 'rs'
    sprintf(rsid, "%08" PRIx32, (uint32_t)strtoul(ptr, NULL, 10));
    bcf_update_info_string(out_hdr, rec, "RSX", rsid);
    return rec;
}

void destroy(void)
{
}
