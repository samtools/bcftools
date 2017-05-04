/* plugins/vcf2bed.c -- convert a VCF file to the BED format. 

   Copyright (c) 2016 

   Author: Gaik Tamazian <mail (at) gtamazian (dot) com>

   The MIT License

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
   DEALINGS IN THE SOFTWARE.
*/
 
#include <stdio.h>
#include <htslib/vcf.h>

static const bcf_hdr_t *hdr_reader = NULL; 

const char *about(void)
{
    return "Convert a VCF file to the BED format.\n";
}

const char *usage(void)
{
    return
        "About: Convert a VCF file to the BED format and print\n" 
        "Usage: bcftools +vcf2bed [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   the plugin provides no options\n"
        "\n"
        "Example:\n"
        "   bcftools +vcf2bed in.vcf > out.bed\n"
        "\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    hdr_reader = in;
    return 1;
}

/* For each variant, produce a line in the BED format containing
 * a sequence ID (CHR), the variant starting and ending positions
 * (POS and POS + length of REF) and the variant name (ID). */
bcf1_t *process(bcf1_t *rec)
{
    bcf_unpack(rec, BCF_UN_STR);
    printf("%s\t%d\t%d\t%s\n", 
        bcf_hdr_id2name(hdr_reader, rec->rid),
        rec->pos, 
        rec->pos + rec->rlen,
        rec->d.id); 
    return NULL;
}

void destroy(void)
{
}

