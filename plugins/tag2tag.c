/*  plugins/fill-AN-AC.c -- fills AN and AC INFO fields.

    Copyright (C) 2014 Genome Research Ltd.

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
#include <getopt.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include "bcftools.h"


#define GP_TO_GL 1

static int mode = 0, drop_source_tag = 0;
static bcf_hdr_t *in_hdr, *out_hdr;
static float *farr = NULL;
static int mfarr = 0;

const char *about(void)
{
    return "Convert between similar tags, such as GL and GP.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Convert between similar tags, such as GL and GP.\n"
        "Usage: bcftools +tag2tag [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
//todo        "       --gl-to-gp    convert FORMAT/GL to FORMAT/GP\n" 
        "       --gp-to-gl    convert FORMAT/GP to FORMAT/GL\n"
        "   -r, --replace     drop the source tag\n"
        "\n"
        "Example:\n"
        "   bcftools +tag2tag in.vcf -- -r --gp-to-gl\n"
        "\n";
}


static void init_header(bcf_hdr_t *hdr, const char *ori, int ori_type, const char *new_hdr_line)
{
    if ( ori )
        bcf_hdr_remove(hdr,ori_type,ori);

    bcf_hdr_append(hdr, new_hdr_line);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    static struct option loptions[] =
    {
        {"replace",0,0,'r'},
        {"gp-to-gl",0,0,1},
        {0,0,0,0}
    };
    char c;
    while ((c = getopt_long(argc, argv, "?hr",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case  1 : mode = GP_TO_GL; break;
            case 'r': drop_source_tag = 1; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( !mode ) mode = GP_TO_GL;

    in_hdr  = in;
    out_hdr = out;

    if ( mode==GP_TO_GL )
        init_header(out_hdr,drop_source_tag?"GP":NULL,BCF_HL_FMT,"##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihoods\">");

    return 0;
}

bcf1_t *process(bcf1_t *rec)
{
    int i, n;
    if ( mode==GP_TO_GL )
    {
        n = bcf_get_format_float(in_hdr,rec,"GP",&farr,&mfarr);
        for (i=0; i<n; i++)
        {
            if ( bcf_float_is_missing(farr[i]) || bcf_float_is_vector_end(farr[i]) ) continue;
            farr[i] = farr[i] ? log(farr[i]) : -99;
        }
        bcf_update_format_float(out_hdr,rec,"GL",farr,n);
        if ( drop_source_tag )
            bcf_update_format_float(out_hdr,rec,"GP",NULL,0);
    }
    return rec;
}

void destroy(void)
{
    free(farr);
}


