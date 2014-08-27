/*  vcfindex.c -- Index bgzip compressed VCF/BCF files for random access.

    Copyright (C) 2014 Genome Research Ltd.

    Author: Shane McCarthy <sm15@sanger.ac.uk>

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
#include <unistd.h>
#include <getopt.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <sys/stat.h>

#define BCF_LIDX_SHIFT    14

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Index bgzip compressed VCF/BCF files for random access.\n");
    fprintf(stderr, "Usage:   bcftools index [options] <in.bcf>|<in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -c, --csi             generate CSI-format index for VCF/BCF files [default]\n");
    fprintf(stderr, "    -f, --force           overwrite index if it already exists\n");
    fprintf(stderr, "    -m, --min-shift INT   set minimal interval size for CSI indices to 2^INT [14]\n");
    fprintf(stderr, "    -t, --tbi             generate TBI-format index for VCF files\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfindex(int argc, char *argv[])
{
    int c, force = 0, tbi = 0;
    int min_shift = BCF_LIDX_SHIFT;

    static struct option loptions[] =
    {
        {"csi",no_argument,NULL,'c'},
        {"tbi",no_argument,NULL,'t'},
        {"force",no_argument,NULL,'f'},
        {"min-shift",required_argument,NULL,'m'},
        {NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "ctfm:", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'c': tbi = 0; break;
            case 't': tbi = 1; min_shift = 0; break;
            case 'f': force = 1; break;
            case 'm': min_shift = atoi(optarg); break;
            default: usage();
        }
    }
    if ( optind==argc ) usage();
    if (tbi && min_shift>0)
    {
        fprintf(stderr, "[E::%s] min-shift option only expected for CSI indices \n", __func__);
        return 1;
    }
    if (min_shift < 0 || min_shift > 30)
    {
        fprintf(stderr, "[E::%s] expected min_shift in range [0,30] (%d)\n", __func__, min_shift);
        return 1;
    }

    char *fname = argv[optind];
    int ftype = hts_file_type(fname);
    if (!ftype || (ftype != FT_BCF_GZ && ftype != FT_VCF_GZ))
    {
        fprintf(stderr, "[E::%s] unknown filetype; expected bgzip compressed VCF or BCF\n", __func__);
        if (!(ftype & FT_GZ))
            fprintf(stderr, "[E::%s] was the VCF/BCF compressed with bgzip?\n", __func__);
        return 1;
    }
    if (tbi && ftype == FT_BCF_GZ)
    {
        fprintf(stderr, "[Warning] TBI-index does not work for BCF files. Generating CSI instead.\n");
        tbi = 0; min_shift = BCF_LIDX_SHIFT;
    }
    if (min_shift == 0 && ftype == FT_BCF_GZ)
    {
        fprintf(stderr, "[E::%s] Require min_shift>0 for BCF files.\n", __func__);
        return 1;
    }
    if (!tbi && ftype == FT_VCF_GZ && min_shift == 0)
    {
        fprintf(stderr, "[Warning] min-shift set to 0 for VCF file. Generating TBI file.\n");
        tbi = 1;
    }

    if (!force)
    {
        // Before complaining about existing index, check if the VCF file isn't newer.
        char *idx_fname = (char*)alloca(strlen(fname) + 5);
        strcat(strcpy(idx_fname, fname), tbi ? ".tbi" : ".csi");
        struct stat stat_tbi, stat_file;
        if ( stat(idx_fname, &stat_tbi)==0 )
        {
            stat(fname, &stat_file);
            if ( stat_file.st_mtime <= stat_tbi.st_mtime )
            {
                fprintf(stderr,"[E::%s] the index file exists. Please use '-f' to overwrite.\n", __func__);
                return 1;
            }
        }
    }

    if (ftype == FT_BCF_GZ)
    {
        if ( bcf_index_build(fname, min_shift) != 0 )
        {
            fprintf(stderr,"[E::%s] bcf_index_build failed for %s\n", __func__, fname);
            return 1;
        }
    }
    else
    {
        if ( tbx_index_build(fname, min_shift, &tbx_conf_vcf) != 0 )
        {
            fprintf(stderr,"[E::%s] tbx_index_build failed for %s\n", __func__, fname);
            return 1;
        }
    }
    return 0;
}
