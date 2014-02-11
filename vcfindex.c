/* The MIT License

   Copyright (c) 2014 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Index bgzip compressed VCF/BCF files for random access.\n");
    fprintf(stderr, "Usage:   bcftools index [options] <in.bcf>|<in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -f, --force           overwrite index if it already exists\n");
    fprintf(stderr, "    -m, --min-shift INT   set the minimal interval size to 1<<INT [14]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Notes: The old tabix (.tbi) index can be invoked by setting -m0.\n");
    fprintf(stderr, "       Otherwise the new coordinate-sorted (.csi) index is created.\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfindex(int argc, char *argv[])
{
    int c, min_shift = 14, force = 0;

    static struct option loptions[] = 
    {
        {"help",0,0,'h'},
        {"force",0,0,'f'},
        {"min-shift",1,0,'m'},
        {0,0,0,0}
    };

    while ((c = getopt_long(argc, argv, "h?fm:", loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'f': force = 1; break;
            case 'm': min_shift = atoi(optarg); break;
            default: usage();
        }
    }
    if ( optind==argc ) usage();
    if (min_shift < 0 || min_shift > 30)
    {
        fprintf(stderr, "[E::%s] expected min_shift in range [0,30] (%d)\n", __func__, min_shift);
        return 1;
    }

    char *fname = argv[optind];
    int ftype = hts_file_type(fname);
    if (!ftype)
    {
        fprintf(stderr, "[E::%s] unknown filetype; expected .vcf.gz or .bcf\n", __func__);
        return 1;
    }

    if (!force)
    {
        char *fn;
        FILE *fp;
        fn = (char*)alloca(strlen(fname) + 5);
        strcat(strcpy(fn, fname), min_shift <= 0 ? ".tbi" : ".csi");
        if ((fp = fopen(fn, "rb")) != 0)
        {
            fclose(fp);
            fprintf(stderr, "[E::%s] the index file exists; use option '-f' to overwrite\n", __func__);
            return 1;
        }
    }

    if (ftype == FT_BCF_GZ)
    {
        if ( bcf_index_build(fname, min_shift) != 0 ) 
        {
            fprintf(stderr,"[E::%s] bcf_index_build failed: Is the BCF compressed?\n", __func__);
            return 1;
        }        
    }
    else if (ftype == FT_VCF_GZ)
    {
        if ( tbx_index_build(fname, min_shift, &tbx_conf_vcf) != 0 )
        {
            fprintf(stderr,"[E::%s] tbx_index_build failed: Is the file bgzip-compressed?\n", __func__);
            return 1;
        }
    }
    return 0;
}
