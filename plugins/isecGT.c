/*
    Copyright (C) 2016-2023 Genome Research Ltd.

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
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <stdarg.h>
#include <stdint.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <inttypes.h>
#include <unistd.h>
#include <errno.h>
#include "bcftools.h"
#include "smpl_ilist.h"

typedef struct
{
    int argc, output_type, regions_is_file, targets_is_file, clevel;
    char **argv, *output_fname, *regions_list, *targets_list;
    int32_t *arr_a, narr_a, *arr_b, narr_b;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr_a, *hdr_b;
    htsFile *out_fh;
    char *index_fn;
    int write_index;
}
args_t;

const char *about(void)
{
    return "Compare two files and set non-identical genotypes to missing.\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: Compare two files and set non-identical genotypes in the first file to missing.\n"
        "\n"
        "Usage: bcftools +isecGT <A.bcf> <B.bcf> [Plugin Options]\n"
        "Plugin options:\n"
        "   -o, --output FILE               Write output to a file [standard output]\n"
        "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n"
        "   -r, --regions REGION            Restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE         Restrict to regions listed in a file\n"
        "   -t, --targets REGION            Similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE         Similar to -R but streams rather than index-jumps\n"
        "   -W, --write-index[=FMT]         Automatically index the output files [off]\n"
        "\n";
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->output_fname = "-";
    args->output_type = FT_VCF;
    args->clevel = -1;
    static struct option loptions[] =
    {
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"targets",required_argument,NULL,'t'},
        {"targets-file",required_argument,NULL,'T'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"write-index",optional_argument,NULL,'W'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "o:O:r:R:t:T:W::",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'o': args->output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': args->output_type = FT_BCF_GZ; break;
                          case 'u': args->output_type = FT_BCF; break;
                          case 'z': args->output_type = FT_VCF_GZ; break;
                          case 'v': args->output_type = FT_VCF; break;
                          default:
                          {
                              args->clevel = strtol(optarg,&tmp,10);
                              if ( *tmp || args->clevel<0 || args->clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                          }
                      }
                      if ( optarg[1] )
                      {
                          args->clevel = strtol(optarg+1,&tmp,10);
                          if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                      }
                      break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; args->regions_is_file = 1; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; args->targets_is_file = 1; break;
            case 'W':
                if (!(args->write_index = write_index_parse(optarg)))
                    error("Unsupported index format '%s'\n", optarg);
                break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }

    if ( optind+2!=argc ) error("%s",usage_text());

    args->sr = bcf_sr_init();
    args->sr->require_index = 1;
    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->sr, args->regions_list, args->regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(args->sr, args->targets_list, args->targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
        args->sr->collapse |= COLLAPSE_BOTH;
    }
    if ( !bcf_sr_add_reader(args->sr,argv[optind]) ) error("Error opening %s: %s\n", argv[optind],bcf_sr_strerror(args->sr->errnum));
    if ( !bcf_sr_add_reader(args->sr,argv[optind+1]) ) error("Error opening %s: %s\n", argv[optind+1],bcf_sr_strerror(args->sr->errnum));
    args->hdr_a = bcf_sr_get_header(args->sr,0);
    args->hdr_b = bcf_sr_get_header(args->sr,1);
    smpl_ilist_t *smpl = smpl_ilist_map(args->hdr_a, args->hdr_b, SMPL_STRICT);
    char wmode[8];
    set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
    args->out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( bcf_hdr_write(args->out_fh, args->hdr_a)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
    if ( init_index2(args->out_fh,args->hdr_a,args->output_fname,
                     &args->index_fn, args->write_index)<0 )
      error("Error: failed to initialise index for %s\n",args->output_fname);

    while ( bcf_sr_next_line(args->sr) )
    {
        if ( !bcf_sr_has_line(args->sr,0) ) continue;
        if ( !bcf_sr_has_line(args->sr,1) )
        {
            if ( bcf_write(args->out_fh, args->hdr_a, bcf_sr_get_line(args->sr,0))!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
            continue;
        }

        bcf1_t *line_a = bcf_sr_get_line(args->sr,0);
        bcf1_t *line_b = bcf_sr_get_line(args->sr,1);
        int ngt_a = bcf_get_genotypes(args->hdr_a, line_a, &args->arr_a, &args->narr_a);
        int ngt_b = bcf_get_genotypes(args->hdr_b, line_b, &args->arr_b, &args->narr_b);
        assert( ngt_a==ngt_b );     // todo
        ngt_a /= smpl->n;
        ngt_b /= smpl->n;
        int i, j, dirty = 0;
        for (i=0; i<smpl->n; i++)
        {
            int32_t *a = args->arr_a + i*ngt_a;
            int32_t *b = args->arr_b + smpl->idx[i]*ngt_b;
            for (j=0; j<ngt_a; j++)
                if ( a[j]!=b[j] ) break;
            if ( j<ngt_a )
            {
                dirty = 1;
                for (j=0; j<ngt_a; j++) a[j] = bcf_gt_missing;
            }
        }
        if ( dirty ) bcf_update_genotypes(args->hdr_a, line_a, args->arr_a, ngt_a*smpl->n);
        if ( bcf_write(args->out_fh, args->hdr_a, line_a)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
    }
    if ( args->write_index )
    {
        if ( bcf_idx_save(args->out_fh)<0 )
        {
            if ( hts_close(args->out_fh)!=0 ) error("Error: close failed .. %s\n", args->output_fname?args->output_fname:"stdout");
            error("Error: cannot write to index %s\n", args->index_fn);
        }
        free(args->index_fn);
    }
    if ( hts_close(args->out_fh)!=0 ) error("Close failed: %s\n",args->output_fname);
    smpl_ilist_destroy(smpl);
    bcf_sr_destroy(args->sr);
    free(args->arr_a);
    free(args->arr_b);
    free(args);
    return 0;
}

