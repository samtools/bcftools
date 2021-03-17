/* 
    Copyright (C) 2017-2021 Genome Research Ltd.

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
#include <strings.h>
#include <getopt.h>
#include <stdarg.h>
#include <unistd.h>
#include <stdint.h>
#include <errno.h>
#include <inttypes.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "vcfbuf.h"
#include "filter.h"

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct
{
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)
    vcfbuf_t *vcfbuf;
    int argc, region_is_file, target_is_file, output_type, verbose, nrm, ntot, print_overlaps, rmdup;
    char **argv, *region, *target, *fname, *output_fname;
    htsFile *out_fh;
    bcf_hdr_t *hdr;
    bcf_srs_t *sr;
}
args_t;

const char *about(void)
{
    return "Remove overlapping variants\n";
}

static const char *usage_text(void)
{
    return 
        "\n"
        "About: Remove overlapping variants.\n"
        "\n"
        "Usage: bcftools +remove-overlaps [Options]\n"
        "Plugin options:\n"
        "   -d, --rm-dup                    remove only duplicate sites and remove them completely\n"
        "   -p, --print-overlaps            do the opposite and print only overlapping sites\n"
        "   -v, --verbose                   print a list of removed sites\n"
        "Standard options:\n"
        "   -e, --exclude EXPR              exclude sites for which the expression is true\n"
        "   -i, --include EXPR              include only sites for which the expression is true\n"
        "   -o, --output FILE               write output to the FILE [standard output]\n"
        "   -O, --output-type b|u|z|v       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n"
        "   -r, --regions REGION            restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE         restrict to regions listed in a file\n"
        "   -t, --targets REGION            similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE         similar to -R but streams rather than index-jumps\n"
        "\n";
}

static void init_data(args_t *args)
{
    args->sr = bcf_sr_init();
    if ( args->region )
    {
        args->sr->require_index = 1;
        if ( bcf_sr_set_regions(args->sr, args->region, args->region_is_file)<0 ) error("Failed to read the regions: %s\n",args->region);
    }
    if ( args->target && bcf_sr_set_targets(args->sr, args->target, args->target_is_file, 0)<0 ) error("Failed to read the targets: %s\n",args->target);
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr = bcf_sr_get_header(args->sr,0);

    args->out_fh = hts_open(args->output_fname,hts_bcf_wmode2(args->output_type,args->output_fname));
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( bcf_hdr_write(args->out_fh, args->hdr)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);

    args->vcfbuf = vcfbuf_init(args->hdr, 0);
    if ( args->rmdup )
        vcfbuf_set_opt(args->vcfbuf,int,VCFBUF_RMDUP,1)
    else
        vcfbuf_set_opt(args->vcfbuf,int,VCFBUF_OVERLAP_WIN,1)

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);
}
static void destroy_data(args_t *args)
{
    if ( args->filter )
        filter_destroy(args->filter);
    if ( hts_close(args->out_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname);
    vcfbuf_destroy(args->vcfbuf);
    bcf_sr_destroy(args->sr);
    free(args);
}
static void flush(args_t *args, int flush_all)
{
    int nbuf = vcfbuf_nsites(args->vcfbuf);
    bcf1_t *rec;
    while ( (rec = vcfbuf_flush(args->vcfbuf, flush_all)) )
    {
        if ( nbuf>2 || (nbuf>1 && flush_all) )
        {
            args->nrm++;
            if ( args->verbose ) printf("%s\t%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
            if ( args->print_overlaps && bcf_write1(args->out_fh, args->hdr, rec)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
            continue;     // skip overlapping variants
        }
        if ( !args->print_overlaps && bcf_write1(args->out_fh, args->hdr, rec)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
    }
}
static void process(args_t *args)
{
    args->ntot++;
    bcf1_t *rec = bcf_sr_get_line(args->sr,0);
    if ( args->filter )
    {
        int ret  = filter_test(args->filter, rec, NULL);
        if ( args->filter_logic==FLT_INCLUDE ) { if ( !ret ) return; }
        else if ( ret ) return;
    }
    bcf_sr_t *sr = bcf_sr_get_reader(args->sr, 0);
    sr->buffer[0] = vcfbuf_push(args->vcfbuf, rec);
    flush(args,0);
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_type  = FT_VCF;
    args->output_fname = "-";
    static struct option loptions[] =
    {
        {"rm-dup",no_argument,NULL,'d'},
        {"print-overlaps",no_argument,NULL,'p'},
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"verbose",no_argument,NULL,'v'},
        {NULL,0,NULL,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "r:R:t:T:o:O:i:e:vpd",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'd': args->rmdup = 1; break;
            case 'p': args->print_overlaps = 1; break;
            case 'v': args->verbose = 1; break;
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'T': args->target_is_file = 1; // fall-through
            case 't': args->target = optarg; break; 
            case 'R': args->region_is_file = 1; // fall-through
            case 'r': args->region = optarg; break; 
            case 'o': args->output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': args->output_type = FT_BCF_GZ; break;
                          case 'u': args->output_type = FT_BCF; break;
                          case 'z': args->output_type = FT_VCF_GZ; break;
                          case 'v': args->output_type = FT_VCF; break;
                          default: error("The output type \"%s\" not recognised\n", optarg);
                      }
                      break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error("%s",usage_text()); }
    }
    else if ( optind+1!=argc ) error("%s",usage_text());
    else args->fname = argv[optind];

    init_data(args);
    
    while ( bcf_sr_next_line(args->sr) ) process(args);
    flush(args,1);

    fprintf(stderr,"Processed/Removed\t%d\t%d\n",args->ntot,args->nrm);

    destroy_data(args);
    return 0;
}


