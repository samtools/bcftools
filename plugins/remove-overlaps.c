/*
    Copyright (C) 2017-2024 Genome Research Ltd.

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
#include <htslib/bgzf.h>
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
    int argc, region_is_file, target_is_file, output_type, verbose, nrm, ntot, clevel;
    int reverse;
    char **argv, *region, *target, *fname, *output_fname, *mark_expr, *mark_tag, *missing_expr;
    htsFile *out_fh;
    BGZF *fh_bgzf;
    bcf_hdr_t *hdr;
    bcf_srs_t *sr;
    char *index_fn;
    int write_index, record_cmd_line;
    kstring_t kstr;
}
args_t;

const char *about(void)
{
    return "Remove, list or mark overlapping variants\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: Remove, list or mark overlapping variants.\n"
        "\n"
        "Usage: bcftools +remove-overlaps [OPTIONS]\n"
        "Plugin options:\n"
        "   -M, --mark-tag TAG          Mark -m sites with INFO/TAG\n"
        "   -m, --mark EXPR             Mark (if also -M is present) or remove sites [overlap]\n"
        "                                   dup       .. all duplicate sites\n"
        "                                   overlap   .. overlapping sites\n"
        "                                   min(QUAL) .. sites with lowest QUAL until overlaps are resolved\n"
        "       --missing EXPR          Value to use for missing tags with -m 'min(QUAL)'\n"
        "                                   0   .. the default\n"
        "                                   DP  .. heuristics, scale maximum QUAL value proportionally to INFO/DP\n"
        "       --reverse               Apply the reverse logic, for example preserve duplicates instead of removing\n"
        "Standard options:\n"
        "   -e, --exclude EXPR          Exclude sites for which the expression is true\n"
        "   -i, --include EXPR          Include only sites for which the expression is true\n"
        "   -o, --output FILE                 Write output to the FILE [standard output]\n"
        "   -O, --output-type u|b|v|z|t[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level\n"
        "                                       t: plain list of sites (chr,pos), tz: compressed list [v]\n"
        "   -r, --regions REGION              Restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE           Restrict to regions listed in a file\n"
        "   -t, --targets REGION              Similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE           Similar to -R but streams rather than index-jumps\n"
        "       --no-version                  Do not append version and command line to the header\n"
        "   -W, --write-index[=FMT]           Automatically index the output files [off]\n"
        "\n";
}

// duh: FT_TAB_TEXT is 0 :-/
static int is_text(int flag)
{
    if ( flag==FT_TAB_TEXT || flag==FT_GZ ) return 1;
    return 0;
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

    char wmode[8];
    set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
    if ( is_text(args->output_type) )
    {
        args->fh_bgzf = bgzf_open(args->output_fname, args->output_type&FT_GZ ? "wg" : "wu");
    }
    else
    {
        args->out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
        if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));

        // todo: allow both INFO vs FILTER?
        if ( args->mark_tag )
        {
            int ret = bcf_hdr_printf(args->hdr, "##INFO=<ID=%s,Type=Flag,Number=0,Description=\"Marked by +remove-overlaps\">",args->mark_tag);
            if ( ret!=0 ) error("Error adding the header tag INFO/%s\n",args->mark_tag);
        }
        if ( args->record_cmd_line ) bcf_hdr_append_version(args->hdr, args->argc, args->argv, "bcftools_remove-overlaps");
        if ( bcf_hdr_write(args->out_fh, args->hdr)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
        if ( init_index2(args->out_fh,args->hdr,args->output_fname,&args->index_fn, args->write_index)<0 )
            error("Error: failed to initialise index for %s\n",args->output_fname);
    }

    args->vcfbuf = vcfbuf_init(args->hdr, 0);
    vcfbuf_set(args->vcfbuf,MARK,args->mark_expr);
    if ( args->missing_expr ) vcfbuf_set(args->vcfbuf,MARK_MISSING_EXPR,args->missing_expr);

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);
}
static void destroy_data(args_t *args)
{
    if ( args->filter )
        filter_destroy(args->filter);
    if ( args->write_index )
    {
        if ( bcf_idx_save(args->out_fh)<0 )
        {
            if ( hts_close(args->out_fh)!=0 ) error("Error: close failed .. %s\n", args->output_fname?args->output_fname:"stdout");
            error("Error: cannot write to index %s\n", args->index_fn);
        }
        free(args->index_fn);
    }
    if ( args->out_fh && hts_close(args->out_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname);
    if ( args->fh_bgzf && bgzf_close(args->fh_bgzf)!=0 ) error("[%s] Error: close failed .. %s\n",__func__,args->output_fname);
    vcfbuf_destroy(args->vcfbuf);
    bcf_sr_destroy(args->sr);
    free(args->kstr.s);
    free(args);
}
static void flush(args_t *args, int flush_all)
{
    bcf1_t *rec;
    while ( (rec=vcfbuf_flush(args->vcfbuf,flush_all)) )
    {
        int keep = vcfbuf_get_val(args->vcfbuf,int,MARK) ? 0 : 1;
        if ( args->reverse ) keep = keep ? 0 : 1;
        if ( !keep )
        {
            args->nrm++;
            if ( !args->mark_tag ) continue;
            bcf_update_info_flag(args->hdr,rec,args->mark_tag,NULL,1);
        }

        int ret;
        if ( is_text(args->output_type) )
        {
            args->kstr.l = 0;
            ksprintf(&args->kstr,"%s\t%"PRIhts_pos"\n",bcf_seqname(args->hdr,rec),rec->pos+1);
            if ( args->kstr.l && bgzf_write(args->fh_bgzf, args->kstr.s, args->kstr.l)!=args->kstr.l )
                error("Failed to write to %s\n", args->output_fname);
        }
        else
        {
            ret = bcf_write1(args->out_fh, args->hdr, rec);
            if ( ret!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
        }
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
    args->mark_expr = "overlap";
    args->clevel = -1;
    args->record_cmd_line = 1;
    static struct option loptions[] =
    {
        {"mark-tag",required_argument,NULL,'M'},
        {"mark",required_argument,NULL,'m'},
        {"reverse",no_argument,NULL,1},
        {"no-version",no_argument,NULL,2},
        {"missing",required_argument,NULL,3},
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"write-index",optional_argument,NULL,'W'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "m:M:r:R:t:T:o:O:i:e:dW::",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'm': args->mark_expr = optarg; break;
            case 'M': args->mark_tag = optarg; break;
            case  1 : args->reverse = 1; break;
            case  2 : args->record_cmd_line = 0; break;
            case  3 : args->missing_expr = optarg; break;
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
                          case 't': args->output_type = FT_TAB_TEXT; break;
                          default:
                          {
                              args->clevel = strtol(optarg,&tmp,10);
                              if ( *tmp || args->clevel<0 || args->clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                          }
                      }
                      if ( optarg[1]=='z' )
                      {
                          optarg++;
                          args->output_type |= FT_GZ;
                      }
                      if ( optarg[1] )
                      {
                          args->clevel = strtol(optarg+1,&tmp,10);
                          if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                      }
                      break;
            case 'W':
                if (!(args->write_index = write_index_parse(optarg)))
                    error("Unsupported index format '%s'\n", optarg);
                break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( args->write_index && is_text(args->output_type) ) error("Cannot index text output\n");
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


