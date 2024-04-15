/*
    Copyright (C) 2017-2023 Genome Research Ltd.

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
/*
    Compress gVCF file by resizing gVCF blocks according to specified criteria.
*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <getopt.h>
#include <stdarg.h>
#include <unistd.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/synced_bcf_reader.h>
#include "bcftools.h"
#include "filter.h"

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define GQ_KEY_NONE NULL
#define GQ_KEY_GQ   "GQ"
#define GQ_KEY_RGQ  "RGQ"

typedef struct
{
    int32_t end, min_dp, gq, pl[3], grp;
    char *gq_key;
    bcf1_t *rec;
}
block_t;
typedef struct
{
    char *expr;     // expression
    int flt_id;     // filter id, -1 for PASS
    filter_t *flt;  // filter
}
grp_t;
typedef struct
{
    filter_t *filter;
    char *filter_str;
    int filter_logic;
    block_t gvcf;
    htsFile *fh_out;
    int ngrp;
    grp_t *grp;
    char *group_by;
    int argc, region_is_file, target_is_file, output_type, trim_alts, clevel;
    int32_t *tmpi, mtmpi, mean_min_dp_reported;
    char **argv, *region, *target, *fname, *output_fname, *keep_tags;
    bcf_hdr_t *hdr_in, *hdr_out;
    bcf_srs_t *sr;
    char *index_fn;
    int write_index;
}
args_t;

const char *about(void)
{
    return "Compress gVCF file by resizing gVCF blocks according to specified criteria.\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: Compress gVCF file by resizing gVCF blocks according to specified criteria.\n"
        "\n"
        "Usage: bcftools +gvcfz [Options]\n"
        "Plugin options:\n"
        "   -a, --trim-alt-alleles          Trim alternate alleles not seen in the genotypes\n"
        "   -e, --exclude <expr>            Exclude sites for which the expression is true\n"
        "   -i, --include <expr>            Include sites for which the expression is true\n"
        "   -g, --group-by EXPR             Group gVCF blocks according to the expression\n"
        "   -o, --output FILE               Write gVCF output to the FILE\n"
        "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n"
        "   -W, --write-index[=FMT]         Automatically index the output files [off]\n"
        "Examples:\n"
        "   # Compress blocks by GQ and DP. Multiple blocks separated by a semicolon can be defined\n"
        "   bcftools +gvcfz input.bcf -g'PASS:GQ>60 & DP<20; PASS:GQ>40 & DP<15; Flt1:QG>20; Flt2:-'\n"
        "\n"
        "   # Compress all non-reference sites into a single block, remove unused alternate alleles\n"
        "   bcftools +gvcfz input.bcf -a -g'PASS:GT!=\"alt\"'\n"
        "\n";
}

static void init_groups(args_t *args)
{
    args->hdr_out = bcf_hdr_dup(args->hdr_in);
    bcf_hdr_printf(args->hdr_out, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">");

    // avoid nested double quotes in FILTER description
    char *hdr_str = strdup(args->group_by);
    char *tmp = hdr_str;
    while (*tmp)
    {
        if ( *tmp=='"' ) *tmp = '\'';
        tmp++;
    }

    char *rmme_str = strdup(args->group_by), *beg = rmme_str;
    while ( *beg )
    {
        while ( *beg && isspace(*beg) ) beg++;
        if ( !beg ) break;
        char *end = beg;
        while ( *end && *end!=':' ) end++;
        if ( *end!=':' ) error("Could not parse the expression: \"%s\"\n", args->group_by);
        *end = 0;
        char *flt = beg;
        beg = ++end;
        while ( *end && *end!=';' ) end++;
        char tmp = *end; *end = 0;
        if ( strcmp(flt,"PASS") )
        {
            bcf_hdr_printf(args->hdr_out, "##FILTER=<ID=%s,Description=\"%s\">", flt, hdr_str);
            if (bcf_hdr_sync(args->hdr_out) < 0)
                error_errno("[%s] Failed to update header", __func__);
        }
        args->ngrp++;
        args->grp = (grp_t*) realloc(args->grp,sizeof(grp_t)*args->ngrp);
        grp_t *grp = args->grp + args->ngrp - 1;
        grp->expr = strdup(beg);
        grp->flt_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, flt);
        if ( !bcf_hdr_idinfo_exists(args->hdr_out, BCF_HL_FLT, grp->flt_id) ) error("Could not initialize the filter \"%s\"\n", flt);
        if ( !strcmp(flt,"PASS") ) grp->flt_id = -1;

        // remove trailing spaces
        beg = grp->expr + strlen(grp->expr); while ( beg >= grp->expr && isspace(*beg) ) { *beg = 0; beg--; }
        beg = grp->expr; while ( *beg && isspace(*beg) ) beg++;

        grp->flt = strcmp("-",beg) ? filter_init(args->hdr_in, grp->expr) : NULL;

        if ( !tmp ) break;
        beg = end + 1;
    }
    free(rmme_str);
    free(hdr_str);
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->ngrp; i++)
    {
        if ( args->grp[i].flt ) filter_destroy(args->grp[i].flt);
        free(args->grp[i].expr);
    }
    free(args->grp);

    if ( args->filter ) filter_destroy(args->filter);
    if ( args->write_index )
    {
        if ( bcf_idx_save(args->fh_out)<0 )
        {
            if ( hts_close(args->fh_out)!=0 ) error("Error: close failed .. %s\n", args->output_fname?args->output_fname:"stdout");
            error("Error: cannot write to index %s\n", args->index_fn);
        }
        free(args->index_fn);
    }
    if ( hts_close(args->fh_out)!=0 ) error("failed to close %s\n", args->output_fname);

    bcf_sr_destroy(args->sr);
    if ( args->hdr_out ) bcf_hdr_destroy(args->hdr_out);
    if ( args->gvcf.rec ) bcf_destroy(args->gvcf.rec);
    free(args->tmpi);
    free(args);
}

static void flush_block(args_t *args, bcf1_t *rec)
{
    block_t *gvcf = &args->gvcf;
    if ( gvcf->grp < 0 ) return;
    if ( rec && gvcf->end - 1 >= rec->pos ) gvcf->end = rec->pos; // NB: end is 1-based, rec->pos is 0-based

    if ( gvcf->rec->pos+1 < gvcf->end && bcf_update_info_int32(args->hdr_out,gvcf->rec,"END",&gvcf->end,1) != 0 )
        error("Could not update INFO/END at %s:%"PRId64"\n", bcf_seqname(args->hdr_out,gvcf->rec),(int64_t) gvcf->rec->pos+1);
    if ( bcf_update_format_int32(args->hdr_out,gvcf->rec,"DP",&gvcf->min_dp,1) != 0 )
        error("Could not update FORMAT/DP at %s:%"PRId64"\n", bcf_seqname(args->hdr_out,gvcf->rec),(int64_t) gvcf->rec->pos+1);
    if ( gvcf->gq_key )
    {
        if ( bcf_update_format_int32(args->hdr_out,gvcf->rec,gvcf->gq_key,&gvcf->gq,1) != 0 )
            error("Could not update FORMAT/%s at %s:%"PRId64"\n", gvcf->gq_key, bcf_seqname(args->hdr_out,gvcf->rec),(int64_t) gvcf->rec->pos+1);
    }
    if ( gvcf->pl[0] >=0 )
    {
        if ( bcf_update_format_int32(args->hdr_out,gvcf->rec,"PL",&gvcf->pl,3) != 0 )
            error("Could not update FORMAT/PL at %s:%"PRId64"\n", bcf_seqname(args->hdr_out,gvcf->rec),(int64_t) gvcf->rec->pos+1);
    }
    if ( gvcf->grp < args->ngrp && args->grp[gvcf->grp].flt_id >= 0 )
        bcf_add_filter(args->hdr_out, gvcf->rec, args->grp[gvcf->grp].flt_id);

    if ( bcf_write(args->fh_out, args->hdr_out, gvcf->rec)!=0 ) error("Failed to write the header\n");

    gvcf->grp = -1;
}
static void process_gvcf(args_t *args)
{
    bcf1_t *rec = bcf_sr_get_line(args->sr,0);

    if ( args->filter )
    {
        int pass = filter_test(args->filter, rec, NULL);
        if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
        if ( !pass ) return;
    }

    if ( rec->n_allele > 2 || (rec->n_allele == 2 && strcmp("<NON_REF>",rec->d.allele[1]) && strcmp("<*>",rec->d.allele[1])) )
    {
        if ( args->trim_alts )
        {
            bcf_unpack(rec, BCF_UN_ALL);
            if ( bcf_trim_alleles(args->hdr_in, rec)<0 )
                error("Error: Could not trim alleles at %s:%"PRId64"\n", bcf_seqname(args->hdr_in, rec),(int64_t)  rec->pos+1);

            // trim the ref allele if necessary
            if ( rec->d.allele[0][1] )
            {
                rec->d.allele[0][1] = 0;
                bcf_update_alleles(args->hdr_in, rec, (const char**)rec->d.allele, 1);
            }

        }
        if ( rec->n_allele > 2 || (rec->n_allele == 2 && strcmp("<NON_REF>",rec->d.allele[1]) && strcmp("<*>",rec->d.allele[1])) )
        {
            // not a gvcf block
            flush_block(args, rec);
            if ( bcf_write(args->fh_out, args->hdr_out, rec)!=0 ) error("Failed to write\n");
            return;
        }
    }

    int ret = bcf_get_info_int32(args->hdr_in,rec,"END",&args->tmpi,&args->mtmpi);
    int32_t end = ret==1 ? args->tmpi[0] : rec->pos + 1;

    char *gq_key = GQ_KEY_GQ;
    ret = bcf_get_format_int32(args->hdr_in,rec,gq_key,&args->tmpi,&args->mtmpi);
    if ( ret!=1 )
    {
        gq_key = GQ_KEY_RGQ;
        if ( ret<1 ) ret = bcf_get_format_int32(args->hdr_in,rec,gq_key,&args->tmpi,&args->mtmpi);
        if ( ret!=1 ) gq_key = GQ_KEY_NONE;
    }
    int32_t gq = ret==1 ? args->tmpi[0] : 0;

    int32_t min_dp = 0;
    if ( bcf_get_format_int32(args->hdr_in,rec,"MIN_DP",&args->tmpi,&args->mtmpi)==1 )
        min_dp = args->tmpi[0];
    else if ( bcf_get_format_int32(args->hdr_in,rec,"DP",&args->tmpi,&args->mtmpi)==1 )
        min_dp = args->tmpi[0];
    else
        error("Expected one FORMAT/MIN_DP or FORMAT/DP value at %s:%"PRId64"\n", bcf_seqname(args->hdr_in,rec),(int64_t) rec->pos+1);

    int32_t pl[3] = {-1,-1,-1};
    ret = bcf_get_format_int32(args->hdr_in,rec,"PL",&args->tmpi,&args->mtmpi);
    if ( ret>3 ) error("Expected three FORMAT/PL values at %s:%"PRId64"\n", bcf_seqname(args->hdr_in,rec),(int64_t) rec->pos+1);
    else if ( ret==3 )
    {
        pl[0] = args->tmpi[0];
        pl[1] = args->tmpi[1];
        pl[2] = args->tmpi[2];
    }

    int i;
    for (i=0; i<args->ngrp; i++)
        if ( !args->grp[i].flt || filter_test(args->grp[i].flt, rec, NULL)==1 ) break;

    if ( args->gvcf.grp != i ) flush_block(args, rec);      // new block
    if ( args->gvcf.grp >= 0 && args->gvcf.rec->rid != rec->rid ) flush_block(args, NULL);  // new chromosome

    if ( args->gvcf.grp >= 0 ) // extend an existing block
    {
        if ( args->gvcf.end < end ) args->gvcf.end = end;
        if ( args->gvcf.gq_key!=GQ_KEY_NONE && gq_key!=GQ_KEY_NONE && args->gvcf.gq > gq ) args->gvcf.gq = gq;
        if ( args->gvcf.min_dp > min_dp ) args->gvcf.min_dp = min_dp;
        if ( args->gvcf.pl[0] > pl[0] ) args->gvcf.pl[0] = pl[0];
        if ( args->gvcf.pl[1] > pl[1] ) args->gvcf.pl[1] = pl[1];
        if ( args->gvcf.pl[2] > pl[2] ) args->gvcf.pl[2] = pl[2];
        return;
    }

    // start a new block
    args->gvcf.rec = bcf_copy(args->gvcf.rec, rec);
    args->gvcf.grp = i;
    args->gvcf.min_dp   = min_dp;
    args->gvcf.end      = end;
    args->gvcf.pl[0]    = pl[0];
    args->gvcf.pl[1]    = pl[1];
    args->gvcf.pl[2]    = pl[2];
    args->gvcf.gq_key   = gq_key;
    if ( gq_key!=GQ_KEY_NONE ) args->gvcf.gq = gq;
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_type  = FT_VCF;
    args->output_fname = "-";
    args->clevel = -1;
    static struct option loptions[] =
    {
        {"trim-alt-alleles",required_argument,0,'a'},
        {"include",required_argument,0,'i'},
        {"exclude",required_argument,0,'e'},
        {"group-by",required_argument,NULL,'g'},
        {"stats",required_argument,NULL,'s'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"write-index",optional_argument,NULL,'W'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "vr:R:t:T:o:O:g:i:e:aW::",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'a': args->trim_alts = 1; break;
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'g': args->group_by = optarg; break;
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
            case 'W':
                if (!(args->write_index = write_index_parse(optarg)))
                    error("Unsupported index format '%s'\n", optarg);
                break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error("%s", usage_text()); }
    }
    else if ( optind+1!=argc ) error("%s", usage_text());
    else args->fname = argv[optind];

    if ( !args->group_by ) error("Missing the -g option\n");

    args->gvcf.rec = bcf_init();
    args->gvcf.grp = -1;            // the block is inactive
    args->sr = bcf_sr_init();
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr_in = bcf_sr_get_header(args->sr,0);
    if ( args->filter_str )
        args->filter = filter_init(args->hdr_in, args->filter_str);
    init_groups(args);
    char wmode[8];
    set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
    args->fh_out = hts_open(args->output_fname ? args->output_fname : "-", wmode);
    if ( bcf_hdr_write(args->fh_out, args->hdr_out)!=0 ) error("Failed to write the header\n");
    if ( init_index2(args->fh_out,args->hdr_out,args->output_fname,
                     &args->index_fn, args->write_index)<0 )
        error("Error: failed to initialise index for %s\n",args->output_fname);
    while ( bcf_sr_next_line(args->sr) ) process_gvcf(args);
    flush_block(args, NULL);

    destroy_data(args);
    return 0;
}


