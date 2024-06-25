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
/*
    Prune sites by missingness, LD

    See calc_ld() in vcfbuf.c for the actual LD calculation

*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <getopt.h>
#include <stdarg.h>
#include <unistd.h>
#include <stdint.h>
#include <errno.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <htslib/hts_os.h>
#include <assert.h>
#include <time.h>
#include "bcftools.h"
#include "vcfbuf.h"
#include "filter.h"

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define LD_ANNOTATE 1
#define LD_SET_MAX  2

typedef struct
{
    filter_t *filter;
    char *filter_str, *af_tag;
    int filter_logic;   // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)
    vcfbuf_t *vcfbuf;
    double ld_max[VCFBUF_LD_N];
    int ld_max_set[VCFBUF_LD_N];
    char *ld_annot[VCFBUF_LD_N], *ld_annot_pos[VCFBUF_LD_N];
    int ld_mask;
    int argc, region_is_file, target_is_file, output_type, ld_filter_id, rand_missing, nsites, ld_win, rseed, clevel;
    char *nsites_mode;
    int keep_sites;
    char **argv, *region, *target, *fname, *output_fname, *ld_filter;
    htsFile *out_fh;
    bcf_hdr_t *hdr;
    bcf_srs_t *sr;
    char *index_fn;
    int write_index;
}
args_t;

const char *about(void)
{
    return "Prune sites by missingness, linkage disequilibrium\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: Prune sites by missingness or linkage disequilibrium.\n"
        "\n"
        "Usage: bcftools +prune [Options]\n"
        "Plugin options:\n"
        "       --AF-tag STR                Use this tag with -n to determine allele frequency\n"
        "   -a, --annotate r2,LD            Add position of an upstream record with the biggest r2/LD value\n"
        "   -e, --exclude EXPR              Exclude sites for which the expression is true\n"
        "   -f, --set-filter STR            Apply soft filter STR instead of discarding the site (only with -m)\n"
        "   -i, --include EXPR              Include only sites for which the expression is true\n"
        "   -k, --keep-sites                Leave sites filtered by -i/-e unchanged instead of discarding them\n"
        "   -m, --max [r2|LD=]FLOAT         Remove sites with r2 or Lewontin's D bigger than FLOAT within the -w window\n"
        "   -n, --nsites-per-win N          Keep at most N sites in the -w window. See also -N, --nsites-per-win-mode\n"
        "   -N, --nsites-per-win-mode STR   Keep sites with biggest AF (\"maxAF\"); sites that come first (\"1st\"); pick randomly (\"rand\") [maxAF]\n"
        "   -o, --output FILE               Write output to the FILE [standard output]\n"
        "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n"
        "       --random-seed INT           Use the provided random seed for reproducibility\n"
        "       --randomize-missing         Replace missing data with randomly assigned genotype based on site's allele frequency\n"
        "   -r, --regions REGION            Restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE         Restrict to regions listed in a file\n"
        "   -t, --targets REGION            Similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE         Similar to -R but streams rather than index-jumps\n"
        "   -w, --window INT[bp|kb|Mb]      The window size of INT sites or INT bp/kb/Mb for the -n/-l options [100kb]\n"
        "   -W, --write-index[=FMT]         Automatically index the output files [off]\n"
        "Examples:\n"
        "   # Discard records with r2 bigger than 0.6 in a window of 1000 sites\n"
        "   bcftools +prune -m 0.6 -w 1000 input.bcf -Ob -o output.bcf\n"
        "\n"
        "   # Set FILTER (but do not discard) records with r2 bigger than 0.4 in the default window of 100kb\n"
        "   bcftools +prune -m 0.4 -f MAX_R2 input.bcf -Ob -o output.bcf\n"
        "\n"
        "   # Filter nothing, only annotate each site with r2 and LD to the previous site\n"
        "   bcftools +prune -w 1 -a r2,LD -f . input.bcf -Ob -o output.bcf\n"
        "\n"
        "   # Annotate INFO field of all records with maximum r2 in a window of 1000 sites\n"
        "   bcftools +prune -m 0.6 -w 1000 -f MAX_R2 input.bcf -Ob -o output.bcf\n"
        "\n"
        "   # Discard records with r2 bigger than 0.6, first removing records with more than 2% of genotypes missing\n"
        "   bcftools +prune -m 0.6 -e'F_MISSING>=0.02' input.bcf -Ob -o output.bcf\n"
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

    char wmode[8];
    set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
    args->out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));

    if ( args->ld_filter && strcmp(".",args->ld_filter) )
    {
        kstring_t str = {0,0,0};
        if ( args->ld_max_set[VCFBUF_LD_IDX_R2] )
        {
            kputs("R2 bigger than ",&str);
            kputd(args->ld_max[VCFBUF_LD_IDX_R2],&str);
        }
        if ( args->ld_max_set[VCFBUF_LD_IDX_LD] )
        {
            if ( str.l ) kputs(" or ",&str);
            kputs("LD bigger than ",&str);
            kputd(args->ld_max[VCFBUF_LD_IDX_LD],&str);
        }
        if ( args->ld_max_set[VCFBUF_LD_IDX_HD] )
        {
            if ( str.l ) kputs(" or ",&str);
            kputs("HD bigger than ",&str);
            kputd(args->ld_max[VCFBUF_LD_IDX_HD],&str);
        }
        bcf_hdr_printf(args->hdr,"##FILTER=<ID=%s,Description=\"An upstream site within %d%s with %s\">",args->ld_filter,
                args->ld_win < 0 ? -args->ld_win/1000 : args->ld_win,
                args->ld_win < 0 ? "kb" : " sites",
                str.s);
        free(str.s);
    }
    if ( args->ld_mask & LD_ANNOTATE )
    {
        if ( args->ld_annot[VCFBUF_LD_IDX_R2] )
        {
            bcf_hdr_printf(args->hdr,"##INFO=<ID=%s,Number=1,Type=Float,Description=\"Pairwise r2 with the %s site\">",args->ld_annot[VCFBUF_LD_IDX_R2],args->ld_annot_pos[VCFBUF_LD_IDX_R2]);
            bcf_hdr_printf(args->hdr,"##INFO=<ID=%s,Number=1,Type=Integer,Description=\"The position of the site for which %s was calculated\">",args->ld_annot_pos[VCFBUF_LD_IDX_R2],args->ld_annot[VCFBUF_LD_IDX_R2]);
        }
        if ( args->ld_annot[VCFBUF_LD_IDX_LD] )
        {
            bcf_hdr_printf(args->hdr,"##INFO=<ID=%s,Number=1,Type=Float,Description=\"Pairwise Lewontin's D' (PMID:19433632) with the %s site\">",args->ld_annot[VCFBUF_LD_IDX_LD],args->ld_annot_pos[VCFBUF_LD_IDX_LD]);
            bcf_hdr_printf(args->hdr,"##INFO=<ID=%s,Number=1,Type=Integer,Description=\"The position of the site for which %s was calculated\">",args->ld_annot_pos[VCFBUF_LD_IDX_LD],args->ld_annot[VCFBUF_LD_IDX_LD]);
        }
        if ( args->ld_annot[VCFBUF_LD_IDX_HD] )
        {
            bcf_hdr_printf(args->hdr,"##INFO=<ID=%s,Number=1,Type=Float,Description=\"Pairwise Ragsdale's \\hat{D} (PMID:31697386) with the %s site\">",args->ld_annot[VCFBUF_LD_IDX_HD],args->ld_annot_pos[VCFBUF_LD_IDX_HD]);
            bcf_hdr_printf(args->hdr,"##INFO=<ID=%s,Number=1,Type=Integer,Description=\"The position of the site for which %s was calculated\">",args->ld_annot_pos[VCFBUF_LD_IDX_HD],args->ld_annot[VCFBUF_LD_IDX_HD]);
        }
    }
    if ( bcf_hdr_write(args->out_fh, args->hdr)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
    if ( init_index2(args->out_fh,args->hdr,args->output_fname,
                     &args->index_fn, args->write_index)<0 )
        error("Error: failed to initialise index for %s\n",args->output_fname);
    args->ld_filter_id = -1;
    if ( args->ld_filter && strcmp(".",args->ld_filter) )
        args->ld_filter_id = bcf_hdr_id2int(args->hdr, BCF_DT_ID, args->ld_filter);

    args->vcfbuf = vcfbuf_init(args->hdr, args->ld_win);
    if ( args->ld_max_set[VCFBUF_LD_IDX_R2] ) vcfbuf_set(args->vcfbuf,LD_MAX_R2,args->ld_max[VCFBUF_LD_IDX_R2]);
    if ( args->ld_max_set[VCFBUF_LD_IDX_LD] ) vcfbuf_set(args->vcfbuf,LD_MAX_LD,args->ld_max[VCFBUF_LD_IDX_LD]);
    if ( args->ld_max_set[VCFBUF_LD_IDX_HD] ) vcfbuf_set(args->vcfbuf,LD_MAX_HD,args->ld_max[VCFBUF_LD_IDX_HD]);
    if ( args->rand_missing || (args->nsites_mode && !strcasecmp(args->nsites_mode,"rand")) )
    {
        fprintf(stderr,"Using random seed: %d\n",args->rseed);
        hts_srand48(args->rseed);
    }
    if ( args->rand_missing ) vcfbuf_set(args->vcfbuf,LD_RAND_MISSING,1);
    if ( args->nsites )
    {
        vcfbuf_set(args->vcfbuf,PRUNE_NSITES,args->nsites);
        vcfbuf_set(args->vcfbuf,PRUNE_NSITES_MODE,args->nsites_mode);
    }
    if ( args->af_tag ) vcfbuf_set(args->vcfbuf,PRUNE_AF_TAG,args->af_tag);

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
    if ( hts_close(args->out_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname);
    vcfbuf_destroy(args->vcfbuf);
    bcf_sr_destroy(args->sr);
    free(args);
}
static void flush(args_t *args, int flush_all)
{
    bcf1_t *rec;
    while ( (rec = vcfbuf_flush(args->vcfbuf, flush_all)) )
        if ( bcf_write1(args->out_fh, args->hdr, rec)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
}
static void process(args_t *args)
{
    int i, filter = 0;
    bcf1_t *rec = bcf_sr_get_line(args->sr,0);
    if ( args->filter )
    {
        int ret  = filter_test(args->filter, rec, NULL);
        if ( args->filter_logic==FLT_INCLUDE ) { if ( !ret ) filter = 1; }
        else if ( ret ) filter = 1;
        if ( filter && !args->keep_sites ) return;
    }
    bcf_sr_t *sr = bcf_sr_get_reader(args->sr, 0);
    if ( args->ld_mask )
    {
        vcfbuf_ld_t ld;
        if ( vcfbuf_ld(args->vcfbuf, rec, &ld) == 0 )
        {
            int pass = 1;
            for (i=0; i<VCFBUF_LD_N; i++)
            {
                if ( !args->ld_max_set[i] || ld.val[i] <= args->ld_max[i] ) continue;
                pass = 0;
                break;
            }
            if ( !pass )
            {
                if ( !args->ld_filter ) return;     // hard filter
                if ( args->ld_filter_id >= 0 )      // soft FILTER
                    bcf_add_filter(args->hdr, rec, args->ld_filter_id);
            }
            for (i=0; i<VCFBUF_LD_N; i++)
            {
                if ( args->ld_annot[i] )
                {
                    int32_t pos = ld.rec[i]->pos+1;
                    bcf_update_info_int32(args->hdr, rec, args->ld_annot_pos[i], &pos, 1);
                }
            }
            for (i=0; i<VCFBUF_LD_N; i++)
            {
                if ( args->ld_annot[i] )
                {
                    float val   = ld.val[i];
                    bcf_update_info_float(args->hdr, rec, args->ld_annot[i], &val, 1);
                }
            }
        }
    }
    if ( filter ) vcfbuf_set(args->vcfbuf,LD_FILTER1,1);
    sr->buffer[0] = vcfbuf_push(args->vcfbuf, rec);
    flush(args,0);
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_type  = FT_VCF;
    args->output_fname = "-";
    args->ld_win = -100e3;
    args->nsites_mode = "maxAF";
    args->rseed = time(NULL);
    args->clevel = -1;
    static struct option loptions[] =
    {
        {"keep-sites",no_argument,NULL,'k'},
        {"randomize-missing",no_argument,NULL,1},
        {"AF-tag",required_argument,NULL,2},
        {"random-seed",required_argument,NULL,3},
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"annotate",required_argument,NULL,'a'},
        {"set-filter",required_argument,NULL,'f'},
        {"max",required_argument,NULL,'m'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"nsites-per-win",required_argument,NULL,'n'},
        {"nsites-per-win-mode",required_argument,NULL,'N'},
        {"window",required_argument,NULL,'w'},
        {"write-index",optional_argument,NULL,'W'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "vr:R:t:T:m:o:O:a:f:i:e:n:N:w:kW::",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case  1 : args->rand_missing = 1; break;
            case  2 : args->af_tag = optarg; break;
            case  3 :
                args->rseed = strtol(optarg,&tmp,10);
                if ( tmp==optarg || *tmp ) error("Could not parse: --random-seed %s\n", optarg);
                break;
            case 'W':
                if (!(args->write_index = write_index_parse(optarg)))
                    error("Unsupported index format '%s'\n", optarg);
                break;
            case 'k': args->keep_sites = 1; break;
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'a':
                {
                    int n, i;
                    char **tag = hts_readlist(optarg,0,&n);
                    if ( !n || !tag ) error("Could not parse: --annotate %s\n",optarg);
                    for (i=0; i<n; i++)
                    {
                        if ( !strcasecmp("R2",tag[i]) )
                        {
                            args->ld_annot_pos[VCFBUF_LD_IDX_R2] = "POS_R2";
                            args->ld_annot[VCFBUF_LD_IDX_R2]     = "R2";
                        }
                        else if ( !strcasecmp("LD",tag[i]) )
                        {
                            args->ld_annot_pos[VCFBUF_LD_IDX_LD] = "POS_LD";
                            args->ld_annot[VCFBUF_LD_IDX_LD]     = "LD";
                        }
                        else if ( !strcasecmp("HD",tag[i]) )
                        {
                            args->ld_annot_pos[VCFBUF_LD_IDX_HD] = "POS_HD";
                            args->ld_annot[VCFBUF_LD_IDX_HD]     = "HD";
                        }
                        else error("The tag \"%s\" is not supported\n",tag[i]);
                        free(tag[i]);
                    }
                    free(tag);
                    args->ld_mask |= LD_ANNOTATE;
                }
                break;
            case 'f': args->ld_filter = optarg; break;
            case 'n':
                args->nsites = strtod(optarg,&tmp);
                if ( tmp==optarg || *tmp ) error("Could not parse: --nsites-per-win %s\n", optarg);
                break;
            case 'N':
                if ( !strcasecmp(optarg,"maxAF") ) args->nsites_mode = optarg;
                else if ( !strcasecmp(optarg,"1st") ) args->nsites_mode = optarg;
                else if ( !strcasecmp(optarg,"rand") ) args->nsites_mode = optarg;
                else error("The mode \"%s\" is not recognised\n",optarg);
                break;
            case 'm':
                if ( !strncasecmp("R2=",optarg,3) )
                {
                    args->ld_max_set[VCFBUF_LD_IDX_R2] = 1;
                    args->ld_max[VCFBUF_LD_IDX_R2] = strtod(optarg+3,&tmp);
                }
                else if ( !strncasecmp("LD=",optarg,3) )
                {
                    args->ld_max_set[VCFBUF_LD_IDX_LD] = 1;
                    args->ld_max[VCFBUF_LD_IDX_LD] = strtod(optarg+3,&tmp);
                }
                else if ( !strncasecmp("HD=",optarg,3) )
                {
                    args->ld_max_set[VCFBUF_LD_IDX_HD] = 1;
                    args->ld_max[VCFBUF_LD_IDX_HD] = strtod(optarg+3,&tmp);
                }
                else
                {
                    args->ld_max_set[VCFBUF_LD_IDX_R2] = 1;
                    args->ld_max[VCFBUF_LD_IDX_R2] = strtod(optarg,&tmp);
                }
                if ( !tmp || *tmp ) error("Could not parse: --max %s\n", optarg);
                args->ld_mask |= LD_SET_MAX;
                break;
            case 'w':
                args->ld_win = strtod(optarg,&tmp);
                if ( !*tmp ) break;
                if ( tmp==optarg ) error("Could not parse: --window %s\n", optarg);
                else if ( !strcasecmp("bp",tmp) ) args->ld_win *= -1;
                else if ( !strcasecmp("kb",tmp) ) args->ld_win *= -1000;
                else if ( !strcasecmp("Mb",tmp) ) args->ld_win *= -1000000;
                else error("Could not parse: --window %s\n", optarg);
                break;
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
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");
    if ( !args->ld_mask && !args->nsites ) error("%sError: Expected pruning (--max,--nsites-per-win) or annotation (--annotate) options\n\n", usage_text());
    if ( args->ld_filter && strcmp(".",args->ld_filter) && !(args->ld_mask & LD_SET_MAX) ) error("The --set-filter option requires --max.\n");
    if ( args->keep_sites && args->nsites ) error("The --keep-sites option cannot be combined with --nsites-per-win\n");

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

    destroy_data(args);
    return 0;
}


