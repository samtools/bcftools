/* The MIT License

    Copyright (C) 2020-2023 Giulio Genovese

    Author: Giulio Genovese <giulio.genovese@gmail.com>

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

#include <getopt.h>
#include <unistd.h>
#include <ctype.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/hts_defs.h>
#include "bcftools.h"
#include "htslib/khash_str2int.h"
#include "regidx.h"

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct
{
    htsFile *fh;        // output file handle
    char *fname;        // output file name
    char *index_fn;
}
subset_t;

typedef struct
{
    char *filter_str;
    int filter_logic;   // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)
    void *hash;
    regidx_t *reg_idx;
    regitr_t *reg_itr;
    int argc, region_is_file, target_is_file, nsites, chunk_cnt, rec_cnt, scatter_is_file, output_type, n_threads, clevel;
    int regions_overlap, targets_overlap;
    char **argv, *region, *target, *scatter, *fname, *prefix, *extra, *output_dir;
    bcf_srs_t *sr;
    kstring_t str;
    subset_t *sets;
    int nsets, msets;
    int record_cmd_line;
    char **hts_opts;
    int nhts_opts;
    bcf_hdr_t *hdr;
    int write_index;
}
args_t;

const char *about(void)
{
    return "Scatter VCF by chunks or regions, creating multiple VCFs.\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: Split VCF by chunks or regions, creating multiple VCFs.\n"
        "\n"
        "Usage: bcftools +scatter [Options]\n"
        "Plugin options:\n"
        "   -e, --exclude EXPR              Exclude sites for which the expression is true\n"
        "   -i, --include EXPR              Select sites for which the expression is true\n"
        "       --no-version                Do not append version and command line to the header\n"
        "   -o, --output DIR                Write output to the directory DIR\n"
        "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n"
        "       --threads INT               Use multithreading with INT worker threads [0]\n"
        "   -r, --regions REGION            Restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE         Restrict to regions listed in a file\n"
        "       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n"
        "   -t, --targets [^]REGION         Similar to -r but streams rather than index-jumps. Exclude regions with \"^\" prefix\n"
        "   -T, --targets-file [^]FILE      Similar to -R but streams rather than index-jumps. Exclude regions with \"^\" prefix\n"
        "       --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n"
        "   -n, --nsites-per-chunk INT      Keep N sites for each chunk\n"
        "   -s, --scatter REGION            Comma-separated list of regions defining variant windows for each output VCF\n"
        "   -S, --scatter-file FILE         Regions listed in a file with an optional second column used to name each output VCF.\n"
        "                                        Variants from multiple regions can be assigned to the same output VCF.\n"
        "   -x, --extra STRING              Output records not overlapping listed regions in separate file\n"
        "   -p, --prefix STRING             Prepend string to output VCF names\n"
        "       --hts-opts LIST             Low-level options to pass to HTSlib, e.g. block_size=32768\n"
        "   -W, --write-index[=FMT]         Automatically index the output files [off]\n"
        "\n"
        "Examples:\n"
        "   # Scatter a VCF file by shards with 10000 variants each\n"
        "   bcftools +scatter input.bcf -Ob -o dir -n 10000 -p part-\n"
        "\n"
        "   # Scatter a VCF file by chromosomes\n"
        "   bcftools +scatter input.bcf -Ob -o dir -s chr1,chr2,chr3,<...>,chr22,chrX -x other\n"
        "\n";
}

void mkdir_p(const char *fmt, ...) HTS_FORMAT(HTS_PRINTF_FMT, 1, 2);

// most of this code was inspired by Petr Danecek's code in regidx.c
#define MAX_COOR_0 REGIDX_MAX   // CSI and hts_itr_query limit, 0-based
int regidx_parse_reg_name(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload, void *usr)
{
    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments

    char *se = ss;
    while ( *se && *se!=':' && !isspace(*se) ) se++;

    *chr_beg = ss;
    *chr_end = se-1;

    if ( *se != ':' )
    {
        *beg = 0;
        *end = MAX_COOR_0;
    }
    else
    {
        ss = se+1;
        *beg = strtod(ss, &se);
        if ( ss==se ) { fprintf(stderr,"Could not parse reg line: %s\n", line); return -2; }
        if ( *beg==0 ) { fprintf(stderr,"Could not parse reg line, expected 1-based coordinate: %s\n", line); return -2; }
        (*beg)--;

        if ( !se[0] || isspace(se[0])) {
            *end = *beg;
        } else if ( se[0] == '-' && (!se[1] || isspace(se[1])) ) {
           *end = MAX_COOR_0;
           se++;
        } else {
            ss = se+1;
            *end = strtod(ss, &se);
            if ( ss==se ) *end = *beg;
            else if ( *end==0 ) { fprintf(stderr,"Could not parse reg line, expected 1-based coordinate: %s\n", line); return -2; }
            else (*end)--;
        }
    }

    ss = se;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !ss[0] ) ss = (char *)line;

    int *idx = (int *)payload;
    args_t *args = (args_t *)usr;
    int ret = khash_str2int_get(args->hash, ss, idx);
    if (ret < 0) {
        hts_expand0(subset_t, args->nsets + 1, args->msets, args->sets);
        args->sets[args->nsets].fname = strdup(ss);
        *idx = khash_str2int_inc(args->hash, args->sets[args->nsets].fname);
        args->nsets++;
    }

    return 0;
}

static void open_set(subset_t *set, args_t *args)
{
    int j;
    args->str.l = 0;
    kputs(args->output_dir, &args->str);
    if ( args->str.s[args->str.l-1] != '/' ) kputc('/', &args->str);
    int k, l = args->str.l;
    if (args->prefix) kputs(args->prefix, &args->str);
    kputs(set->fname, &args->str);
    for (k=l; k<args->str.l; k++) if ( isspace(args->str.s[k]) ) args->str.s[k] = '_';
    if ( args->output_type & FT_BCF ) kputs(".bcf", &args->str);
    else if ( args->output_type & FT_GZ ) kputs(".vcf.gz", &args->str);
    else kputs(".vcf", &args->str);

    char wmode[8];
    set_wmode(wmode,args->output_type,args->str.s,args->clevel);
    set->fh = hts_open(args->str.s, wmode);
    if ( set->fh == NULL ) error("[%s] Error: cannot write to \"%s\": %s\n", __func__, args->str.s, strerror(errno));
    if ( args->n_threads > 0)
        hts_set_opt(set->fh, HTS_OPT_THREAD_POOL, args->sr->p);
    if ( args->hts_opts )
    {
        hts_opt *opts = NULL;
        for (j=0; j<args->nhts_opts; j++)
            if ( hts_opt_add(&opts, args->hts_opts[j]) ) error("Could not set the HTS option \"%s\"\n", args->hts_opts[j]);
        if ( hts_opt_apply(set->fh, opts) ) error("Could not apply the HTS options\n");
        hts_opt_free(opts);
    }
    if ( !args->hdr )
    {
        args->hdr = bcf_sr_get_header(args->sr, 0);
        if ( args->record_cmd_line ) bcf_hdr_append_version(args->hdr, args->argc, args->argv, "bcftools_plugin");
    }
    if ( bcf_hdr_write(set->fh, args->hdr)!=0 ) error("[%s] Error: cannot write the header to %s\n", __func__, args->str.s);
    if ( init_index2(set->fh,args->hdr,args->str.s,&set->index_fn,
                     args->write_index)<0 )
        error("Error: failed to initialise index for %s\n",args->str.s);
}

static void init_data(args_t *args)
{
    args->sr = bcf_sr_init();
    if ( args->region )
    {
        args->sr->require_index = 1;
        bcf_sr_set_opt(args->sr, BCF_SR_REGIONS_OVERLAP, args->regions_overlap);
        if ( bcf_sr_set_regions(args->sr, args->region, args->region_is_file)<0 )
            error("Failed to read the regions: %s\n", args->region);
    }
    if ( args->target )
    {
        bcf_sr_set_opt(args->sr, BCF_SR_TARGETS_OVERLAP, args->targets_overlap);
        if ( bcf_sr_set_targets(args->sr, args->target, args->target_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->target);
    }
    if ( bcf_sr_set_threads(args->sr, args->n_threads)<0 ) error("Failed to create threads\n");
    if ( !bcf_sr_add_reader(args->sr, args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));

    mkdir_p("%s/", args->output_dir);

    if (args->nsites) {
        args->nsets = 1;
        hts_expand0(subset_t, args->nsets, args->msets, args->sets);
    } else {
        args->hash = khash_str2int_init();
        if (args->scatter_is_file) {
            args->reg_idx = regidx_init(args->scatter, regidx_parse_reg_name, NULL, sizeof(int), args);
            if ( !args->reg_idx ) exit(-1);
        } else {
            char *ss;
            for (ss = args->scatter; *ss; ss++) if (*ss == ',') *ss = '\n';
            args->reg_idx = regidx_init_string(args->scatter, regidx_parse_reg_name, NULL, sizeof(int), args);
            for (ss = args->scatter; *ss; ss++) if (*ss == '\n') *ss = ',';
        }
        args->reg_itr = regitr_init(args->reg_idx);
        if (args->extra) {
            hts_expand0(subset_t, args->nsets + 1, args->msets, args->sets);
            args->sets[args->nsets].fname = strdup(args->extra);
            args->nsets++;
        }
        int i;
        for (i=0; i<args->nsets; i++) open_set(&args->sets[i], args);
    }
}

static void destroy_data(args_t *args)
{
    if (args->scatter) {
        khash_str2int_destroy(args->hash);
        regidx_destroy(args->reg_idx);
        regitr_destroy(args->reg_itr);
    }
    free(args->str.s);
    int i;
    for (i=0; i<args->nsets; i++)
    {
        subset_t *set = &args->sets[i];
        if (set->fname)
        {
            if ( args->write_index )
            {
                if ( bcf_idx_save(set->fh)<0 )
                {
                    if ( hts_close(set->fh)!=0 ) error("Error: close failed .. %s\n", set->fname);
                    error("Error: cannot write to index %s\n", set->index_fn);
                }
                free(set->index_fn);
            }
            if ( hts_close(set->fh)!=0 ) error("Error: close failed .. %s\n", set->fname);
            free(set->fname);
        }
    }
    for (i=0; i<args->nhts_opts; i++) free(args->hts_opts[i]);
    free(args->hts_opts);
    free(args->sets);
    bcf_sr_destroy(args->sr);
    free(args);
}

static void process(args_t *args)
{
    bcf1_t *rec = bcf_sr_get_line(args->sr, 0);
    subset_t *set;

    if (args->nsites) {
        subset_t *set = &args->sets[0];
        if (!args->rec_cnt) {
            args->str.l = 0;
            kputw(args->chunk_cnt, &args->str);
            set->fname = strdup(args->str.s);
            open_set(set, args);
            args->rec_cnt = 0;
        }
        if ( bcf_write(set->fh, args->hdr, rec)!=0 ) error("[%s] Error: failed to write the record\n", __func__);
        args->rec_cnt++;
        if (args->rec_cnt == args->nsites) {
          args->rec_cnt = 0;
          int err = 0;
          if ( args->write_index && bcf_idx_save(set->fh) < 0)
              err = 1;
          if ( hts_close(set->fh)!=0 ) error("Error: close failed .. %s\n", set->fname);
          if (err)
              error("Error: cannot write to index %s\n", set->index_fn);
          free(set->fname);
          set->fname = NULL;
          args->chunk_cnt++;
        }
    } else {
        if ( regidx_overlap(args->reg_idx, bcf_hdr_id2name(args->hdr, rec->rid), rec->pos, rec->pos, args->reg_itr) ) {
            while (regitr_overlap(args->reg_itr)) {
                int idx = regitr_payload(args->reg_itr, int);
                set = &args->sets[idx];
                if ( bcf_write(set->fh, args->hdr, rec)!=0 ) error("[%s] Error: failed to write the record\n", __func__);
            }
        } else if (args->extra) {
            set = &args->sets[args->nsets-1];
            if ( bcf_write(set->fh, args->hdr, rec)!=0 ) error("[%s] Error: failed to write the record\n", __func__);
        }
    }
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1, sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_type  = FT_VCF;
    args->record_cmd_line = 1;
    args->regions_overlap = 1;
    args->targets_overlap = 0;
    args->clevel = -1;
    static struct option loptions[] =
    {
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"no-version",no_argument,NULL,1},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"threads",required_argument,NULL,2},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"regions-overlap",required_argument,NULL,3},
        {"targets", required_argument, NULL,'t'},
        {"targets-file", required_argument, NULL,'T'},
        {"targets-overlap",required_argument,NULL,4},
        {"nsites-per-chunk",required_argument,NULL,'n'},
        {"scatter",required_argument,NULL,'s'},
        {"scatter-file",required_argument,NULL,'S'},
        {"extra",required_argument,NULL,'x'},
        {"prefix",required_argument,NULL,'p'},
        {"hts-opts",required_argument,NULL,5},
        {"write-index",optional_argument,NULL,'W'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "e:i:o:O:r:R:t:T:n:s:S:x:p:W::h?", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case  1 : args->record_cmd_line = 0; break;
            case 'o': args->output_dir = optarg; break;
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
            case  2 : args->n_threads = strtol(optarg, 0, 0); break;
            case 'r': args->region = optarg; break;
            case 'R': args->region = optarg; args->region_is_file = 1;  break;
            case  3 :
                args->regions_overlap = parse_overlap_option(optarg);
                if ( args->regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case 't': args->target = optarg; break;
            case 'T': args->target = optarg; args->target_is_file = 1; break;
            case  4 :
                args->targets_overlap = parse_overlap_option(optarg);
                if ( args->targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case 'n':
                args->nsites = strtod(optarg, &tmp);
                if ( tmp==optarg || *tmp ) error("Could not parse: --nsites-per-chunk %s\n", optarg);
                if ( args->nsites <= 0 ) error("Positive integer required: --nsites-per-chunk %s\n", optarg);
                break;
            case 's': args->scatter = optarg; break;
            case 'S': args->scatter = optarg; args->scatter_is_file = 1;  break;
            case 'x': args->extra = optarg;  break;
            case 'p': args->prefix = optarg;  break;
            case  5 : args->hts_opts = hts_readlist(optarg, 0, &args->nhts_opts); break;
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

    if ( !args->output_dir ) error("Missing the -o option\n");
    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given\n");
    if ( !args->nsites && !args->scatter ) error("Missing either the -n or one of the -s or -S options\n");
    if ( args->nsites && args->scatter ) error("Only one of -n or either -s or -S can be given\n");
    if ( args->nsites && args->extra ) error("Cannot use -x together with -n\n");

    init_data(args);

    while ( bcf_sr_next_line(args->sr) ) process(args);

    destroy_data(args);
    return 0;
}
