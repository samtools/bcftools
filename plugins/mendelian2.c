/* The MIT License

   Copyright (c) 2015-2024 Genome Research Ltd.

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
#include <string.h>
#include <strings.h>
#include <getopt.h>
#include <math.h>
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/khash_str2int.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>     // for isatty
#include "bcftools.h"
#include "regidx.h"
#include "filter.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define MODE_ANNOTATE  (1<<0)   // -m a
#define MODE_COUNT     (1<<1)   // -m c
#define MODE_DELETE    (1<<2)   // -m d
#define MODE_LIST_ERR  (1<<3)   // -m e
#define MODE_DROP_ERR  (1<<4)   // -m E
#define MODE_LIST_GOOD (1<<5)   // -m g
#define MODE_LIST_MISS (1<<6)   // -m m
#define MODE_DROP_MISS (1<<7)   // -m M
#define MODE_LIST_SKIP (1<<8)   // -m s
#define MODE_DROP_SKIP (1<<9)   // -m S
#define LIST_MODES (MODE_LIST_ERR|MODE_LIST_GOOD|MODE_LIST_MISS|MODE_LIST_SKIP)

#define iMOM 0
#define iDAD 1
#define iKID 2

typedef struct
{
    int nno_gt,         // number of rows with no FMT/GT
        nnot_diploid,   // FMT/GT2 not diploid
        nfail,          // number of -i/-e filters failed
        nmiss,          // number of genotypes with a missing allele in the trio
        ngood,          // number of good genotypes (after any -i/-e filters applied)
        nmerr,          // number of mendelian errors
        ngood_alt,      // number of error-free non-ref genotypes
        nrule;          // number of genotypes with no rule to apply
}
stats_t;

typedef struct
{
    int idx[3];         // VCF sample index of mom,dad,kid, in that order
    stats_t stats;      // per trio stats collected over all sites
    int fail, has_merr; // per site filtering and mendeian error info
    int sex_id;         // implicitly defined by rules (typically 1X=0, 2X=1), see str2sex_id
}
trio_t;

typedef struct
{
   uint8_t
        ploidy:4,       // 0,1,2
        sex_id:2,       // id given by str2sex_id mapping
        inherits:2;     // one of the M,F,MF,. inheritance modes given as combination of 1<<iDAD,1<<iMOM bitmasks
}
rule_t;

typedef struct _args_t
{
    int argc, filter_logic, regions_is_file, targets_is_file, output_type, record_cmd_line, clevel;
    int regions_overlap, targets_overlap;
    char *filter_str;
    filter_t *filter;
    char **argv, *ped_fname, *pfm, *output_fname, *fname, *regions, *targets, *rules_str, *rules_fname;
    htsFile *out_fh;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr, *hdr_out;
    trio_t *trio;
    int ntrio;
    regitr_t *itr;
    regidx_t *rules;            // inheritance rules per region, unlisted assumed to have "MF" inheritance for both sexes
    void *str2sex_id;           // implicitly defined by rules, the default is 0="1X" (male) and 1="2X" (female)
    int nsex_id;                // number of sexes (expect 2 for human
    rule_t *rule;               // the current site's inheritance for both sexes
    int mode;
    int32_t *gt_arr, nmerr;
    int ngt_arr;
    stats_t stats;              // common per-site and per-sample stats
    int nref_only, nmany_als;   // per-site stats
    char *index_fn;
    int write_index;
}
args_t;

const char *about(void)
{
    return "Count Mendelian consistent / inconsistent genotypes.\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: Count Mendelian consistent / inconsistent genotypes.\n"
        "Usage: bcftools +mendelian2 [Options]\n"
        "Common Options:\n"
        "   -e, --exclude EXPR              Exclude sites for which the expression is true\n"
        "   -i, --include EXPR              Include sites for which the expression is true\n"
        "   -o, --output FILE               Output file name [stdout]\n"
        "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n"
        "   -r, --regions REG               Restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE         Restrict to regions listed in a file\n"
        "       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n"
        "   -t, --targets REG               Similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE         Similar to -R but streams rather than index-jumps\n"
        "       --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n"
        "       --no-version                Do not append version and command line to the header\n"
        "   -W, --write-index[=FMT]         Automatically index the output files [off]\n"
        "\n"
        "Options:\n"
        "   -m, --mode c|[adeEgmMS]         Output mode, the default is `-m c`. Multiple modes can be combined in VCF/BCF\n"
        "                                   output and the drop modes take precedence.\n"
        "                                   Text stats output:\n"
        "                                       c .. print counts, a text summary with the number of errors per trio\n"
        "                                   VCF/BCF output:\n"
        "                                       a .. add INFO/MERR annotation with the number of inconsistent trios\n"
        "                                       d .. delete genotypes in inconsistent trios (set to \"./.\")\n"
        "                                       e .. output sites with at least one erroneous trio\n"
        "                                       E .. drop sites with at least one erroneous trio\n"
        "                                       g .. output sites with at least one good trio (i.e. non-missing and consistent)\n"
        "                                       m .. output sites with missing genotypes in at least one trio\n"
        "                                       M .. drop sites with missing genotypes in at least one trio\n"
        "                                       S .. drop sites skipped for various reasons when collecting stats\n"
        "   -p, --pfm [1X:|2X:]P,F,M        Sample names of child (the proband), father, mother; \"1X:\" for male pattern of chrX inheritance [2X:]\n"
        "   -P, --ped FILE                  PED file with the columns: <ignored>,proband,father,mother,sex(1:male,2:female)\n"
        "       --rules ASSEMBLY[?]         Predefined inheritance rules, \"list\" to print available settings, \"list?\" for detailed information\n"
        "       --rules-file FILE           Inheritance rules, run with `--rules list?` for examples\n"
        "\n"
        "Example:\n"
        "   # Print number of good, erroneous and missing genotypes\n"
        "   bcftools +mendelian2 in.vcf -p 1X:Child,Father,Mother -c\n"
        "\n";
}

typedef struct
{
    const char *alias, *about, *rules;
}
rules_predef_t;

static rules_predef_t rules_predefs[] =
{
    { .alias = "GRCh37",
      .about = "Human Genome reference assembly GRCh37 / hg19, both chr naming conventions",
      .rules =
            "   # Unlisted regions inherit from both parents (MF)\n"
            "   1X  X:1-60000               M\n"
            "   1X  X:2699521-154931043     M\n"
            "   1X  Y:1-59373566            F\n"
            "   2X  Y:1-59373566            .\n"
            "   1X  MT:1-16569              M\n"
            "   2X  MT:1-16569              M\n"
            "\n"
            "   1X  chrX:1-60000            M\n"
            "   1X  chrX:2699521-154931043  M\n"
            "   1X  chrY:1-59373566         F\n"
            "   2X  chrY:1-59373566         .\n"
            "   1X  chrM:1-16569            M\n"
            "   2X  chrM:1-16569            M\n"
    },
    { .alias = "GRCh38",
      .about = "Human Genome reference assembly GRCh38 / hg38, both chr naming conventions",
      .rules =
            "   # Unlisted regions inherit from both parents (MF)\n"
            "   1X  X:1-9999                M\n"
            "   1X  X:2781480-155701381     M\n"
            "   1X  Y:1-57227415            F\n"
            "   2X  Y:1-57227415            .\n"
            "   1X  MT:1-16569              M\n"
            "   2X  MT:1-16569              M\n"
            "\n"
            "   1X  chrX:1-9999             M\n"
            "   1X  chrX:2781480-155701381  M\n"
            "   1X  chrY:1-57227415         F\n"
            "   2X  chrY:1-57227415         .\n"
            "   1X  chrMT:1-16569           M\n"
            "   2X  chrMT:1-16569           M\n"
    },
    {
        .alias = NULL,
        .about = NULL,
        .rules = NULL,
    }
};

static int parse_rules(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload, void *usr)
{
    args_t *args = (args_t*)usr;

    // e.g
    //  1X  Y:1-59373566  F     # son
    //  2X  Y:1-59373566  .     # daughter

    if ( line[0]=='#' ) return -1;  // skip comment lines

    // eat any leading spaces
    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) return -1;      // skip empty lines

    // sex id, e.g. 1X or 2X
    char keep, *tmp, *se = ss;
    while ( *se && !isspace(*se) ) se++;
    if ( !*se ) error("Could not parse the sex ID in the region line: %s\n", line);
    keep = *se;
    *se = 0;
    int sex_id = 0;
    if ( khash_str2int_get(args->str2sex_id,ss,&sex_id)!=0 )
    {
        tmp = strdup(ss);
        if ( khash_str2int_set(args->str2sex_id,tmp,args->nsex_id) < 0 ) error("Could not insert into a hash\n");
        sex_id = args->nsex_id++;
    }
    *se = keep;
    while ( *se && isdigit(*se) ) se++;
    while ( *se && isspace(*se) ) se++;
    ss = se;

    // chromosome name, beg, end
    while ( se[1] && !isspace(se[1]) ) se++;
    while ( se > ss && isdigit(*se) ) se--;
    if ( *se!='-' ) error("Could not parse the region: %s\n",line);
    *end = strtol(se+1, &tmp, 10) - 1;
    if ( tmp==se+1 ) error("Could not parse the region: %s\n",line);
    while ( se > ss && *se!=':' ) se--;
    *beg = strtol(se+1, &tmp, 10) - 1;
    if ( tmp==se+1 ) error("Could not parse the region: %s\n",line);

    *chr_beg = ss;
    *chr_end = se-1;

    // skip region
    while ( *ss && !isspace(*ss) ) ss++;
    while ( *ss && isspace(*ss) ) ss++;

    rule_t *rule = (rule_t*) payload;
    rule->sex_id = sex_id;
    rule->inherits = 0;
    rule->ploidy   = 0;

    // alleles inherited from mother (M), father (F), both (MF), none (.)
    while ( *ss && !isspace(*ss) )
    {
        if ( *ss=='M' ) { rule->inherits |= 1<<iMOM; rule->ploidy++; }
        else if ( *ss=='F' ) { rule->inherits |= 1<<iDAD; rule->ploidy++; }
        else if ( *ss=='.' ) { rule->inherits = 0; rule->ploidy = 0; }
        else error("Could not parse the region: %s\n", line);
        ss++;
        if ( rule->ploidy > 2 ) error("Todo: not ready for ploidy > 2: %s\n", line);
    }
    return 0;
}
regidx_t *init_rules(args_t *args, char *alias)
{
    const rules_predef_t *rules = rules_predefs;
    if ( !alias ) alias = "GRCh37";

    int detailed = 0, len = strlen(alias);
    if ( alias[len-1]=='?' ) { detailed = 1; alias[len-1] = 0; }

    while ( rules->alias && strcasecmp(alias,rules->alias) ) rules++;

    if ( !rules->alias )
    {
        fprintf(stderr,"\nPRE-DEFINED INHERITANCE RULES\n\n");
        fprintf(stderr," * Columns are: SEX_ID CHROM:BEG-END INHERITED_FROM\n");
        fprintf(stderr," * Coordinates are 1-based inclusive.\n\n");
        rules = rules_predefs;
        while ( rules->alias )
        {
            fprintf(stderr,"%s\n   .. %s\n\n", rules->alias,rules->about);
            if ( detailed )
                fprintf(stderr,"%s\n", rules->rules);
            rules++;
        }
        fprintf(stderr,"Run as --rules <alias> (e.g. --rules GRCh37).\n");
        fprintf(stderr,"To see the detailed ploidy definition, append a question mark (e.g. --rules GRCh37?).\n");
        fprintf(stderr,"\n");
        exit(-1);
    }
    else if ( detailed )
    {
        fprintf(stderr,"%s", rules->rules);
        exit(-1);
    }
    return regidx_init_string(rules->rules, parse_rules, NULL, sizeof(rule_t), args);
}

static int cmp_trios(const void *_a, const void *_b)
{
    trio_t *a = (trio_t *) _a;
    trio_t *b = (trio_t *) _b;
    int i;
    int amin = a->idx[0];
    for (i=1; i<3; i++)
        if ( amin > a->idx[i] ) amin = a->idx[i];
    int bmin = b->idx[0];
    for (i=1; i<3; i++)
        if ( bmin > b->idx[i] ) bmin = b->idx[i];
    if ( amin < bmin ) return -1;
    if ( amin > bmin ) return 1;
    return 0;
}
static void parse_ped(args_t *args, char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);

    kstring_t str = {0,0,0};
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Empty file: %s\n", fname);

    int moff = 0, *off = NULL, mtrio = 0;
    do
    {
        // familyID    sampleID paternalID maternalID sex   phenotype   population relationship   siblings   secondOrder   thirdOrder   children    comment
        // BB03    HG01884 HG01885 HG01956 2   0   ACB child   0   0   0   0
        int ncols = ksplit_core(str.s,0,&moff,&off);
        if ( ncols<4 ) error("Could not parse the ped file: %s\n", str.s);

        int dad = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[2]]);
        if ( dad < 0 ) continue;
        int mom = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[3]]);
        if ( mom < 0 ) continue;
        int kid = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[1]]);
        if ( kid < 0 ) continue;

        int sex = 0;
        if ( ncols>=5 )
        {
            char *tmp;
            sex = strtol(&str.s[off[4]],&tmp,10);
            if ( tmp==&str.s[off[4]] || *tmp ) error("Could not parse the PED file, the 5th column should be numeric: %s\n",str.s);
            if ( sex!=1 && sex!=2 ) sex = 0;
        }

        args->ntrio++;
        hts_expand0(trio_t,args->ntrio,mtrio,args->trio);
        trio_t *trio = &args->trio[args->ntrio-1];
        trio->idx[iDAD] = dad;
        trio->idx[iMOM] = mom;
        trio->idx[iKID] = kid;
        if ( khash_str2int_get(args->str2sex_id, sex==1 ? "1X" : "2X", &trio->sex_id)<0 )
            error("Missing the sex \"%s\", it's not in the rules :-/\n",sex==1 ? "1X" : "2X");
    }
    while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );

    // sort the sample by index so that they are accessed more or less sequentially
    qsort(args->trio,args->ntrio,sizeof(trio_t),cmp_trios);

    // check for duplicates
    int i;
    for (i=1; i<args->ntrio; i++)
    {
        trio_t *ta = &args->trio[i-1];
        trio_t *tb = &args->trio[i];
        if ( ta->idx[0]==tb->idx[0] && ta->idx[1]==tb->idx[1] && ta->idx[2]==tb->idx[2] )
            error("Error: duplicate trio entries detected in the PED file: %s\n",fname);
    }

    fprintf(stderr,"Identified %d complete trio%s in the VCF file\n", args->ntrio,args->ntrio==1?"":"s");

    free(str.s);
    free(off);
    if ( hts_close(fp)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,fname);
}

static void init_data(args_t *args)
{
    args->sr = bcf_sr_init();
    if ( args->regions )
    {
        args->sr->require_index = 1;
        bcf_sr_set_opt(args->sr,BCF_SR_REGIONS_OVERLAP,args->regions_overlap);
        if ( bcf_sr_set_regions(args->sr, args->regions, args->regions_is_file)<0 ) error("Failed to read the regions: %s\n",args->regions);
    }
    if ( args->targets )
    {
        bcf_sr_set_opt(args->sr,BCF_SR_TARGETS_OVERLAP,args->targets_overlap);
        if ( bcf_sr_set_targets(args->sr, args->targets, args->targets_is_file, 0)<0 ) error("Failed to read the targets: %s\n",args->targets);
    }
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr = bcf_sr_get_header(args->sr,0);

    if ( args->filter_str ) args->filter = filter_init(args->hdr, args->filter_str);

    args->str2sex_id = khash_str2int_init();
    if ( args->rules_fname )
        args->rules = regidx_init(args->rules_fname, parse_rules, NULL, sizeof(rule_t), args);
    else
        args->rules = init_rules(args, args->rules_str);
    if ( !args->rules ) error("Could not parse the Mendelian rules\n");
    args->itr  = regitr_init(args->rules);
    args->rule = (rule_t*) malloc(sizeof(*args->rule)*args->nsex_id);

    int i, n = 0;
    char **list;
    if ( args->pfm )
    {
        args->ntrio = 1;
        args->trio  = (trio_t*) calloc(1,sizeof(trio_t));
        list = hts_readlist(args->pfm, 0, &n);
        if ( n!=3 ) error("Expected three sample names with -t\n");
        const int ped_idx[3] = {2,1,0};  // sample order is different on the command line (P,F,M) and in the code (M,F,P)
        for (i=0; i<3; i++)
            args->trio[0].idx[i] = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, list[ped_idx[i]]);
        if ( args->trio[0].idx[iKID] < 0 )
        {
            if ( strlen(list[0])>3 && !strncasecmp(list[0],"1X:",3) )
            {
                args->trio[0].idx[iKID] = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, list[0]+3);
                args->trio[0].sex_id = iDAD;
            }
            else if ( strlen(list[0])>3 && !strncasecmp(list[0],"2X:",3) )
            {
                args->trio[0].idx[iKID] = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, list[0]+3);
                args->trio[0].sex_id = iMOM;
            }
        }
        for (i=0; i<3; i++)
        {
            if ( args->trio[0].idx[i] < 0 ) error("The sample is not present: %s\n", list[ped_idx[i]]);
            free(list[ped_idx[i]]);
        }
        free(list);
    }
    else
    {
        parse_ped(args,args->ped_fname);
        if ( !args->ntrio ) error("No complete trio present\n");
    }

    args->hdr_out = bcf_hdr_dup(args->hdr);
    if ( args->mode&MODE_ANNOTATE )
        bcf_hdr_append(args->hdr_out,"##INFO=<ID=MERR,Number=1,Type=Integer,Description=\"Number of trios with Mendelian errors\">");
    if ( args->record_cmd_line )
        bcf_hdr_append_version(args->hdr_out, args->argc, args->argv, "bcftools_trio-dnm2");

    if ( args->mode!=MODE_COUNT )
    {
        char wmode[8];
        set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
        args->out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
        if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
        if ( bcf_hdr_write(args->out_fh, args->hdr_out)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
        if ( init_index2(args->out_fh,args->hdr_out,args->output_fname,
                         &args->index_fn, args->write_index)<0 )
          error("Error: failed to initialise index for %s\n",args->output_fname);
    }
}

static void destroy_data(args_t *args)
{
    if ( args->filter ) filter_destroy(args->filter);
    regidx_destroy(args->rules);
    regitr_destroy(args->itr);
    khash_str2int_destroy_free(args->str2sex_id);
    free(args->trio);
    free(args->gt_arr);
    free(args->rule);
    if ( args->out_fh )
    {
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
    }
    bcf_hdr_destroy(args->hdr_out);
    bcf_sr_destroy(args->sr);
    free(args);
}

static int test_filters(args_t *args, bcf1_t *rec)
{
    uint8_t *smpl_pass;
    int i,j, pass_site = filter_test(args->filter, rec, (const uint8_t**) &smpl_pass);
    if ( args->filter_logic & FLT_EXCLUDE )
    {
        if ( pass_site )
        {
            if ( !smpl_pass ) return 0;     // no samples, -e mode, the expression failed
            pass_site = 0;
            for (i=0; i<args->ntrio; i++)
            {
                int pass_trio = 1;
                for (j=0; j<3; j++)
                {
                    int idx = args->trio[i].idx[j];
                    if ( smpl_pass[idx] ) { pass_trio = 0; break; }     // with -e one sample passes, the whole trio fails
                }
                args->trio[i].fail = pass_trio ? 0 : 1;
                if ( pass_trio ) pass_site = 1;
            }
            return pass_site;
        }
        for (i=0; i<args->ntrio; i++) args->trio[i].fail = 0;
        return 1;
    }
    if ( !pass_site ) return 0;
    if ( smpl_pass )
    {
        pass_site = 0;
        for (i=0; i<args->ntrio; i++)
        {
            int pass_trio = 1;
            for (j=0; j<3; j++)
            {
                int idx = args->trio[i].idx[j];
                if ( !smpl_pass[idx] ) { pass_trio = 0; break; }
            }
            args->trio[i].fail = pass_trio ? 0 : 1;
            if ( pass_trio ) pass_site = 1;
        }
        return pass_site;
    }
    for (i=0; i<args->ntrio; i++) args->trio[i].fail = 0;
    return 1;
}

static int parse_gt(int32_t *gt, int ngt, uint64_t *a, uint64_t *b)
{
    *a = *b = 0;

    if ( gt[0]==bcf_gt_missing || gt[0]==bcf_int32_vector_end ) return 0;
    *a |= 1<<bcf_gt_allele(gt[0]);

    if ( ngt==1 || gt[1]==bcf_int32_vector_end ) return 1;

    if ( gt[1]==bcf_gt_missing ) return 0;
    *b |= 1<<bcf_gt_allele(gt[1]);

    return 2;
}
static void delete_gt(int32_t *gt, int ngt)
{
    int i;
    for (i=0; i<ngt; i++)
    {
        if ( gt[i]==bcf_int32_vector_end ) return;
        gt[i] = bcf_gt_missing;
    }
}

#define HAS_GOOD 1
#define HAS_MERR 2
#define HAS_MISS 4
static int collect_stats(args_t *args, bcf1_t *rec)
{
    args->nmerr = 0;

    int ret = 0;
    int ngt = bcf_get_genotypes(args->hdr, rec, &args->gt_arr, &args->ngt_arr);
    if ( ngt<0 )
    {
        args->stats.nno_gt++;
        return ret;
    }
    if ( ngt!=2*bcf_hdr_nsamples(args->hdr) && ngt!=bcf_hdr_nsamples(args->hdr) )
    {
        args->stats.nnot_diploid++;
        return ret;
    }
    ngt /= bcf_hdr_nsamples(args->hdr);

    int i,j,itr_set = regidx_overlap(args->rules, bcf_seqname(args->hdr,rec),rec->pos,rec->pos+rec->rlen-1, args->itr);
    for (i=0; i<args->nsex_id; i++)
    {
        args->rule[i].sex_id   = i;
        args->rule[i].inherits = (1<<iMOM)|(1<<iDAD);
        args->rule[i].ploidy   = 2;
    }
    while ( itr_set && regitr_overlap(args->itr) )
    {
        rule_t *rule = &regitr_payload(args->itr,rule_t);
        args->rule[rule->sex_id] = *rule;
    }

    for (i=0; i<args->ntrio; i++)
    {
        trio_t *trio = &args->trio[i];
        trio->has_merr = 0;
        if ( trio->fail )
        {
            trio->stats.nfail++;
            continue;
        }
        rule_t *rule = &args->rule[trio->sex_id];
        if ( !rule->inherits ) { trio->stats.nrule++; continue; }
        uint64_t kid1, kid2, parent, mom, dad;
        int nal = parse_gt(&args->gt_arr[ngt*trio->idx[iKID]],ngt,&kid1,&kid2);
        if ( nal < rule->ploidy ) { ret |= HAS_MISS; trio->stats.nmiss++; continue; }
        if ( nal==1 )
        {
            for (j=0; j<args->nsex_id; j++)
            {
                if ( rule->inherits & (1<<j) ) break;
            }
            nal = parse_gt(&args->gt_arr[ngt*trio->idx[j]],ngt,&parent,&parent);
            if ( !nal ) { ret |= HAS_MISS; trio->stats.nmiss++; continue; }
            if ( parent&kid1 )
            {
                ret |= HAS_GOOD;
                trio->stats.ngood++;
                if ( parent!=1 || parent!=kid1 ) trio->stats.ngood_alt++;
                continue;
            }
            ret |= HAS_MERR;
            trio->stats.nmerr++;
            trio->has_merr = 1;
            args->nmerr++;
            continue;
        }
        int nal_mom = parse_gt(&args->gt_arr[ngt*trio->idx[iMOM]],ngt,&mom,&mom);
        int nal_dad = parse_gt(&args->gt_arr[ngt*trio->idx[iDAD]],ngt,&dad,&dad);
        if ( (kid1&dad && kid2&mom) || (kid1&mom && kid2&dad) )
        {
            // both children's alleles phased
            ret |= HAS_GOOD;
            trio->stats.ngood++;
            if ( dad!=1 || mom!=1 || (kid1|kid2)!=1 ) trio->stats.ngood_alt++;
            continue;
        }
        if ( !nal_mom || !nal_dad ) { ret |= HAS_MISS; trio->stats.nmiss++; }       // one or both parents missing
        if ( !nal_mom && !nal_dad ) continue;                                       // both parents missing
        if ( !nal_mom && ((kid1|kid2)&dad) ) continue;                              // one parent missing but the kid is consistent with the other
        if ( !nal_dad && ((kid1|kid2)&mom) ) continue;
        ret |= HAS_MERR;
        trio->stats.nmerr++;
        trio->has_merr = 1;
        args->nmerr++;
    }

    if ( !(args->mode&MODE_DELETE) || !(ret&HAS_MERR) ) return ret;
    for (i=0; i<args->ntrio; i++)
    {
        trio_t *trio = &args->trio[i];
        if ( !trio->has_merr ) continue;
        delete_gt(&args->gt_arr[ngt*trio->idx[iKID]],ngt);
        delete_gt(&args->gt_arr[ngt*trio->idx[iDAD]],ngt);
        delete_gt(&args->gt_arr[ngt*trio->idx[iMOM]],ngt);
    }
    bcf_update_genotypes(args->hdr_out, rec, args->gt_arr, ngt*bcf_hdr_nsamples(args->hdr));
    return ret;
}

static int process_record(args_t *args, bcf1_t *rec)    // returns 1 to print rec, 0 to suppress printing
{
    int skip_site = 0;
    if ( rec->n_allele==1 || bcf_get_variant_types(rec)==VCF_REF )
    {
        args->nref_only++;
        skip_site = 1;
    }
    else if ( rec->n_allele > 64 )  // we use uint64_t bitmask in collect_stats()
    {
        args->nmany_als++;
        skip_site = 1;
    }
    else if ( args->filter && !test_filters(args,rec) )
    {
        args->stats.nfail++;
        skip_site = 1;
    }
    if ( skip_site ) return args->mode&MODE_DROP_SKIP ? 1 : 0;

    int ret = collect_stats(args, rec);
    if ( ret&HAS_MERR ) args->stats.nmerr++;
    if ( ret&HAS_MISS ) args->stats.nmiss++;
    if ( ret&HAS_GOOD ) args->stats.ngood++;
    if ( args->mode&MODE_COUNT ) return 0;
    if ( args->mode&MODE_DROP_ERR && ret&HAS_MERR ) return 0;
    if ( args->mode&MODE_DROP_MISS && ret&HAS_MISS ) return 0;

    if ( args->mode&MODE_ANNOTATE )
    {
        bcf_update_info_int32(args->hdr_out,rec,"MERR",&args->nmerr,1);
    }

    if ( args->mode&LIST_MODES )
    {
        if ( args->mode&MODE_LIST_ERR && ret&HAS_MERR ) return 1;
        if ( args->mode&MODE_LIST_MISS && ret&HAS_MISS ) return 1;
        if ( args->mode&MODE_LIST_GOOD && ret&HAS_GOOD ) return 1;
        return 0;
    }
    return 1;
}

static void print_stats(args_t *args)
{
    if ( !(args->mode & MODE_COUNT) ) return;

    FILE *log_fh = stderr;
    if ( args->mode==MODE_COUNT )
    {
        log_fh = strcmp("-",args->output_fname) ? fopen(args->output_fname,"w") : stdout;
        if ( !log_fh ) error("Error: cannot write to %s\n", args->output_fname);
    }

    fprintf(log_fh,"# Summary stats\n");
    fprintf(log_fh,"sites_ref_only\t%d\t# sites skipped because there was no ALT allele\n", args->nref_only);
    fprintf(log_fh,"sites_many_als\t%d\t# skipped because of too many ALT alleles\n", args->nmany_als);
    fprintf(log_fh,"sites_fail\t%d\t# skipped because of failed -i/-e filter\n", args->stats.nfail);
    fprintf(log_fh,"sites_no_GT\t%d\t# skipped because of absent FORMAT/GT field\n", args->stats.nno_gt);
    fprintf(log_fh,"sites_not_diploid\t%d\t# skipped because FORMAT/GT not formatted diploid\n", args->stats.nnot_diploid);
    fprintf(log_fh,"sites_missing\t%d\t# number of sites with at least one trio GT missing\n", args->stats.nmiss);
    fprintf(log_fh,"sites_merr\t%d\t# number of sites with at least one Mendelian error\n", args->stats.nmerr);
    fprintf(log_fh,"sites_good\t%d\t# number of sites with at least one good trio\n", args->stats.ngood);

    int i;
    fprintf(log_fh,"# Per-trio stats, each column corresponds to one trio. List of trios is below.\n");
    fprintf(log_fh,"# The meaning of per-trio stats is the same as described above, ngood_alt is\n");
    fprintf(log_fh,"# the number of good genotypes with at least one non-reference allele, and is\n");
    fprintf(log_fh,"# included in the ngood counter\n");
    fprintf(log_fh,"ngood");
    for (i=0; i<args->ntrio; i++) fprintf(log_fh,"\t%d",args->trio[i].stats.ngood);
    fprintf(log_fh,"\n");

    fprintf(log_fh,"ngood_alt");
    for (i=0; i<args->ntrio; i++) fprintf(log_fh,"\t%d",args->trio[i].stats.ngood_alt);
    fprintf(log_fh,"\n");

    fprintf(log_fh,"nmerr");
    for (i=0; i<args->ntrio; i++) fprintf(log_fh,"\t%d",args->trio[i].stats.nmerr);
    fprintf(log_fh,"\n");

    fprintf(log_fh,"nmissing");
    for (i=0; i<args->ntrio; i++) fprintf(log_fh,"\t%d",args->trio[i].stats.nmiss);
    fprintf(log_fh,"\n");

    fprintf(log_fh,"nfail");
    for (i=0; i<args->ntrio; i++) fprintf(log_fh,"\t%d",args->trio[i].stats.nfail);
    fprintf(log_fh,"\n");

    fprintf(log_fh,"# List of trios. Their ids are in the same order as the values listed in the stats lines above. For\n");
    fprintf(log_fh,"# example, the values for the first trio (id=1) and the third trio (id=3) are in the 2nd and the 4th\n");
    fprintf(log_fh,"# column and their stats can be obtained with the unix command\n");
    fprintf(log_fh,"#     cat stats.txt | grep ^n | cut -f1,2,4\n");
    fprintf(log_fh,"# TRIO\t[2]id\t[3]child\t[4]father\t[5]mother\n");
    for (i=0; i<args->ntrio; i++)
    {
        trio_t *trio = &args->trio[i];
        fprintf(log_fh,"TRIO\t%d\t%s\t%s\t%s\n", i+1,
                bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,trio->idx[iKID]),
                bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,trio->idx[iDAD]),
                bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,trio->idx[iMOM]));
    }

    if ( log_fh!=stderr && log_fh!=stdout && fclose(log_fh) ) error("Error: close failed for %s\n", args->output_fname);
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_fname = "-";
    args->record_cmd_line = 1;
    args->regions_overlap = 1;
    args->targets_overlap = 0;
    args->clevel = -1;

    static struct option loptions[] =
    {
        {"pfm",1,0,'p'},
        {"ped",1,0,'P'},
        {"mode",1,0,'m'},
        {"rules",1,0,1},
        {"rules-file",1,0,2},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"regions-overlap",required_argument,NULL,14},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"targets-overlap",required_argument,NULL,15},
        {"include",required_argument,0,'i'},
        {"exclude",required_argument,0,'e'},
        {"write-index",optional_argument,NULL,'W'},
        {0,0,0,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "?hp:P:m:o:O:i:e:t:T:r:R:W::",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 't': args->targets = optarg; break;
            case 'T': args->targets = optarg; args->targets_is_file = 1; break;
            case 'r': args->regions = optarg; break;
            case 'R': args->regions = optarg; args->regions_is_file = 1; break;
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
                      };
                      if ( optarg[1] )
                      {
                          args->clevel = strtol(optarg+1,&tmp,10);
                          if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                      }
                      break;
            case 'm':
                tmp = optarg;
                while ( *tmp )
                {
                    if ( *tmp=='+' || *tmp=='g' ) args->mode |= MODE_LIST_GOOD;
                    else if ( *tmp=='x' || *tmp=='e' ) args->mode |= MODE_LIST_ERR;
                    else if ( *tmp=='a' ) args->mode |= MODE_ANNOTATE;
                    else if ( *tmp=='d' ) args->mode |= MODE_DELETE;
                    else if ( *tmp=='c' ) args->mode |= MODE_COUNT;
                    else if ( *tmp=='E' ) args->mode |= MODE_DROP_ERR;
                    else if ( *tmp=='u' || *tmp=='m' ) args->mode |= MODE_LIST_MISS;
                    else if ( *tmp=='M' ) args->mode |= MODE_DROP_MISS;
                    else if ( *tmp=='s' ) args->mode |= MODE_LIST_SKIP;
                    else if ( *tmp=='S' ) args->mode |= MODE_DROP_SKIP;
                    else error("The argument \"%c\" not recognised: --mode %s\n", *tmp,optarg);
                    tmp++;
                }
                break;
            case 'P': args->ped_fname = optarg; break;
            case 'p': args->pfm = optarg; break;
            case  1 : args->rules_str = optarg; break;
            case  2 : args->rules_fname = optarg; break;
            case 'W':
                if (!(args->write_index = write_index_parse(optarg)))
                    error("Unsupported index format '%s'\n", optarg);
                break;
            case 'h':
            case '?':
            default: error("%s",usage_text()); break;
        }
    }

    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error("%s", usage_text()); }
    }
    else if ( optind+1!=argc ) error("%s", usage_text());
    else args->fname = argv[optind];

    if ( !args->ped_fname && !args->pfm ) error("Missing the -p or -P option\n");
    if ( args->ped_fname && args->pfm ) error("Expected only -p or -P option, not both\n");
    if ( !args->mode ) args->mode = MODE_COUNT;

    init_data(args);

    while ( bcf_sr_next_line(args->sr) )
    {
        bcf1_t *rec = bcf_sr_get_line(args->sr,0);
        if ( !process_record(args, rec) ) continue;
        if ( bcf_write(args->out_fh, args->hdr_out, rec)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
    }

    print_stats(args);
    destroy_data(args);

    return 0;
}

