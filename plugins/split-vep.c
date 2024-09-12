/* The MIT License

   Copyright (c) 2019-2024 Genome Research Ltd.

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
#include <unistd.h>     // for isatty
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/khash_str2int.h>
#include <regex.h>
#include "../bcftools.h"
#include "../filter.h"
#include "../convert.h"
#include "../cols.h"


// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define SELECT_TR_ALL       0
#define SELECT_TR_WORST     1
#define SELECT_TR_EXPR      2
#define SELECT_CSQ_ANY     -1

#define TR_OP_EQ  0   // =
#define TR_OP_NE  1   // !=
#define TR_OP_RE  2   // ~
#define TR_OP_NR  3   // !~

#define GENES_RESTRICT   0
#define GENES_PRIORITIZE 1

typedef struct
{
    regex_t *regex;
    char *type;
    int bcf_ht_type;
}
col2type_t;

typedef struct
{
    char *field;    // the name of the VEP field, e.g. Consequence,Gene,etc.
    char *tag;      // the name of the VCF tag: the annot_t.field with the -p prefix
    int idx;        // 0-based index within the VEP annotation string
    int type;       // annotation type, one of the BCF_HT_* types
    kstring_t str;  // annotation value, ready to pass to bcf_update_info_*
}
annot_t;

// for the use with --select TR
typedef struct
{
    char *field;        // field name to use
    int idx;            // field index
    int op;             // one of the TR_OP_* operators
    regex_t *regex;     // to be used with TR_OP_RE and TR_OP_NR
    char *value;        // to be used with TR_OP_EQ and TR_OP_NE
}
select_tr_t;

typedef struct
{
    convert_t *convert;
    filter_t *filter;
    int argc, filter_logic, regions_is_file, targets_is_file, list_hdr, record_cmd_line, clevel;
    int regions_overlap, targets_overlap;
    kstring_t kstr;
    char *filter_str,
        *vep_tag,       // the --annotation INFO tag to process
        *vep_format;    // the format of VEP fields, used for error reporting only
    char **argv, *output_fname, *fname, *regions, *targets, *format_str;
    int output_type;
    htsFile *fh_vcf;
    BGZF *fh_bgzf;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr, *hdr_out;
    int nfield;         // number of all available VEP fields
    char **field;       // list of all available VEP fields
    int nannot;         // number of requested fields
    annot_t *annot;     // requested fields
    int nscale;         // number of items in the severity scale
    char **scale;       // severity scale (list)
    int ncsq_str;       // the length of csq_str allocated by bcf_get_info_string()
    char *csq_str;      // the current bcf_get_info_string() result
    int csq_idx;        // the index of the Consequence field; for the --select CSQ option
    char *severity,     // the --severity scale option
        *select,        // the --select option
        *column_str,    // the --columns option
        *column_types,  // the --columns-types option
        *annot_prefix,  // the --annot-prefix option
        *genes_fname,   // the --gene-list option
        *gene_fields_str;   // the --gene-list-fields option
    int *gene_fields, ngene_fields; // the indices of --gene-list-fields names
    void *field2idx,    // VEP field name to index, used in initialization
        *csq2severity,  // consequence type to severity score
        *genes;         // hashed --gene-list genes
    cols_t *cols_tr,    // the current CSQ tag split into transcripts
        **cols_csq;     // the current CSQ transcripts split into fields
    int ncols_csq,mcols_csq;    // the number of cols_csq transcripts and the high-water mark
    int min_severity, max_severity;     // ignore consequences outside this severity range
    int drop_sites;                     // the -x, --drop-sites option
    int select_tr;                      // one of SELECT_TR_*
    select_tr_t tr_expr;                // with SELECT_TR_EXPR, see --select <FIELD><OP><VALUE>
    uint8_t *smpl_pass;                 // for filtering at sample level, used with -f
    int duplicate;              // the -d, --duplicate option is set
    char *all_fields_delim;     // the -A, --all-fields option is set
    float *farr;                // helper arrays for bcf_update_* functions
    int32_t *iarr, *list_tr;
    int niarr,miarr, nfarr,mfarr, mlist_tr;
    col2type_t *column2type;
    int ncolumn2type;
    int raw_vep_request;        // raw VEP tag requested and will need subsetting
    int allow_undef_tags;
    int genes_mode;             // --gene-list +FILE, one of GENES_* mode, prioritize or restrict
    int print_header;
    char *index_fn;
    int write_index;
}
args_t;

args_t args;

const char *about(void)
{
    return "Query structured annotations such as INFO/CSQ produced by VEP of INFO/BCSQ produced by bctools/csq.\n";
}

static const char *default_severity(void)
{
    return
        "# Default consequence substrings ordered in ascending order by severity.\n"
        "# Consequences with the same severity can be put on the same line in arbitrary order.\n"
        "# See also https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html\n"
        "intergenic\n"
        "feature_truncation feature_elongation\n"
        "regulatory\n"
        "TF_binding_site TFBS\n"
        "downstream upstream\n"
        "non_coding_transcript non_coding\n"
        "intron NMD_transcript\n"
        "non_coding_transcript_exon\n"
        "5_prime_utr 3_prime_utr\n"
        "coding_sequence mature_miRNA\n"
        "stop_retained start_retained synonymous\n"
        "incomplete_terminal_codon\n"
        "splice_region\n"
        "missense inframe protein_altering\n"
        "transcript_amplification\n"
        "exon_loss\n"
        "disruptive\n"
        "start_lost stop_lost stop_gained frameshift\n"
        "splice_acceptor splice_donor\n"
        "transcript_ablation\n";
}
static const char *default_column_types(void)
{
    return
        "# Default CSQ subfield types, unlisted fields are type String.\n"
        "# Note that the name search is done using regular expressions, with\n"
        "# \"^\" and \"$\" appended automatically\n"
        // These positions are in fact not integers but strings such as 8586-8599/9231
        //      "cDNA_position              Integer\n"
        //      "CDS_position               Integer\n"
        //      "Protein_position           Integer\n"
        "DISTANCE                   Integer\n"
        "STRAND                     Integer\n"
        "TSL                        Integer\n"
        "GENE_PHENO                 Integer\n"
        "HGVS_OFFSET                Integer\n"
        ".*_POPS                    String\n"   // e.g. MAX_AF_POPS
        "AF                         Float\n"
        ".*_AF                      Float\n"
        "MAX_AF_.*                  Float\n"
        "MOTIF_POS                  Integer\n"
        "MOTIF_SCORE_CHANGE         Float\n"
        "existing_InFrame_oORFs     Integer\n"
        "existing_OutOfFrame_oORFs  Integer\n"
        "existing_uORFs             Integer\n"
        "SpliceAI_pred_DP_.*        Integer\n"
        "SpliceAI_pred_DS_.*        Float\n";
}
static const char *usage_text(void)
{
    return
        "\n"
        "About: Query structured annotations such INFO/BCSQ or CSQ created by bcftools/csq or VEP. For more\n"
        "   more information and pointers see http://samtools.github.io/bcftools/howtos/plugin.split-vep.html\n"
        "Usage: bcftools +split-vep [Plugin Options]\n"
        "Plugin options:\n"
        "   -a, --annotation STR            INFO annotation to parse, CSQ by default or BCSQ when not present [CSQ]\n"
        "   -A, --all-fields DELIM          Output all fields replacing the -a tag (\"%CSQ\" by default) in the -f\n"
        "                                     filtering expression using the output field delimiter DELIM. This can be\n"
        "                                     \"tab\", \"space\" or an arbitrary string.\n"
        "   -c, --columns [LIST|-][:TYPE]   Extract the fields listed either as 0-based indexes or names, \"-\" to extract all\n"
        "                                     fields. See --columns-types for the defaults. Supported types are String/Str,\n"
        "                                     Integer/Int and Float/Real. Unlisted fields are set to String. Existing header\n"
        "                                     definitions will not be overwritten, remove first with `bcftools annotate -x`\n"
        "       --columns-types -|FILE      Pass \"-\" to print the default -c types or FILE to override the presets\n"
        "   -d, --duplicate                 Output per transcript/allele consequences on a new line rather rather than\n"
        "                                     as comma-separated fields on a single line\n"
        "   -f, --format STR                Create non-VCF output; similar to `bcftools query -f` but drops lines w/o consequence\n"
        "   -g, --gene-list [+]FILE         Consider only features listed in FILE, or prioritize if FILE is prefixed with \"+\"\n"
        "       --gene-list-fields LIST     Fields to match against by the -g list, by default gene names [SYMBOL,Gene,gene]\n"
        "   -H, --print-header              Print header, -HH to omit column indices\n"
        "   -l, --list                      Parse the VCF header and list the annotation fields\n"
        "   -p, --annot-prefix STR          Before doing anything else, prepend STR to all CSQ fields to avoid tag name conflicts\n"
        "   -s, --select TR:CSQ             Select transcripts to extract by type and/or consequence severity. (See also -S and -x.)\n"
        "                                     TR, transcript:   all,worst,primary,pick,mane,EXPRESSION [all]\n"
        "                                     CSQ, consequence: any,missense,missense+,etc [any]\n"
        "                                   Where the transcript EXPRESSION is of the form <FIELD><OPERATOR><VALUE>\n"
        "                                     FIELD:    field name (e.g. \"CANONICAL\")\n"
        "                                     OPERATOR: string comparison (=,!=), regex matching (~,!~)\n"
        "                                     VALUE:    required string value (e.g. \"YES\")\n"
        "                                   The TR presets are defined as follows\n"
        "                                     primary:  CANONICAL=YES\n"
        "                                     pick:     PICK=1\n"
        "                                     mane:     MANE_SELECT!=\"\"\n"
        "   -S, --severity -|FILE           Pass \"-\" to print the default severity scale or FILE to override\n"
        "                                     the default scale\n"
        "   -u, --allow-undef-tags          Print \".\" for undefined tags\n"
        "   -x, --drop-sites                Drop sites without consequences after -s is applied (the default with -f)\n"
        "   -X, --keep-sites                Do not drop sites without consequences (the default without -f)\n"
        "Common options:\n"
        "   -e, --exclude EXPR              Exclude sites and samples for which the expression is true\n"
        "   -i, --include EXPR              Include sites and samples for which the expression is true\n"
        "       --no-version                Do not append version and command line to the header\n"
        "   -o, --output FILE               Output file name [stdout]\n"
        "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n"
        "   -r, --regions REG               Restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE         Restrict to regions listed in a file\n"
        "       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n"
        "   -t, --targets REG               Similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE         Similar to -R but streams rather than index-jumps\n"
        "       --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n"
        "   -W, --write-index[=FMT]         Automatically index the output files [off]\n"
        "\n"
        "Examples:\n"
        "   # List available fields of the INFO/CSQ or INFO/BCSQ annotation\n"
        "   bcftools +split-vep -l file.vcf.gz\n"
        "\n"
        "   # List available fields of INFO/BCSQ when both INFO/CSQ and INFO/BCSQ are present\n"
        "   bcftools +split-vep -l file.vcf.gz -a BCSQ\n"
        "\n"
        "   # List the default severity scale\n"
        "   bcftools +split-vep -S -\n"
        "\n"
        "   # Extract Consequence, IMPACT and gene SYMBOL of the most severe consequence into\n"
        "   # INFO annotations starting with the prefix \"vep\". For brevity, the columns can\n"
        "   # be given also as 0-based indexes\n"
        "   bcftools +split-vep -c Consequence,IMPACT,SYMBOL -s worst -p vep file.vcf.gz\n"
        "   bcftools +split-vep -c 1-3 -s worst -p vep file.vcf.gz\n"
        "\n"
        "   # Same as above but use the text output of the \"bcftools query\" format\n"
        "   bcftools +split-vep -s worst -f '%CHROM %POS %Consequence %IMPACT %SYMBOL\\n' file.vcf.gz\n"
        "\n"
        "   # Print all subfields (tab-delimited) in place of %CSQ, each consequence on a new line\n"
        "   bcftools +split-vep -f '%CHROM %POS %CSQ\\n' -d -A tab file.vcf.gz\n"
        "\n"
        "   # Extract gnomAD_AF subfield into a new INFO/gnomAD_AF annotation of Type=Float so that\n"
        "   # numeric filtering can be used.\n"
        "   bcftools +split-vep -c gnomAD_AF:Float file.vcf.gz -i'gnomAD_AF<0.001'\n"
        "\n"
        "   # Similar to above, but add the annotation only if the consequence severity is missense\n"
        "   # or equivalent. In order to drop sites with different consequences completely, provide\n"
        "   # the -x option. See the online documentation referenced above for more examples.\n"
        "   bcftools +split-vep -c gnomAD_AF:Float -s :missense    file.vcf.gz\n"
        "   bcftools +split-vep -c gnomAD_AF:Float -s :missense -x file.vcf.gz\n"
        "\n"
        "   See also http://samtools.github.io/bcftools/howtos/plugin.split-vep.html\n"
        "\n";
}

static void destroy_annot(args_t *args)
{
    int i;
    for (i=0; i<args->nannot; i++)
    {
        annot_t *ann = &args->annot[i];
        free(ann->field);
        free(ann->tag);
        free(ann->str.s);
    }
    free(args->annot);
    args->annot = NULL;
    args->nannot = 0;
}

static void expand_csq_expression(args_t *args, kstring_t *str)
{
    if ( !args->all_fields_delim ) return;

    str->l = 0;
    kputc('%',str);
    kputs(args->vep_tag,str);
    char *ptr = strstr(args->format_str,str->s);
    if ( !ptr ) return;
    char *end = ptr + str->l, tmp = *end;
    if ( isalnum(tmp) || tmp=='_' || tmp=='.' ) return;
    *end = 0;

    str->l = 0;
    kputsn(args->format_str, ptr - args->format_str, str);

    int i;
    for (i=0; i<args->nfield; i++)
    {
        if ( i>0 ) kputs(args->all_fields_delim, str);
        kputc('%', str);
        kputs(args->field[i], str);
    }

    *end = tmp;
    kputs(end, str);

    free(args->format_str);
    args->format_str = str->s;
    str->l = str->m = 0;
    str->s = NULL;
}

static void init_column2type(args_t *args)
{
    int i, ntype = 0;
    char **type = NULL;
    if ( args->column_types && strcmp("-",args->column_types) )
        type = hts_readlines(args->column_types, &ntype);
    else
    {
        char *str = strdup(default_column_types());
        char *beg = str;
        while ( *beg )
        {
            char *end = beg;
            while ( *end && *end!='\n' ) end++;
            char tmp = *end;
            *end = 0;
            ntype++;
            type = (char**) realloc(type,sizeof(*type)*ntype);
            type[ntype-1] = strdup(beg);
            if ( !tmp ) break;
            beg = end+1;
        }
        free(str);
    }
    if ( !type || !ntype ) error("Failed to parse the column types\n");
    kstring_t tmp = {0,0,0};
    for (i=0; i<ntype; i++)
    {
        if ( type[i][0]=='#' ) continue;
        tmp.l = 0;
        kputc('^',&tmp);
        char *ptr = type[i];
        while ( *ptr && !isspace(*ptr) ) ptr++;
        if ( !*ptr ) error("Error: failed to parse the column type \"%s\"\n",type[i]);
        kputsn(type[i],ptr-type[i],&tmp);
        kputc('$',&tmp);
        while ( *ptr && isspace(*ptr) ) ptr++;
        if ( !*ptr ) error("Error: failed to parse the column type \"%s\"\n",type[i]);
        args->ncolumn2type++;
        args->column2type = (col2type_t*) realloc(args->column2type,sizeof(*args->column2type)*args->ncolumn2type);
        col2type_t *ct = &args->column2type[args->ncolumn2type-1];
        ct->regex = (regex_t *) malloc(sizeof(regex_t));
        if ( regcomp(ct->regex, tmp.s, REG_NOSUB) )
                error("Error: fail to compile the column type regular expression \"%s\": %s\n", tmp.s,type[i]);
        ct->bcf_ht_type = -1;
        if ( !strcmp(ptr,"Float") ) ct->bcf_ht_type = BCF_HT_REAL;
        else if ( !strcmp(ptr,"Integer") ) ct->bcf_ht_type = BCF_HT_INT;
        else if ( !strcmp(ptr,"Flag") ) ct->bcf_ht_type = BCF_HT_FLAG;
        else if ( !strcmp(ptr,"String") ) ct->bcf_ht_type = BCF_HT_STR;
        if ( ct->bcf_ht_type==-1 ) error("Error: the column type \"%s\" is not supported: %s\n",ptr,type[i]);
        ct->type = strdup(ptr);
    }
    free(tmp.s);
    if ( !args->ncolumn2type ) error("Failed to parse the column types\n");
    for (i=0; i<ntype; i++) free(type[i]);
    free(type);
}
static void destroy_column2type(args_t *args)
{
    int i;
    for (i=0; i<args->ncolumn2type; i++)
    {
        if ( args->column2type[i].regex ) regfree(args->column2type[i].regex);
        free(args->column2type[i].regex);
        free(args->column2type[i].type);
    }
    free(args->column2type);
    args->ncolumn2type = 0;
    args->column2type = NULL;
}
static const char *get_column_type(args_t *args, char *field, int *type)
{
    if ( !args->column2type ) init_column2type(args);
    int i;
    for (i=0; i<args->ncolumn2type; i++)
    {
        int match = regexec(args->column2type[i].regex, field, 0,NULL,0) ? 0 : 1;
        if ( match )
        {
            *type = args->column2type[i].bcf_ht_type;
            return args->column2type[i].type;
        }
    }
    *type = BCF_HT_STR;
    return "String";
}

// Is the tag "field" present in the -f query string, such as '%CHROM %POS %field\n'?
static int query_has_field(char *fmt, char *field, kstring_t *str)
{
    str->l = 0;
    kputc('%',str);
    kputs(field,str);
    char end, *ptr = fmt;
    while ( ptr )
    {
        ptr = strstr(ptr,str->s);
        if ( !ptr ) return 0;
        end = ptr[str->l];
        if ( isalnum(end) || end=='_' || end=='.' )
        {
            ptr++;
            continue;
        }
        break;
    }
    return 1;
}
/**
   The valid_tag array was generated with
        perl -le '@v = (split(//,q[_.]),"a"..."z","A"..."Z","0"..."9"); @a = (0) x 256; foreach $c (@v) { $a[ord($c)] = 1; } print join(", ",@a)' | fold -w 48
*/
static const uint8_t valid_tag[256] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1,
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
static char *sanitize_field_name(const char *str)
{
    if ( !strcmp(str,"1000G") ) return strdup(str);
    char *tmp;
    if ( str[0]=='.' || (str[0]>='0' && str[0]<='9') )
    {
        // the field starts with an invalid character, prefix with underscore
        int len = 1 + strlen(str);
        tmp = (char*)malloc(len+1);
        tmp[0] = '_';
        memcpy(tmp+1,str,len);
    }
    else tmp = strdup(str);
    char *out = tmp;
    while ( *tmp )
    {
        if ( !valid_tag[(uint8_t)*tmp] ) *tmp = '_';
        tmp++;
    }
    return out;
}
char *strdup_annot_prefix(args_t *args, const char *str)
{
    char *out;
    if ( !args->annot_prefix )
    {
        out = sanitize_field_name(str);
        return out;
    }
    int str_len = strlen(str);
    int prefix_len = strlen(args->annot_prefix);
    char *tmp = calloc(str_len+prefix_len+1,1);
    memcpy(tmp,args->annot_prefix,prefix_len);
    memcpy(tmp+prefix_len,str,str_len);
    out = sanitize_field_name(tmp);
    free(tmp);
    return out;
}
static void init_gene_list(args_t *args)
{
    int i,j,ngene = 0;
    args->genes_mode = args->genes_fname[0]=='+' ? GENES_PRIORITIZE : GENES_RESTRICT;
    char *fname = args->genes_fname[0]=='+' ? args->genes_fname+1 : args->genes_fname;
    char **gene = hts_readlines(fname, &ngene);
    if ( !ngene || !gene ) error("Could not read the file provided with --gene-list %s\n",fname);
    args->genes = khash_str2int_init();
    for (i=0; i<ngene; i++)
    {
        char *tmp = strdup(gene[i]);
        if ( khash_str2int_set(args->genes,tmp,1) < 0 ) free(tmp);  // duplicate gene name
        free(gene[i]);
    }
    free(gene);

    if ( !args->gene_fields_str ) args->gene_fields_str = "SYMBOL,Gene,gene";
    gene = hts_readlist(args->gene_fields_str,0,&ngene);
    args->gene_fields = (int*)malloc(ngene*sizeof(*args->gene_fields));
    for (i=0,j=0; i<ngene; i++)
    {
        char *tmp = strdup_annot_prefix(args, gene[i]);
        if ( khash_str2int_get(args->field2idx,tmp,&args->gene_fields[j])==0 ) j++;
        free(tmp);
        free(gene[i]);
    }
    free(gene);
    args->ngene_fields = j;
    if ( !args->ngene_fields ) error("None of the \"%s\" fields is present in INFO/%s\n",args->gene_fields_str,args->vep_tag);
}

// The program was requested to create a text output as with `bcftools query -f`. For this we need to
// determine the fields to be extracted from the formatting expression. Any subfields not listed in the
// header will be added to column_str so that they can be added to the header
static void parse_format_str(args_t *args)
{
    int i;
    kstring_t str = {0,0,0};

    // Special case: -A was given, extract all fields, for this the -a tag (%CSQ) must be present
    if ( args->all_fields_delim ) expand_csq_expression(args, &str);

    for (i=0; i<args->nfield; i++)
    {
        if ( !query_has_field(args->format_str,args->field[i],&str) ) continue;

        int tag_id = bcf_hdr_id2int(args->hdr, BCF_DT_ID, args->field[i]);
        if ( bcf_hdr_idinfo_exists(args->hdr,BCF_HL_INFO,tag_id) )
            fprintf(stderr,"Note: ambiguous key %%%s; using the %s subfield of %s, not the INFO/%s tag\n", args->field[i],args->field[i],args->vep_tag,args->field[i]);

        int olen = args->column_str ? strlen(args->column_str) : 0;
        int nlen = strlen(args->field[i]);
        args->column_str = (char*)realloc(args->column_str, olen + nlen + 2);
        if ( olen )
        {
            memcpy(args->column_str+olen,",",1);
            olen++;
        }
        memcpy(args->column_str+olen,args->field[i],nlen);
        args->column_str[olen+nlen] = 0;
    }
    if ( query_has_field(args->format_str,args->vep_tag,&str) ) args->raw_vep_request = 1;
    free(str.s);
}

// The program was requested to extract one or more columns via -c. It can contain names,  0-based indexes or ranges of indexes
static void parse_column_str(args_t *args)
{
    if ( args->nannot ) return; // already called from parse_filter_str
    int i,j;
    int *column = NULL;
    int *types  = NULL;
    if ( !strcmp("-",args->column_str) )    // all subfields
    {
        free(args->column_str);
        kstring_t str = {0,0,0};
        ksprintf(&str,"0-%d",args->nfield-1);
        args->column_str = str.s;
    }
    char *ep = args->column_str;
    while ( *ep )
    {
        char *tp, *bp = ep;
        while ( *ep && *ep!=',' ) ep++;
        char keep = *ep;
        *ep = 0;
        int type = -1;
        int idx_beg, idx_end;
        if ( !strcmp("-",bp) )
        {
            kstring_t str = {0,0,0};
            ksprintf(&str,"0-%d",args->nfield-1);
            if ( keep ) ksprintf(&str,",%s",ep+1);
            free(args->column_str);
            args->column_str = str.s;
            ep = str.s;
            continue;
        }
        char *tmp = strdup_annot_prefix(args, bp);  // replace characters disallowed in tag names into underscores
        if ( khash_str2int_get(args->field2idx, bp, &idx_beg)==0 || khash_str2int_get(args->field2idx, tmp, &idx_beg)==0 )
        {
            // either the original or sanitized version of the tag exists
            idx_end = idx_beg;
        }
        else if ( (tp=strrchr(bp,':')) )    // notice this requests the last occurrence of ':'
        {
            // there is a colon in the original string, expecting type specification
            *tp = 0;
            if ( khash_str2int_get(args->field2idx, bp, &idx_beg)!=0 )
            {
                // even removing the part after the last colon does not give a valid tag
                error("No such column: \"%s\"\n", bp);
            }
            idx_end = idx_beg;
            *tp = ':';
            if ( !strcasecmp(tp+1,"string") || !strcasecmp(tp+1,"str") ) type = BCF_HT_STR;
            else if ( !strcasecmp(tp+1,"float") || !strcasecmp(tp+1,"real") ) type = BCF_HT_REAL;
            else if ( !strcasecmp(tp+1,"integer") || !strcasecmp(tp+1,"int") ) type = BCF_HT_INT;
            else if ( !strcasecmp(tp+1,"flag") ) type = BCF_HT_FLAG;
            else error("The type \"%s\" nor the column \"%s\" recognised\n", tp+1,bp);
        }
        else
        {
            char *mp;
            idx_beg = strtol(bp,&mp,10);
            if ( !*mp ) idx_end = idx_beg;
            else if ( *mp=='-' )
                idx_end = strtol(mp+1,&mp,10);
            if ( *mp )
            {
                if ( *mp==':' )
                {
                    idx_end = idx_beg;
                    if ( !strcasecmp(mp+1,"string") || !strcasecmp(mp+1,"str") ) type = BCF_HT_STR;
                    else if ( !strcasecmp(mp+1,"float") || !strcasecmp(mp+1,"real") ) type = BCF_HT_REAL;
                    else if ( !strcasecmp(mp+1,"integer") || !strcasecmp(mp+1,"int") ) type = BCF_HT_INT;
                    else if ( !strcasecmp(mp+1,"flag") ) type = BCF_HT_FLAG;
                    else error("The type \"%s\" nor the column \"%s\" recognised\n", mp+1,bp);
                }
                else if ( !strcmp(bp,args->vep_tag) )
                {
                    free(tmp);
                    args->raw_vep_request = 1;
                    if ( !keep ) break;
                    ep++;
                    continue;
                }
                else
                    error("No such column: \"%s\"\n", bp);
            }
        }
        free(tmp);

        i = args->nannot;
        args->nannot += idx_end - idx_beg + 1;
        column = (int*)realloc(column,args->nannot*sizeof(*column));
        types  = (int*)realloc(types,args->nannot*sizeof(*types));
        for (j=idx_beg; j<=idx_end; j++)
        {
            if ( j >= args->nfield ) error("The index is too big: %d\n", j);
            column[i] = j;
            types[i]  = type;
            i++;
        }
        if ( !keep ) break;
        ep++;
    }

    // Prune duplicates
    for (i=0; i<args->nannot; i++)
    {
        for (j=0; j<i; j++)
            if ( !strcmp(args->field[column[i]],args->field[column[j]]) ) break;
        if ( i==j ) continue;           // unique tag, no action needed
        args->nannot--;
        if ( i==args->nannot ) break;   // the last one is to be skipped, we are done
        memmove(&column[i],&column[i+1],sizeof(*column)*(args->nannot-i));
        i--;
    }

    // Now initizalize each annotation, add each column to the VCF header, and reconstruct
    // the column_str in case it will be needed later
    free(args->column_str);
    kstring_t str = {0,0,0};
    args->annot = (annot_t*)calloc(args->nannot,sizeof(*args->annot));
    for (i=0; i<args->nannot; i++)
    {
        annot_t *ann = &args->annot[i];
        ann->type = types[i];
        ann->idx = j = column[i];
        ann->field = strdup(args->field[j]);
        ann->tag = strdup(args->field[j]);
        const char *type = "String";
        if ( ann->type==BCF_HT_REAL ) type = "Float";
        else if ( ann->type==BCF_HT_INT ) type = "Integer";
        else if ( ann->type==BCF_HT_FLAG ) type = "Flag";
        else if ( ann->type==BCF_HT_STR ) type = "String";
        else if ( ann->type==-1 ) type = get_column_type(args, args->field[j], &ann->type);
        bcf_hdr_printf(args->hdr_out, "##INFO=<ID=%s,Number=.,Type=%s,Description=\"The %s field from INFO/%s\">", ann->tag,type,ann->field,args->vep_tag);
        if ( str.l ) kputc(',',&str);
        kputs(ann->tag,&str);
    }
    args->column_str = str.s;
    if ( args->raw_vep_request && args->select_tr==SELECT_TR_ALL ) args->raw_vep_request = 0;
    if ( args->raw_vep_request )
    {
        args->nannot++;
        args->annot = (annot_t*)realloc(args->annot,args->nannot*sizeof(*args->annot));
        annot_t *ann = &args->annot[args->nannot-1];
        ann->type  = BCF_HT_STR;
        ann->idx   = -1;
        ann->field = strdup(args->vep_tag);
        ann->tag   = strdup(args->vep_tag);
        memset(&ann->str,0,sizeof(ann->str));
    }
    free(column);
    free(types);
    destroy_column2type(args);

    if ( bcf_hdr_sync(args->hdr_out)<0 )
        error_errno("[%s] Failed to update header", __func__);
}

// Init filters. When neither -c,--columns nor -f,--format contains a VEP subfield used in the
// filtering expression, and the subfield is not defined in the VCF header, then the filtering
// will throw an error. We can be smart, detect these failures and such undefined subfields
// as if the user passed them via the -c option.
static void parse_filter_str(args_t *args)
{
    args->filter = filter_parse(args->hdr_out, args->filter_str);
    if ( !args->filter ) error(NULL);     // this type of error would have been reported
    int ret = filter_status(args->filter);
    if ( ret!=FILTER_OK )
    {
        if ( ret!=FILTER_ERR_UNKN_TAGS ) error(NULL); // this type of error would have been reported

        // add the undefined tags to the -c string
        int ntags, i,j;
        const char **tags = filter_list_undef_tags(args->filter, &ntags);
        kstring_t str;
        str.s = args->column_str;
        str.l = str.m = str.s ? strlen(str.s) : 0;
        for (i=0; i<ntags; i++)
        {
            if ( khash_str2int_get(args->field2idx,tags[i],&j)!=0 )
                error("Error: the tag \"%s\" is not defined in the VCF header or in INFO/%s\n",tags[i],args->vep_tag);
            if ( str.l ) kputc(',',&str);
            kputs(tags[i],&str);
        }
        args->column_str = str.s;
        parse_column_str(args);
        filter_destroy(args->filter);
        args->filter = filter_init(args->hdr_out, args->filter_str);
    }
    int ntags, i;
    const char **tags = filter_list_used_tags(args->filter, &ntags);
    for (i=0; i<ntags; i++)
        if ( !strncmp("INFO/",tags[i],5) && !strcmp(tags[i]+5,args->vep_tag) ) args->raw_vep_request = 1;
}

void init_select_tr_expr(args_t *args, const char *cnst_expr)
{
    char *expr = strdup(cnst_expr);
    char *ptr = expr;
    while ( *ptr )
    {
        if ( *ptr=='=' )
        {
            *ptr = 0;
            args->tr_expr.field = strdup_annot_prefix(args,expr);
            *ptr = '=';
            int is_quoted = ( ptr[1]=='"' && ptr[strlen(ptr)-1]=='"' ) ? 1 : 0;
            args->tr_expr.value = strdup(ptr+1+is_quoted);
            if ( is_quoted ) args->tr_expr.value[strlen(args->tr_expr.value)-1] = 0;
            args->tr_expr.op = TR_OP_EQ;
            break;
        }
        else if ( *ptr=='~' )
        {
            *ptr = 0;
            args->tr_expr.field = strdup_annot_prefix(args,expr);
            *ptr = '~';
            int is_quoted = ( ptr[1]=='"' && ptr[strlen(ptr)-1]=='"' ) ? 1 : 0;
            args->tr_expr.value = strdup(ptr+1+is_quoted);
            if ( is_quoted ) args->tr_expr.value[strlen(args->tr_expr.value)-1] = 0;
            args->tr_expr.regex = (regex_t *) malloc(sizeof(regex_t));
            if ( regcomp(args->tr_expr.regex, args->tr_expr.value, REG_NOSUB) )
                error("Error: fail to compile the regular expression \"%s\"\n", args->tr_expr.value);
            args->tr_expr.op = TR_OP_RE;
        }
        else if ( *ptr=='!' && ptr[1]=='=' )
        {
            *ptr = 0;
            args->tr_expr.field = strdup_annot_prefix(args,expr);
            *ptr = '!';
            int is_quoted = ( ptr[2]=='"' && ptr[strlen(ptr)-1]=='"' ) ? 1 : 0;
            args->tr_expr.value = strdup(ptr+2+is_quoted);
            if ( is_quoted ) args->tr_expr.value[strlen(args->tr_expr.value)-1] = 0;
            args->tr_expr.op = TR_OP_NE;
            break;
        }
        else if ( *ptr=='!' && ptr[1]=='~' )
        {
            *ptr = 0;
            args->tr_expr.field = strdup_annot_prefix(args,expr);
            *ptr = '!';
            int is_quoted = ( ptr[2]=='"' && ptr[strlen(ptr)-1]=='"' ) ? 1 : 0;
            args->tr_expr.value = strdup(ptr+2+is_quoted);
            if ( is_quoted ) args->tr_expr.value[strlen(args->tr_expr.value)-1] = 0;
            args->tr_expr.regex = (regex_t *) malloc(sizeof(regex_t));
            if ( regcomp(args->tr_expr.regex, args->tr_expr.value, REG_NOSUB) )
                error("Error: fail to compile the regular expression \"%s\"\n", args->tr_expr.value);
            args->tr_expr.op = TR_OP_NR;
            break;
        }
        ptr++;
    }
    if ( !args->tr_expr.field ) error("Could not parse the expression: -s %s\n", cnst_expr);
    if ( khash_str2int_get(args->field2idx,args->tr_expr.field,&args->tr_expr.idx)!=0 )
        error("The field \"%s\" was requested via \"%s\" but it is not present in INFO/%s: %s\n",args->tr_expr.field,expr,args->vep_tag,args->vep_format);
    free(expr);
    args->select_tr = SELECT_TR_EXPR;
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
    args->hdr_out = bcf_hdr_dup(args->hdr);

    // Check which one to use: BCSQ or CSQ
    if ( !args->vep_tag )
    {
        int has_CSQ  = bcf_hdr_idinfo_exists(args->hdr,BCF_HL_INFO,bcf_hdr_id2int(args->hdr,BCF_DT_ID,"CSQ"));
        int has_BCSQ = bcf_hdr_idinfo_exists(args->hdr,BCF_HL_INFO,bcf_hdr_id2int(args->hdr,BCF_DT_ID,"BCSQ"));
        if ( has_CSQ && has_BCSQ ) fprintf(stderr,"Warning: both INFO/CSQ and INFO/BCSQ exist, using INFO/CSQ\n");
        if ( !has_CSQ && !has_BCSQ ) error("Error: Neither INFO/CSQ nor INFO/BCSQ was found in the header\n");
        if ( has_CSQ ) args->vep_tag = "CSQ";
        else if ( has_BCSQ ) args->vep_tag = "BCSQ";
    }

    // Parse the header CSQ line, must contain Description with "Format: ..." declaration
    bcf_hrec_t *hrec = bcf_hdr_get_hrec(args->hdr, BCF_HL_INFO, NULL, args->vep_tag, NULL);
    if ( !hrec ) error("The tag INFO/%s not found in the header\n", args->vep_tag);
    int ret = bcf_hrec_find_key(hrec, "Description");
    if ( ret<0 ) error("No \"Description\" field was found for the tag INFO/%s in the header\n", args->vep_tag);
    char *format = strstr(hrec->vals[ret], "Format: ");
    if ( !format ) error("Expected \"Format: \" substring in the header INFO/%s/Description, found: %s\n", args->vep_tag,hrec->vals[ret]);
    format += 8;
    args->vep_format = strdup(format);
    char *ep = args->vep_format + strlen(args->vep_format) - 1;
    if ( *ep=='"' )  *ep = 0;
    ep = format;
    while ( *ep )
    {
        char *bp = ep;
        while ( *ep && *ep!='|' && *ep!='"' && *ep!='(' ) ep++; // don't include bracket in "pos(1-based)"
        char tmp = *ep;
        *ep = 0;
        args->nfield++;
        args->field = (char**)realloc(args->field,args->nfield*sizeof(*args->field));
        args->field[args->nfield-1] = strdup_annot_prefix(args,bp);
        if ( !tmp ) break;
        *ep = tmp;
        while ( *ep && *ep!='|' ) ep++;     // skip bracket in "pos(1-based)"
        if ( *ep ) ep++;
    }
    if ( !args->nfield ) error("Could not parse Description of INFO/%s: %s\n", args->vep_tag,hrec->vals[ret]);
    args->field2idx = khash_str2int_init();
    int i;
    for (i=0; i<args->nfield; i++)
    {
        if ( khash_str2int_has_key(args->field2idx, args->field[i]) )
        {
            fprintf(stderr,"Warning: duplicate INFO/%s key \"%s\"\n", args->vep_tag,args->field[i]);
            continue;
        }
        khash_str2int_set(args->field2idx, args->field[i], i);
    }

    // Severity scale
    kstring_t str = {0,0,0};
    args->csq2severity = khash_str2int_init();
    int severity = 0;
    if ( args->severity )
    {
        kstring_t tmp = {0,0,0};
        htsFile *fp = hts_open(args->severity,"r");
        if ( !fp ) error("Cannot read %s\n", args->severity);
        while ( hts_getline(fp, KS_SEP_LINE, &tmp) > 0 )
        {
            kputs(tmp.s, &str);
            kputc('\n', &str);
        }
        free(tmp.s);
    }
    else
        kputs(default_severity(),&str);
    ep = str.s;
    while ( *ep )
    {
        if ( *ep=='#' )
        {
            while ( *ep && *ep!='\n' ) { *ep = tolower(*ep); ep++; }
            if ( !*ep ) break;
            ep++;
            continue;
        }
        char *bp = ep;
        while ( *ep && !isspace(*ep) ) { *ep = tolower(*ep); ep++; }
        char tmp = *ep;
        *ep = 0;
        args->nscale++;
        args->scale = (char**) realloc(args->scale,args->nscale*sizeof(*args->scale));
        args->scale[args->nscale-1] = strdup(bp);
        if ( !khash_str2int_has_key(args->csq2severity,args->scale[args->nscale-1]) )
            khash_str2int_set(args->csq2severity,args->scale[args->nscale-1], severity);
        if ( !tmp ) break;
        if ( tmp=='\n' ) severity++;
        ep++;
        while ( *ep && isspace(*ep) ) ep++;
    }

    // Transcript and consequence selection
    if ( !args->select ) args->select = "all:any";
    cols_t *cols = cols_split(args->select, NULL, ':');
    char *sel_tr  = cols->off[0][0] ? cols->off[0] : "all";
    char *sel_csq = cols->n==2 && cols->off[1][0] ? cols->off[1] : "any";

    // ... transcript selection
    if ( !strcasecmp(sel_tr,"all") ) args->select_tr = SELECT_TR_ALL;
    else if ( !strcasecmp(sel_tr,"worst") ) args->select_tr = SELECT_TR_WORST;
    else if ( !strcasecmp(sel_tr,"primary") ) init_select_tr_expr(args,"CANONICAL=YES");
    else if ( !strcasecmp(sel_tr,"pick") ) init_select_tr_expr(args,"PICK=1");
    else if ( !strcasecmp(sel_tr,"mane") ) init_select_tr_expr(args,"MANE_SELECT!=\"\"");
    else init_select_tr_expr(args,sel_tr);

    // ... consequence selection
    if ( !strcasecmp(sel_csq,"any") ) { args->min_severity = args->max_severity = SELECT_CSQ_ANY; }     // to avoid unnecessary lookups
    else
    {
        int len = strlen(sel_csq);
        int severity, modifier = '=';
        if ( sel_csq[len-1]=='+' ) { modifier = '+'; sel_csq[len-1] = 0; }
        else if ( sel_csq[len-1]=='-' ) { modifier = '-'; sel_csq[len-1] = 0; }
        if ( khash_str2int_get(args->csq2severity, sel_csq, &severity)!=0 )
            error("Error: the consequence \"%s\" is not recognised. Run \"bcftools +split-vep -S ?\" to see the default list.\n", sel_csq);
        if ( modifier=='=' ) { args->min_severity = severity; args->max_severity = severity; }
        else if ( modifier=='+' ) { args->min_severity = severity; args->max_severity = INT_MAX; }
        else if ( modifier=='-' ) { args->min_severity = 0; args->max_severity = severity; }
    }
    cols_destroy(cols);

    // The "Consequence" column to determine severity for filtering. The name of this column is hardwired for now, both VEP and bt/csq use the same name
    char *tmp = strdup_annot_prefix(args,"Consequence");
    if ( khash_str2int_get(args->field2idx,tmp,&args->csq_idx)!=0 )
        error("The field \"Consequence\" is not present in INFO/%s: %s\n", args->vep_tag,hrec->vals[ret]);
    free(tmp);

    if ( args->format_str ) parse_format_str(args);    // Text output, e.g. bcftools +split-vep -f '%Consequence\n'
    if ( args->filter_str ) parse_filter_str(args);
    if ( args->column_str ) parse_column_str(args);    // The --columns option was given, update the header
    if ( args->format_str )
    {
        if ( !args->column_str && !args->select ) error("Error: No %s field selected in the formatting expression and -s not given: a typo?\n",args->vep_tag);
        args->convert = convert_init(args->hdr_out, NULL, 0, args->format_str);
        if ( !args->convert ) error("Could not parse the expression: %s\n", args->format_str);
        if ( args->allow_undef_tags ) convert_set_option(args->convert, allow_undef_tags, 1);
        if ( args->print_header>1 ) convert_set_option(args->convert, no_hdr_indices, 1);
        convert_set_option(args->convert, force_newline, 1);
    }
    if ( args->genes_fname ) init_gene_list(args);

    int max_unpack = BCF_UN_SHR;
    if ( args->convert ) max_unpack |= convert_max_unpack(args->convert);
    if ( args->filter ) max_unpack |= filter_max_unpack(args->filter);
    if ( !args->format_str ) max_unpack |= BCF_UN_FMT;      // don't drop FMT fields on VCF input when VCF/BCF is output
    args->sr->max_unpack = max_unpack;
    if ( args->convert && (max_unpack & BCF_UN_FMT) )
        convert_set_option(args->convert, subset_samples, &args->smpl_pass);

    free(str.s);
}
static void destroy_data(args_t *args)
{
    free(args->list_tr);
    if ( args->tr_expr.regex ) regfree(args->tr_expr.regex);
    free(args->tr_expr.regex);
    free(args->tr_expr.field);
    free(args->tr_expr.value);
    free(args->vep_format);
    free(args->farr);
    free(args->iarr);
    free(args->kstr.s);
    free(args->column_str);
    free(args->format_str);
    cols_destroy(args->cols_tr);
    int i;
    for (i=0; i<args->mcols_csq; i++)
        cols_destroy(args->cols_csq[i]);
    free(args->cols_csq);
    for (i=0; i<args->nscale; i++) free(args->scale[i]);
    free(args->scale);
    for (i=0; i<args->nfield; i++) free(args->field[i]);
    free(args->field);
    destroy_annot(args);
    free(args->gene_fields);
    if ( args->genes ) khash_str2int_destroy_free(args->genes);
    if ( args->field2idx ) khash_str2int_destroy(args->field2idx);
    if ( args->csq2severity ) khash_str2int_destroy(args->csq2severity);
    bcf_sr_destroy(args->sr);
    bcf_hdr_destroy(args->hdr_out);
    free(args->csq_str);
    if ( args->filter ) filter_destroy(args->filter);
    if ( args->convert ) convert_destroy(args->convert);
    if ( args->fh_vcf )
    {
        if ( args->write_index )
        {
            if ( bcf_idx_save(args->fh_vcf)<0 )
            {
                if ( hts_close(args->fh_vcf)!=0 ) error("Error: close failed .. %s\n", args->output_fname?args->output_fname:"stdout");
                error("Error: cannot write to index %s\n", args->index_fn);
            }
            free(args->index_fn);
        }
        if ( hts_close(args->fh_vcf)!=0 ) error("Error: close failed .. %s\n",args->output_fname);
    }
    if ( args->fh_bgzf && bgzf_close(args->fh_bgzf)!=0 ) error("Error: close failed .. %s\n",args->output_fname);
    free(args);
}
static void list_header(args_t *args)
{
    int i;
    for (i=0; i<args->nfield; i++) printf("%d\t%s\n", i,args->field[i]);
}

static void csq_to_severity(args_t *args, char *csq, int *min_severity, int *max_severity, int exact_match)
{
    *min_severity = INT_MAX;
    *max_severity = -1;
    char *ep = csq;
    while ( *ep )
    {
        char *bp = ep;
        while ( *ep && *ep!='&' ) { *ep = tolower(*ep); ep++; }
        char tmp = *ep;
        *ep = 0;

        int i, severity = -1;
        if ( khash_str2int_get(args->csq2severity, bp, &severity)!=0 )
        {
            for (i=0; i<args->nscale; i++)
                if ( strstr(bp,args->scale[i]) ) break;

            if ( i!=args->nscale )
                khash_str2int_get(args->csq2severity, args->scale[i], &severity);
            else
                severity = args->nscale + 1;

            args->nscale++;
            args->scale = (char**) realloc(args->scale,args->nscale*sizeof(*args->scale));
            args->scale[args->nscale-1] = strdup(bp);
            khash_str2int_set(args->csq2severity,args->scale[args->nscale-1], severity);
            if ( i==args->nscale )
                fprintf(stderr,"Note: assigning a (high) severity score to a new consequence, use -S to override: %s -> %d\n",args->scale[args->nscale-1],args->nscale);

            if ( khash_str2int_get(args->csq2severity, bp, &severity)!=0 ) error("FIXME: failed to look up the consequence \"%s\"\n", bp);
        }
        if ( exact_match < 0 )
        {
            if ( *min_severity > severity ) *min_severity = severity;
            if ( *max_severity < severity ) *max_severity = severity;
        }
        else
        {
            if ( severity==exact_match )
            {
                *min_severity = *max_severity = severity;
                *ep = tmp;
                return;
            }
        }

        if ( !tmp ) break;
        *ep = tmp;
        ep++;
    }
}

static int csq_severity_pass(args_t *args, char *csq)
{
    if ( args->min_severity==args->max_severity && args->min_severity==SELECT_CSQ_ANY ) return 1;

    int min_severity, max_severity, exact_match = args->min_severity==args->max_severity ? args->min_severity : -1;
    csq_to_severity(args, csq, &min_severity, &max_severity, exact_match);
    if ( max_severity < args->min_severity ) return 0;
    if ( min_severity > args->max_severity ) return 0;
    return 1;
}

static int get_matching_transcript(args_t *args, bcf1_t *rec, int *list)
{
    int i, n = 0;
    for (i=0; i<args->ncols_csq; i++)
    {
        cols_t *cols_csq = args->cols_csq[i];
        if ( args->tr_expr.idx >= cols_csq->n )
            error("Too few columns at %s:%"PRId64" .. %d (field %s) >= %d\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1,args->tr_expr.idx,args->tr_expr.field,cols_csq->n);
        if ( args->tr_expr.op==TR_OP_EQ )
        {
            if ( !strcmp(args->tr_expr.value,cols_csq->off[args->tr_expr.idx]) ) list[n++] = i;
        }
        else if ( args->tr_expr.op==TR_OP_NE )
        {
            if ( strcmp(args->tr_expr.value,cols_csq->off[args->tr_expr.idx]) ) list[n++] = i;
        }
        else if ( args->tr_expr.op==TR_OP_RE )
        {
            int match = regexec(args->tr_expr.regex,cols_csq->off[args->tr_expr.idx], 0,NULL,0) ? 0 : 1;
            if ( match ) list[n++] = i;
        }
        else if ( args->tr_expr.op==TR_OP_NR )
        {
            int match = regexec(args->tr_expr.regex,cols_csq->off[args->tr_expr.idx], 0,NULL,0) ? 0 : 1;
            if ( !match ) list[n++] = i;
        }
    }
    return n;
}
static int get_worst_transcript(args_t *args, bcf1_t *rec)
{
    int i, max_severity = -1, imax_severity = 0;
    for (i=0; i<args->ncols_csq; i++)
    {
        cols_t *cols_csq = args->cols_csq[i];
        if ( args->csq_idx >= cols_csq->n )
            error("Too few columns at %s:%"PRId64" .. %d (Consequence) >= %d\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1,args->csq_idx,cols_csq->n);
        char *csq = cols_csq->off[args->csq_idx];

        int min, max;
        csq_to_severity(args, csq, &min, &max, -1);
        if ( max_severity < max ) { imax_severity = i; max_severity = max; }
    }
    return imax_severity;
}
static void annot_reset(annot_t *annot, int nannot)
{
    int i;
    for (i=0; i<nannot; i++) annot[i].str.l = 0;
}
static void annot_append(annot_t *ann, char *value)
{
    if ( ann->str.l ) kputc(',',&ann->str);
    kputs(value, &ann->str);
}
static inline void parse_array_real(annot_t *ann, char *str, float **arr, int *marr, int *narr)
{
    static int warned_type_err = 0;
    char *bp = str, *ep;
    float *ptr = *arr;
    int i, n = 1, m = *marr;
    for (i=0; *bp; bp++)
        if ( *bp == ',' ) n++;

    hts_expand(float*,n,m,ptr);

    i = 0;
    bp = str;
    while ( *bp )
    {
        ptr[i] = strtod(bp, &ep);
        if ( bp==ep )
        {
            if ( !warned_type_err && (ep[0]!='.' || (ep[1]!=',' && ep[1])) )
            {
                fprintf(stderr,"Warning: Could not parse, not a numeric list %s=\"%s\", check the -c and --columns-types options.\n"
                        "         This message is printed only once.\n",ann->tag,str);
                warned_type_err = 1;
            }
            bcf_float_set_missing(ptr[i]);
        }
        i++;
        while ( *ep && *ep!=',' ) ep++;
        bp = *ep ? ep + 1 : ep;
    }
    *narr = i;
    *marr = m;
    *arr  = ptr;
}
static inline void parse_array_int32(annot_t *ann, char *str, int **arr, int *marr, int *narr)
{
    static int warned_type_err = 0;
    char *bp = str, *ep;
    int32_t *ptr = *arr;
    int i, n = 1, m = *marr;
    for (i=0; *bp; bp++)
        if ( *bp == ',' ) n++;

    hts_expand(int32_t*,n,m,ptr);

    i = 0;
    bp = str;
    while ( *bp )
    {
        ptr[i] = strtol(bp, &ep, 10);
        if ( bp==ep )
        {
            if ( !warned_type_err && (ep[0]!='.' || (ep[1]!=',' && ep[1])) )
            {
                fprintf(stderr,"Warning: Could not parse, not a numeric list %s=\"%s\", check the -c and --columns-types options.\n"
                        "         This message is printed only once.\n",ann->tag,str);
                warned_type_err = 1;
            }
            ptr[i] = bcf_int32_missing;
        }
        i++;
        while ( *ep && *ep!=',' ) ep++;
        bp = *ep ? ep + 1 : ep;
    }
    *narr = i;
    *marr = m;
    *arr  = ptr;
}
static void filter_and_output(args_t *args, bcf1_t *rec, int severity_pass, int all_missing)
{
    int i, updated = 0;
    for (i=0; i<args->nannot; i++)
    {
        annot_t *ann = &args->annot[i];
        if ( !ann->str.l ) continue;
        if ( ann->type==BCF_HT_REAL )
        {
            parse_array_real(ann,ann->str.s,&args->farr,&args->mfarr,&args->nfarr);
            bcf_update_info_float(args->hdr_out,rec,ann->tag,args->farr,args->nfarr);
        }
        else if ( ann->type==BCF_HT_INT )
        {
            parse_array_int32(ann,ann->str.s,&args->iarr,&args->miarr,&args->niarr);
            bcf_update_info_int32(args->hdr_out,rec,ann->tag,args->iarr,args->niarr);
        }
        else
            bcf_update_info_string(args->hdr_out,rec,ann->tag,ann->str.s);
        updated++;
    }
    if ( args->filter )
    {
        int pass = filter_test(args->filter, rec, (const uint8_t**) &args->smpl_pass);
        if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
        if ( !pass ) return;
    }
    if ( args->format_str )
    {
        if ( args->nannot )
        {
            if ( args->drop_sites && (!updated || all_missing) ) return;         // the standard case: using -f to print the CSQ subfields, skipping if missing
        }
        else
        {
            if ( !severity_pass ) return;   // request to print only non-CSQ tags at sites that pass severity
        }

        args->kstr.l = 0;
        convert_line(args->convert, rec, &args->kstr);
        if ( args->kstr.l && bgzf_write(args->fh_bgzf, args->kstr.s, args->kstr.l)!=args->kstr.l )
            error("Failed to write to %s\n", args->output_fname);
        return;
    }
    if ( bcf_write(args->fh_vcf, args->hdr_out,rec)!=0 )
        error("Failed to write to %s\n", args->output_fname);
}
static void restrict_csqs_to_genes(args_t *args)
{
    int i,j, nhit = 0;

    // use iarr to mark transcripts from interesting genes
    hts_expand(int32_t,args->cols_tr->n,args->miarr,args->iarr);
    for (i=0; i<args->cols_tr->n; i++)
    {
        cols_t *cols_csq = args->cols_csq[i];
        for (j=0; j<args->ngene_fields; j++)
        {
            int idx = args->gene_fields[j];
            if ( idx >= cols_csq->n ) continue;
            if ( khash_str2int_has_key(args->genes,cols_csq->off[idx]) ) break;
        }
        if ( j < args->ngene_fields )
            args->iarr[i] = 1, nhit++;
        else
            args->iarr[i] = 0;
    }
    if ( !nhit )
    {
        if ( args->genes_mode==GENES_RESTRICT ) args->ncols_csq = 0;    // no gene of interest
        return;
    }

    // remove all uninteresting genes by putting all cols_csq hits first
    i = 0, j = args->cols_tr->n - 1;
    while ( i < j )
    {
        if ( args->iarr[i]==1 ) { i++; continue; }
        if ( args->iarr[j]==0 ) { j--; continue; }
        cols_t *tmp = args->cols_csq[i];
        args->cols_csq[i] = args->cols_csq[j];
        args->cols_csq[j] = tmp;
        char *tmp2 = args->cols_tr->off[i];
        args->cols_tr->off[i] = args->cols_tr->off[j];
        args->cols_tr->off[j] = tmp2;
        i++;
        j--;
    }
    args->ncols_csq = nhit;
}

// Split the VEP annotation by transcript and by field, then check if the number of subfields looks alright.
// Unfortunately, we cannot enforce the number of subfields to match the header definition because that can
// be variable: `bcftools csq` outputs different number of fields for different consequence types.
// So we need to distinguish between this reasonable case and incorrectly formatted consequences such
// as those reported for LoF_info subfield here https://github.com/Ensembl/ensembl-vep/issues/1351.
static void split_csq_fields(args_t *args, bcf1_t *rec, int csq_str_len)
{
    int i;
    args->cols_tr = cols_split(args->csq_str, args->cols_tr, ',');
    if ( args->cols_tr->n > args->mcols_csq )
    {
        args->cols_csq = (cols_t**)realloc(args->cols_csq,args->cols_tr->n*sizeof(*args->cols_csq));
        for (i=args->mcols_csq; i<args->cols_tr->n; i++) args->cols_csq[i] = NULL;
        args->mcols_csq = args->cols_tr->n;
    }
    args->ncols_csq = args->cols_tr->n;
    int nfield_diff = 0, need_fix = 0;
    for (i=0; i<args->cols_tr->n; i++)
    {
        args->cols_csq[i] = cols_split(args->cols_tr->off[i], args->cols_csq[i], '|');
        if ( args->csq_idx >= args->cols_csq[i]->n ) need_fix = 1;
        if ( nfield_diff < abs(args->cols_csq[i]->n - args->nfield) ) nfield_diff = args->cols_csq[i]->n - args->nfield;
    }
    if ( !csq_str_len ) return;     // called 2nd time, don't attempt to fix
    if ( !need_fix ) return;

    static int warned = 0;
    if ( !warned )
        fprintf(stderr,"Warning: The number of INFO/%s subfields at %s:%"PRIhts_pos" does not match the header definition,\n"
                       "         expected %d subfields, found as %s as %d. (This warning is printed only once.)\n",
                       args->vep_tag,bcf_seqname(args->hdr,rec),rec->pos+1,args->nfield,nfield_diff>0?"much":"few",args->nfield+nfield_diff);
    warned = 1;

    // One known failure mode is LoF_info subfield which can contain commas. Work around this by relying on
    // the number of pipe delimiters matching the header definition, replacing offending commas with slash
    // characters. This assumes that there is never a comma in the first subfield, otherwise it would be
    // impossible to distinguish between A|A,A,B,B|B and A|A,A,A,B|B
    int npipe = 0;
    while ( csq_str_len > 0 )
    {
        csq_str_len--;
        if ( args->csq_str[csq_str_len]=='|' ) { npipe++; continue; }
        if ( args->csq_str[csq_str_len]!=',' ) continue;
        if ( npipe && !(npipe % (args->nfield-1)) ) { npipe = 0; continue; }    // wholesome number of pipes encountered, this is a valid comma
        args->csq_str[csq_str_len] = '/';
    }
    split_csq_fields(args,rec,0);   // try once more
}
static void process_record(args_t *args, bcf1_t *rec)
{
    int i,len = bcf_get_info_string(args->hdr,rec,args->vep_tag,&args->csq_str,&args->ncsq_str);
    if ( len<=0 )
    {
        if ( !args->drop_sites )
        {
            annot_reset(args->annot, args->nannot);
            filter_and_output(args,rec,1,1);
        }
        return;
    }

    split_csq_fields(args,rec,len);

    // restrict to -g genes
    if ( args->genes ) restrict_csqs_to_genes(args);

    // select transcripts of interest
    int nlist_tr = 0;
    hts_expand(int32_t,args->ncols_csq,args->mlist_tr,args->list_tr);
    if ( !args->ncols_csq )
        nlist_tr = 0;           // no transcripts left after -g applied
    else if ( args->select_tr==SELECT_TR_EXPR )
    {
        nlist_tr = get_matching_transcript(args,rec,args->list_tr);
    }
    else if ( args->select_tr==SELECT_TR_WORST )
    {
        nlist_tr = 1;
        args->list_tr[0] = get_worst_transcript(args,rec);
    }
    else
    {
        for (i=0; i<args->ncols_csq; i++) args->list_tr[i] = i;
        nlist_tr = args->ncols_csq;
    }
    annot_reset(args->annot, args->nannot);
    int severity_pass = 0;  // consequence severity requested via the -s option (BCF record may be output but not annotated)
    int all_missing   = 1;  // transcripts with all requested annotations missing will be discarded if -f was given
    static int too_few_fields_warned = 0;
    int j,ilist;
    for (ilist=0; ilist<nlist_tr; ilist++)
    {
        i = args->list_tr[ilist];
        cols_t *cols_csq = args->cols_csq[i];
        if ( args->csq_idx >= cols_csq->n )
            error("Too few columns at %s:%"PRId64" .. %d (Consequence) >= %d\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1,args->csq_idx,cols_csq->n);

        char *csq = cols_csq->off[args->csq_idx];
        if ( !csq_severity_pass(args, csq) ) continue;
        severity_pass = 1;

        for (j=0; j<args->nannot; j++)
        {
            annot_t *ann = &args->annot[j];
            if ( ann->idx >= cols_csq->n )
            {
                if ( !too_few_fields_warned )
                {
                    fprintf(stderr, "Warning: fewer %s fields than expected at %s:%"PRId64", filling with dots. This warning is printed only once.\n", args->vep_tag,bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
                    too_few_fields_warned = 1;
                }
                annot_append(ann, ".");
                continue;
            }

            char *ann_str = NULL;
            if ( ann->idx==-1 ) ann_str = args->cols_tr->off[i];
            else if ( *cols_csq->off[ann->idx] ) ann_str = cols_csq->off[ann->idx];
            if ( ann_str )
            {
                annot_append(ann, ann_str);
                all_missing = 0;
            }
            else
                annot_append(ann, "."); // missing value
        }

        if ( args->duplicate )
        {
            filter_and_output(args, rec, severity_pass, all_missing);
            annot_reset(args->annot, args->nannot);
            all_missing   = 1;
            severity_pass = 0;
        }
    }
    if ( !severity_pass && args->drop_sites ) return;
    if ( !args->duplicate )
        filter_and_output(args, rec, severity_pass, all_missing);
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_fname = "-";
    args->output_type  = FT_VCF;
    args->record_cmd_line = 1;
    args->regions_overlap = 1;
    args->targets_overlap = 0;
    args->clevel = -1;
    static struct option loptions[] =
    {
        {"drop-sites",no_argument,0,'x'},
        {"keep-sites",no_argument,0,'X'},
        {"all-fields",no_argument,0,'A'},
        {"duplicate",no_argument,0,'d'},
        {"format",required_argument,0,'f'},
        {"gene-list",required_argument,0,'g'},
        {"gene-list-fields",required_argument,0,5},
        {"annotation",required_argument,0,'a'},
        {"annot-prefix",required_argument,0,'p'},
        {"columns",required_argument,0,'c'},
        {"columns-types",required_argument,0,1},
        {"select",required_argument,0,'s'},
        {"severity",required_argument,0,'S'},
        {"list",no_argument,0,'l'},
        {"include",required_argument,0,'i'},
        {"exclude",required_argument,0,'e'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"regions-overlap",required_argument,NULL,3},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"targets-overlap",required_argument,NULL,4},
        {"no-version",no_argument,NULL,2},
        {"allow-undef-tags",no_argument,0,'u'},
        {"write-index",optional_argument,NULL,'W'},
        {NULL,0,NULL,0}
    };
    int c, drop_sites = -1;
    char *tmp;
    while ((c = getopt_long(argc, argv, "o:O:i:e:r:R:t:T:lS:s:c:p:a:f:dA:xXuHg:W::",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case  2 : args->record_cmd_line = 0; break;
            case  1 : args->column_types = optarg; break;
            case 'A':
                if ( !strcasecmp(optarg,"tab") ) args->all_fields_delim = "\t";
                else if ( !strcasecmp(optarg,"space") ) args->all_fields_delim = " ";
                else args->all_fields_delim = optarg;
                break;
            case 'H': args->print_header++; break;
            case 'x': drop_sites = 1; break;
            case 'X': drop_sites = 0; break;
            case 'd': args->duplicate = 1; break;
            case 'f': args->format_str = strdup(optarg); break;
            case 'g': args->genes_fname = optarg; break;
            case 'a': args->vep_tag = optarg; break;
            case 'p': args->annot_prefix = optarg; break;
            case 'c': args->column_str = strdup(optarg); break;
            case 'S': args->severity = optarg; break;
            case 's': args->select = optarg; break;
            case 'l': args->list_hdr = 1; break;
            case 'u': args->allow_undef_tags = 1; break;
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
                      }
                      if ( optarg[1] )
                      {
                          args->clevel = strtol(optarg+1,&tmp,10);
                          if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                      }
                      break;
            case  3 :
                args->regions_overlap = parse_overlap_option(optarg);
                if ( args->regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case  4 :
                args->targets_overlap = parse_overlap_option(optarg);
                if ( args->targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case  5 : args->gene_fields_str = optarg; break;
            case 'W':
                if (!(args->write_index = write_index_parse(optarg)))
                    error("Unsupported index format '%s'\n", optarg);
                break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( args->print_header && !args->format_str ) error("Error: the -H header printing is supported only with -f\n");
    if ( args->all_fields_delim && !args->format_str ) error("Error: the -A option must be used with -f\n");
    if ( args->severity && (!strcmp("?",args->severity) || !strcmp("-",args->severity)) ) error("%s", default_severity());
    if ( args->column_types && !strcmp("-",args->column_types) ) error("%s", default_column_types());
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error("%s", usage_text()); }
    }
    else if ( optind+1!=argc ) error("%s", usage_text());
    else args->fname = argv[optind];

    init_data(args);

    if ( args->list_hdr )
        list_header(args);
    else
    {
        if ( !args->format_str && !args->column_str )
        {
            if ( args->min_severity==SELECT_CSQ_ANY && args->max_severity==SELECT_CSQ_ANY )
            {
                if ( args->select_tr!=SELECT_TR_EXPR )
                    error("Error: none of the -c,-f,-s options was given, why not use \"bcftools view\" instead?\n");
                if ( drop_sites==-1 )
                {
                    // Note that this is just to compensate for a common user confusion and should not really be here:
                    // the -s option is supposed to select from multiple transcripts, not to decide whether to keep or
                    // drop a site. However, it is tempting to run `bcftools +split-vep -s mane` and expect on output
                    // a VCF with sites impating MANE transcripts only.
                    drop_sites = 1;
                }
                else if ( !drop_sites )
                    error("Error: the option -X has no effect without -c,-f, why not use \"bcftools view\" instead?\n");
            }
        }
        if ( drop_sites==-1 ) drop_sites = args->format_str ? 1 : 0;
        args->drop_sites = drop_sites;

        if ( args->format_str )
        {
            args->fh_bgzf = bgzf_open(args->output_fname, args->output_type&FT_GZ ? "wg" : "wu");
            if ( args->print_header )
            {
                args->kstr.l = 0;
                convert_header(args->convert,&args->kstr);
                if ( args->kstr.l && bgzf_write(args->fh_bgzf, args->kstr.s, args->kstr.l)!=args->kstr.l )
                    error("Failed to write to %s\n", args->output_fname);
            }
        }
        else
        {
            char wmode[8];
            set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
            args->fh_vcf = hts_open(args->output_fname ? args->output_fname : "-", wmode);
            if ( args->record_cmd_line ) bcf_hdr_append_version(args->hdr_out, args->argc, args->argv, "bcftools_split-vep");
            if ( bcf_hdr_write(args->fh_vcf, args->hdr_out)!=0 ) error("Failed to write the header to %s\n", args->output_fname);
            if ( init_index2(args->fh_vcf,args->hdr,args->output_fname,
                             &args->index_fn, args->write_index)<0 )
                error("Error: failed to initialise index for %s\n",args->output_fname);
        }
        while ( bcf_sr_next_line(args->sr) )
            process_record(args, bcf_sr_get_line(args->sr,0));
    }

    destroy_data(args);

    return 0;
}
