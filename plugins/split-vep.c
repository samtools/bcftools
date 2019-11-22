/* The MIT License

   Copyright (c) 2019 Genome Research Ltd.

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
#include "../bcftools.h"
#include "../filter.h"
#include "../convert.h"
#include "../cols.h"


// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define SELECT_TR_ALL       0
#define SELECT_TR_WORST     1
#define SELECT_TR_PRIMARY   2
#define SELECT_CSQ_ANY      -1

typedef struct
{
    char *field;    // the name of the VEP field, e.g. Consequence,Gene,etc.
    char *tag;      // the name of the VCF tag: the annot_t.field with the -p prefix
    int idx;        // 0-based index within the VEP annotation string
    int type;       // annotation type, one of the BCF_HT_* types
    kstring_t str;  // annotation value, ready to pass to bcf_update_info_*
}
annot_t;

typedef struct
{
    convert_t *convert;
    filter_t *filter;
    int argc, filter_logic, regions_is_file, targets_is_file, list_hdr;
    kstring_t kstr;
    char *filter_str,
        *vep_tag;       // the --annotation INFO tag to process
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
    int csq_idx,        // the index of the Consequence field; for the --select CSQ option
        primary_id;     // the index of the CANONICAL field; for the --select TR option
    char *severity,     // the --severity scale option
        *select,        // the --select option
        *column_str,    // the --columns option
        *annot_prefix;  // the --annot-prefix option
    void *field2idx,    // VEP field name to index, used in initialization
        *csq2severity;  // consequence type to severity score
    cols_t *cols_tr,    // the current CSQ tag split into transcripts
        *cols_csq;      // the current CSQ transcript split into fields
    int min_severity, max_severity;     // ignore consequences outside this severity range
    int drop_sites;                     // the -x, --drop-sites option
    int select_tr;                      // one of SELECT_TR_*
    uint8_t *smpl_pass;                 // for filtering at sample level, used with -f
    int duplicate;              // the -d, --duplicate option is set
    char *all_fields_delim;     // the -A, --all-fields option is set
    float *farr;                // helper arrays for bcf_update_* functions
    int32_t *iarr;
    int niarr,miarr, nfarr,mfarr;
}
args_t;

args_t args;

const char *about(void)
{
    return "Query structured annotations such as the CSQ created by VEP.\n";
}

static const char *default_severity(void)
{
    return
        "# Default consequence substrings ordered in ascending order by severity.\n"
        "# Consequences with the same severity can be put on the same line in arbitrary order.\n"
        "intergenic\n"
        "downstream upstream\n"
        "intron\n"
        "non_coding\n"
        "regulatory\n"
        "5_prime_utr 3_prime_utr\n"
        "stop_retained start_retained synonymous\n"
        "splice_region\n"
        "coding_sequence\n"
        "missense\n"
        "inframe\n"
        "exon_loss\n"
        "disruptive\n"
        "splice_acceptor splice_donor\n"
        "start_lost stop_lost stop_gained frameshift\n";
}
static const char *usage_text(void)
{
    return 
        "\n"
        "About: Query structured annotations such INFO/CSQ created by bcftools/csq or VEP. For more\n"
        "   more information and pointers see http://samtools.github.io/bcftools/howtos/plugin.split-vep.html\n"
        "Usage: bcftools +split-vep [Plugin Options]\n"
        "Plugin options:\n"
        "   -a, --annotation STR        INFO annotation to parse [CSQ]\n"
        "   -A, --all-fields DELIM      Output all fields replacing the -a tag (\"%CSQ\" by default) in the -f\n"
        "                                 filtering expression using the output field delimiter DELIM. This can be\n"
        "                                 \"tab\", \"space\" or an arbitrary string.\n"
        "   -c, --columns LIST[:type]   Extract the fields listed either as indexes or names. The default type\n"
        "                                 of the new annotation is String but can be also Integer/Int or Float/Real.\n"
        "   -d, --duplicate             Output per transcript/allele consequences on a new line rather rather than\n"
        "                                 as comma-separated fields on a single line\n"
        "   -f, --format <string>       Formatting expression for non-VCF/BCF output, same as `bcftools query -f`\n"
        "   -l, --list                  Parse the VCF header and list the annotation fields\n"
        "   -p, --annot-prefix          Prefix of INFO annotations to be created after splitting the CSQ string\n"
        "   -s, --select TR:CSQ         Select transcripts to extract by type and/or consequence. (See also the -x switch.)\n"
        "                                 TR, transcript:   worst,primary(*),all        [all]\n"
        "                                 CSQ, consequence: any,missense,missense+,etc  [any]\n"
        "                                 (*) Primary transcripts have the field \"CANONICAL\" set to \"YES\"\n"
        "   -S, --severity -|FILE       Pass \"-\" to print the default severity scale or FILE to override\n"
        "                                 the default scale\n"
        "   -x, --drop-sites            Drop sites with none of the consequences matching the severity specified by -s.\n"
        "                                  This switch is intended for use with VCF/BCF output (i.e. -f not given).\n"
        "Common options:\n"
        "   -e, --exclude EXPR          Exclude sites and samples for which the expression is true\n"
        "   -i, --include EXPR          Include sites and samples for which the expression is true\n"
        "   -o, --output FILE           Output file name [stdout]\n"
        "   -O, --output-type b|u|z|v   b: compressed BCF, u: uncompressed BCF, z: compressed VCF or text, v: uncompressed VCF or text [v]\n"
        "   -r, --regions REG           Restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE     Restrict to regions listed in a file\n"
        "   -t, --targets REG           Similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE     Similar to -R but streams rather than index-jumps\n"
        "\n"
        "Examples:\n"
        "   # List available fields of the INFO/CSQ annotation\n"
        "   bcftools +split-vep -l file.vcf.gz\n"
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
        "   # or equivalent. In order to drop sites with different consequences completely, we add\n"
        "   # the -x switch. See the online documentation referenced above for more examples.\n"
        "   bcftools +split-vep -c gnomAD_AF:Float -s :missense    file.vcf.gz\n"
        "   bcftools +split-vep -c gnomAD_AF:Float -s :missense -x file.vcf.gz\n"
        "\n";
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

static void init_data(args_t *args)
{
    args->sr = bcf_sr_init();
    if ( args->regions )
    {
        args->sr->require_index = 1;
        if ( bcf_sr_set_regions(args->sr, args->regions, args->regions_is_file)<0 ) error("Failed to read the regions: %s\n",args->regions);
    }
    if ( args->targets && bcf_sr_set_targets(args->sr, args->targets, args->targets_is_file, 0)<0 ) error("Failed to read the targets: %s\n",args->targets);
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr = bcf_sr_get_header(args->sr,0);
    args->hdr_out = bcf_hdr_dup(args->hdr);

    // Parse the header CSQ line, must contain Description with "Format: ..." declaration
    bcf_hrec_t *hrec = bcf_hdr_get_hrec(args->hdr, BCF_HL_INFO, NULL, args->vep_tag, NULL);
    if ( !hrec ) error("The tag INFO/%s not found in the header\n", args->vep_tag);
    int ret = bcf_hrec_find_key(hrec, "Description");
    if ( ret<0 ) error("No \"Description\" field was found for the tag INFO/%s in the header\n", args->vep_tag);
    char *format = strstr(hrec->vals[ret], "Format: ");
    if ( !format ) error("Expected \"Format: \" substring in the header INFO/%s/Description, found: %s\n", args->vep_tag,hrec->vals[ret]);
    format += 8;
    char *ep = format;
    while ( *ep )
    {
        char *bp = ep;
        while ( *ep && *ep!='|' ) ep++;
        char tmp = *ep;
        *ep = 0;
        args->nfield++;
        args->field = (char**)realloc(args->field,args->nfield*sizeof(*args->field));
        args->field[args->nfield-1] = strdup(bp);
        if ( !tmp ) break;
        ep++;
    }
    if ( !args->nfield ) error("Could not parse Description of INFO/%s: %s\n", args->vep_tag,hrec->vals[ret]);
    int len = strlen(args->field[args->nfield-1]);
    if ( args->field[args->nfield-1][len-1]=='"' ) args->field[args->nfield-1][len-1] = 0;    // remove the trailing doublequote character
    args->field2idx = khash_str2int_init();
    int i,j;
    for (i=0; i<args->nfield; i++)
    {
        if ( khash_str2int_has_key(args->field2idx, args->field[i]) )
        {
            fprintf(stderr,"Warning: duplicate INFO/%s key \"%s\"\n", args->vep_tag,args->field[i]);
            continue;
        }
        khash_str2int_set(args->field2idx, args->field[i], i);
    }

    // Create a text output as with `bcftools query -f`. For this we need to determine the fields to be extracted
    // from the formatting expression
    kstring_t str = {0,0,0};
    if ( args->format_str && !args->column_str )
    {
        // Special case: -A was given, extract all fields, for this the -a tag (%CSQ) must be present
        if ( args->all_fields_delim ) expand_csq_expression(args, &str);

        for (i=0; i<args->nfield; i++)
        {
            str.l = 0;
            kputc('%',&str);
            kputs(args->field[i],&str);
            char end, *ptr = args->format_str;
            while ( ptr )
            {
                ptr = strstr(ptr,str.s);
                if ( !ptr ) break;
                end = ptr[str.l];
                if ( isalnum(end) || end=='_' || end=='.' )
                {
                    ptr++;
                    continue;
                }
                break;
            }
            if ( !ptr ) continue;
            ptr[str.l] = 0;
            int tag_id = bcf_hdr_id2int(args->hdr, BCF_DT_ID, ptr+1);
            if ( bcf_hdr_idinfo_exists(args->hdr,BCF_HL_INFO,tag_id) )
                fprintf(stderr,"Note: ambigous key %s, using the %s subfield of %s, not the INFO/%s tag\n", ptr,ptr+1,args->vep_tag,ptr+1);

            int olen = args->column_str ? strlen(args->column_str) : 0;
            int nlen = strlen(ptr) - 1;
            args->column_str = (char*)realloc(args->column_str, olen + nlen + 2);
            if ( olen )
            {
                memcpy(args->column_str+olen,",",1);
                olen++;
            }
            memcpy(args->column_str+olen,ptr+1,nlen);
            args->column_str[olen+nlen] = 0;

            ptr[str.l] = end;
        }
    }

    // The "Consequence" column to look up severity, its name is hardwired for now
    if ( khash_str2int_get(args->field2idx,"Consequence",&args->csq_idx)!=0 )
        error("The field \"Consequence\" is not present in INFO/%s: %s\n", args->vep_tag,hrec->vals[ret]);

    // Columns to extract: given as names, 0-based indexes or ranges of indexes
    if ( args->column_str )
    {
        int *column = NULL;
        int *types  = NULL;
        ep = args->column_str;
        while ( *ep )
        {
            char *tp, *bp = ep;
            while ( *ep && *ep!=',' ) ep++;
            char tmp = *ep;
            *ep = 0;
            int type = BCF_HT_STR;
            int idx_beg, idx_end;
            if ( khash_str2int_get(args->field2idx, bp, &idx_beg)==0 )
                idx_end = idx_beg;
            else if ( (tp=strrchr(bp,':')) )
            {
                *tp = 0;
                if ( khash_str2int_get(args->field2idx, bp, &idx_beg)!=0 )
                {
                    *tp = ':';
                    error("No such column: \"%s\"\n", bp);
                }
                idx_end = idx_beg;
                *tp = ':';
                if ( !strcasecmp(tp+1,"string") ) type = BCF_HT_STR;
                else if ( !strcasecmp(tp+1,"float") || !strcasecmp(tp+1,"real") ) type = BCF_HT_REAL;
                else if ( !strcasecmp(tp+1,"integer") || !strcasecmp(tp+1,"int") ) type = BCF_HT_INT;
                else if ( !strcasecmp(tp+1,"flag") ) type = BCF_HT_FLAG;
                else error("The type \"%s\" (or column \"%s\"?) not recognised\n", tp+1,bp);
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
                        if ( !strcasecmp(mp+1,"string") ) type = BCF_HT_STR;
                        else if ( !strcasecmp(mp+1,"float") || !strcasecmp(mp+1,"real") ) type = BCF_HT_REAL;
                        else if ( !strcasecmp(mp+1,"integer") || !strcasecmp(mp+1,"int") ) type = BCF_HT_INT;
                        else if ( !strcasecmp(mp+1,"flag") ) type = BCF_HT_FLAG;
                        else error("The type \"%s\" (or column \"%s\"?) not recognised\n", mp+1,bp);
                    }
                    else 
                        error("No such column: \"%s\"\n", bp);
                }
            }

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
            if ( !tmp ) break;
            ep++;
        }
        args->annot = (annot_t*)calloc(args->nannot,sizeof(*args->annot));
        int len = args->annot_prefix ? strlen(args->annot_prefix) : 0;
        for (i=0; i<args->nannot; i++)
        {
            annot_t *ann = &args->annot[i];
            ann->type = types[i];
            ann->idx = j = column[i];
            ann->field = strdup(args->field[j]);
            int clen = strlen(args->field[j]);
            ann->tag = (char*)malloc(clen+len+1);
            if ( len ) memcpy(ann->tag,args->annot_prefix,len);
            memcpy(ann->tag+len,ann->field,clen);
            ann->tag[len+clen] = 0;
            args->kstr.l = 0;
            char *type = "String";
            if ( ann->type==BCF_HT_REAL ) type = "Float";
            else if ( ann->type==BCF_HT_INT ) type = "Integer";
            else if ( ann->type==BCF_HT_FLAG ) type = "Flag";
            ksprintf(&args->kstr,"##INFO=<ID=%%s,Number=.,Type=%s,Description=\"The %%s field from INFO/%%s\">",type);
            bcf_hdr_printf(args->hdr_out, args->kstr.s, ann->tag,ann->field,args->vep_tag);
        }
        free(column);
        free(types);

        if ( bcf_hdr_sync(args->hdr_out)<0 )
            error_errno("[%s] Failed to update header", __func__);
    }
    if ( args->format_str )
    {
        if ( !args->column_str && !args->select ) error("Error: No %s field selected in the formatting expression and -s not given: a typo?\n",args->vep_tag);
        args->convert = convert_init(args->hdr_out, NULL, 0, args->format_str);
        if ( !args->convert ) error("Could not parse the expression: %s\n", args->format_str);
    }
    if ( args->filter_str )
    {
        int max_unpack = args->convert ? convert_max_unpack(args->convert) : 0;
        args->filter = filter_init(args->hdr_out, args->filter_str);
        max_unpack |= filter_max_unpack(args->filter);
        args->sr->max_unpack = max_unpack;
        if ( max_unpack & BCF_UN_FMT )
            convert_set_option(args->convert, subset_samples, &args->smpl_pass);
    }

    // Severity scale
    args->csq2severity = khash_str2int_init();
    int severity  = 0;
    str.l = 0;
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
    free(str.s);

    // Transcript and/or consequence selection
    if ( !args->select ) args->select = "all:any";
    cols_t *cols = cols_split(args->select, NULL, ':');
    char *sel_tr  = cols->off[0][0] ? cols->off[0] : "all";
    char *sel_csq = cols->n==2 && cols->off[1][0] ? cols->off[1] : "any";
    if ( !strcasecmp(sel_tr,"all") ) args->select_tr = SELECT_TR_ALL;
    else if ( !strcasecmp(sel_tr,"worst") ) args->select_tr = SELECT_TR_WORST;
    else if ( !strcasecmp(sel_tr,"primary") ) args->select_tr = SELECT_TR_PRIMARY;
    else error("Error: the transcript selection key \"%s\" is not recognised.\n", sel_tr);
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

    // The 'CANONICAL' column to look up severity, its name is hardwired for now
    if ( args->select_tr==SELECT_TR_PRIMARY && khash_str2int_get(args->field2idx,"CANONICAL",&args->primary_id)!=0 )
        error("The primary transcript was requested but the field \"CANONICAL\" is not present in INFO/%s: %s\n",args->vep_tag,hrec->vals[ret]);
}
static void destroy_data(args_t *args)
{
    free(args->farr);
    free(args->iarr);
    free(args->kstr.s);
    free(args->column_str);
    free(args->format_str);
    cols_destroy(args->cols_csq);
    cols_destroy(args->cols_tr);
    int i;
    for (i=0; i<args->nscale; i++) free(args->scale[i]);
    free(args->scale);
    for (i=0; i<args->nfield; i++) free(args->field[i]);
    free(args->field);
    for (i=0; i<args->nannot; i++)
    {
        annot_t *ann = &args->annot[i];
        free(ann->field);
        free(ann->tag);
        free(ann->str.s);
    }
    free(args->annot);
    if ( args->field2idx ) khash_str2int_destroy(args->field2idx);
    if ( args->csq2severity ) khash_str2int_destroy(args->csq2severity);
    bcf_sr_destroy(args->sr);
    bcf_hdr_destroy(args->hdr_out);
    free(args->csq_str);
    if ( args->filter ) filter_destroy(args->filter);
    if ( args->convert ) convert_destroy(args->convert);
    if ( args->fh_vcf && hts_close(args->fh_vcf)!=0 ) error("Error: close failed .. %s\n",args->output_fname);
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

static int get_primary_transcript(args_t *args, bcf1_t *rec, cols_t *cols_tr)    // modifies args->cols_csq!
{
    int i;
    for (i=0; i<cols_tr->n; i++)
    {
        args->cols_csq = cols_split(cols_tr->off[i], args->cols_csq, '|');
        if ( args->primary_id >= args->cols_csq->n )
            error("Too few columns at %s:%"PRId64" .. %d (Consequence) >= %d\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1,args->primary_id,args->cols_csq->n);
        if ( !strcmp("YES",args->cols_csq->off[args->primary_id]) ) return i;
    }
    return -1;
}
static int get_worst_transcript(args_t *args, bcf1_t *rec, cols_t *cols_tr)     // modifies args->cols_csq!
{
    int i, max_severity = -1, imax_severity = 0;
    for (i=0; i<cols_tr->n; i++)
    {
        args->cols_csq = cols_split(cols_tr->off[i], args->cols_csq, '|');
        if ( args->csq_idx >= args->cols_csq->n )
            error("Too few columns at %s:%"PRId64" .. %d (Consequence) >= %d\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1,args->csq_idx,args->cols_csq->n);
        char *csq = args->cols_csq->off[args->csq_idx];

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
static inline void parse_array_real(char *str, float **arr, int *marr, int *narr)
{
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
            bcf_float_set_missing(ptr[i]);
        i++;
        while ( *ep && *ep!=',' ) ep++;
        bp = *ep ? ep + 1 : ep;
    }
    *narr = i;
    *marr = m;
    *arr  = ptr;
}
static inline void parse_array_int32(char *str, int **arr, int *marr, int *narr)
{
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
            ptr[i] = bcf_int32_missing;
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
            parse_array_real(ann->str.s,&args->farr,&args->mfarr,&args->nfarr);
            bcf_update_info_float(args->hdr_out,rec,ann->tag,args->farr,args->nfarr);
        }
        else if ( ann->type==BCF_HT_INT )
        {
            parse_array_int32(ann->str.s,&args->iarr,&args->miarr,&args->niarr);
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
            if ( !updated || all_missing ) return;         // the standard case: using -f to print the CSQ subfields, skipping if missing
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
static void process_record(args_t *args, bcf1_t *rec)
{
    int len = bcf_get_info_string(args->hdr,rec,args->vep_tag,&args->csq_str,&args->ncsq_str);
    if ( len<=0 ) return;

    args->cols_tr = cols_split(args->csq_str, args->cols_tr, ',');

    int i,j, itr_min = 0, itr_max = args->cols_tr->n - 1;
    if ( args->select_tr==SELECT_TR_PRIMARY )
    {
        itr_min = itr_max = get_primary_transcript(args, rec, args->cols_tr);
        if ( itr_min<0 ) itr_max = itr_min - 1;
    }
    else if ( args->select_tr==SELECT_TR_WORST )
        itr_min = itr_max = get_worst_transcript(args, rec, args->cols_tr);

    annot_reset(args->annot, args->nannot);
    int severity_pass = 0;  // consequence severity requested via the -s option (BCF record may be output but not annotated)
    int all_missing   = 1;  // transcripts with all requested annotations missing will be discarded if -f was given
    static int too_few_fields_warned = 0;
    for (i=itr_min; i<=itr_max; i++)
    {
        args->cols_csq = cols_split(args->cols_tr->off[i], args->cols_csq, '|');
        if ( args->csq_idx >= args->cols_csq->n )
            error("Too few columns at %s:%"PRId64" .. %d (Consequence) >= %d\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1,args->csq_idx,args->cols_csq->n);

        char *csq = args->cols_csq->off[args->csq_idx];
        if ( !csq_severity_pass(args, csq) ) continue;
        severity_pass = 1;

        for (j=0; j<args->nannot; j++)
        {
            annot_t *ann = &args->annot[j];
            if ( ann->idx >= args->cols_csq->n )
            {
                if ( !too_few_fields_warned )
                {
                    fprintf(stderr, "Warning: fewer %s fields than expected at %s:%"PRId64", filling with dots. This warning is printed only once.\n", args->vep_tag,bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
                    too_few_fields_warned = 1;
                }
                annot_append(ann, ".");
                continue;
            }

            if ( !*args->cols_csq->off[ann->idx] )
                annot_append(ann, "."); // missing value
            else
            {
                annot_append(ann, args->cols_csq->off[ann->idx]);
                all_missing = 0;
            }
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
    args->vep_tag = "CSQ";
    static struct option loptions[] =
    {
        {"drop-sites",no_argument,0,'x'},
        {"all-fields",no_argument,0,'A'},
        {"duplicate",no_argument,0,'d'},
        {"format",required_argument,0,'f'},
        {"annotation",required_argument,0,'a'},
        {"annot-prefix",required_argument,0,'p'},
        {"columns",required_argument,0,'c'},
        {"select",required_argument,0,'s'},
        {"severity",required_argument,0,'S'},
        {"list",no_argument,0,'l'},
        {"include",required_argument,0,'i'},
        {"exclude",required_argument,0,'e'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {NULL,0,NULL,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "o:O:i:e:r:R:t:T:lS:s:c:p:a:f:dA:x",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'A':
                if ( !strcasecmp(optarg,"tab") ) args->all_fields_delim = "\t";
                else if ( !strcasecmp(optarg,"space") ) args->all_fields_delim = " ";
                else args->all_fields_delim = optarg;
                break;
            case 'x': args->drop_sites = 1; break;
            case 'd': args->duplicate = 1; break;
            case 'f': args->format_str = strdup(optarg); break;
            case 'a': args->vep_tag = optarg; break;
            case 'p': args->annot_prefix = optarg; break;
            case 'c': args->column_str = strdup(optarg); break;
            case 'S': args->severity = optarg; break;
            case 's': args->select = optarg; break;
            case 'l': args->list_hdr = 1; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
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
                          default: error("The output type \"%s\" not recognised\n", optarg);
                      }
                      break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( args->drop_sites && args->format_str ) error("Error: the -x behavior is the default (and only supported) with -f\n");
    if ( args->all_fields_delim && !args->format_str ) error("Error: the -A option must be used with -f\n");
    if ( args->severity && (!strcmp("?",args->severity) || !strcmp("-",args->severity)) ) error("%s", default_severity());
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
                error("Error: none of the -c,-f,-s options was given, why not use \"bcftools view\" instead?\n");
            else if ( !args->drop_sites )
                error("Error: when the -s option is used without -x, everything is printed; why not use \"bcftools view\" instead?\n");
        }

        if ( args->format_str )
            args->fh_bgzf = bgzf_open(args->output_fname, args->output_type&FT_GZ ? "wg" : "wu");
        else
        {
            args->fh_vcf = hts_open(args->output_fname, hts_bcf_wmode(args->output_type));
            if ( bcf_hdr_write(args->fh_vcf, args->hdr_out)!=0 ) error("Failed to write the header to %s\n", args->output_fname);
        }
        while ( bcf_sr_next_line(args->sr) )
            process_record(args, bcf_sr_get_line(args->sr,0));
    }

    destroy_data(args);

    return 0;
}
