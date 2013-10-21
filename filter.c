#include <ctype.h>
#include <stdlib.h>
#include <errno.h>
#include "filter.h"
#include "bcftools.h"

typedef struct _token_t
{
    // read-only values, same for all VCF lines
    int tok_type;       // one of the TOK_* keys below
    char *key;          // set only for string constants, otherwise NULL
    char *tag;          // for debugging and printout only, VCF tag name
    double threshold;   // filtering threshold
    int hdr_id;         // BCF header lookup ID
    int idx;            // 1-based index to VCF vector
    void (*setter)(bcf1_t *, struct _token_t *);
    int (*comparator)(struct _token_t *, struct _token_t *, int op_type, bcf1_t *);

    // modified on filter evaluation at each VCF line
    double num_value;   // in case str_value is set, num_value is the string's length
    char *str_value;
    int pass;           // -1 not applicable, 0 fails, >0 pass
    int missing_value;
}
token_t;

struct _filter_t 
{
    bcf_hdr_t *hdr;
    char *str;
    int nfilters;
    token_t *filters, **flt_stack;  // filtering input tokens (in RPN) and evaluation stack
};


#define TOK_VAL  0
#define TOK_LFT  1       // (
#define TOK_RGT  2       // )
#define TOK_LE   3       // less or equal
#define TOK_LT   4       // less than
#define TOK_EQ   5       // equal
#define TOK_BT   6       // bigger than
#define TOK_BE   7       // bigger or equal
#define TOK_NE   8       // not equal
#define TOK_OR   9       // |
#define TOK_AND  10      // &
#define TOK_ADD  11      // +
#define TOK_SUB  12      // -
#define TOK_MULT 13      // *
#define TOK_DIV  14      // /

//                        ( ) [ < = > ] ! | & + - * /
static int op_prec[] = {0,1,1,5,5,5,5,5,5,2,3,6,6,7,7};

static int filters_next_token(char **str, int *len)
{
    char *tmp = *str;
    while ( *tmp && isspace(*tmp) ) tmp++;
    *str = tmp;
    *len = 0;

    while ( tmp[0] )
    {
        if ( tmp[0]=='"' ) break;
        if ( isspace(tmp[0]) ) break;
        if ( tmp[0]=='<' ) break;
        if ( tmp[0]=='>' ) break;
        if ( tmp[0]=='=' ) break;
        if ( tmp[0]=='!' ) break;
        if ( tmp[0]=='&' ) break;
        if ( tmp[0]=='|' ) break;
        if ( tmp[0]=='(' ) break;
        if ( tmp[0]==')' ) break;
        if ( tmp[0]=='+' ) break;
        if ( tmp[0]=='*' ) break;
        if ( tmp[0]=='-' ) break;
        if ( tmp[0]=='/' ) break;
        tmp++;
    }
    if ( tmp > *str )
    {
        *len = tmp - (*str);
        return TOK_VAL;
    }
    if ( tmp[0]=='"' )
    {
        tmp++;
        while ( *tmp && tmp[0]!='"' ) tmp++;
        if ( !*tmp ) return -1;     // missing quotes
        *len = tmp - (*str) + 1;
        return TOK_VAL;
    }
    if ( tmp[0]=='!' )
    {
        if ( tmp[1]=='=' ) { (*str) += 2; return TOK_NE; }
    }
    if ( tmp[0]=='<' )
    {
        if ( tmp[1]=='=' ) { (*str) += 2; return TOK_LE; }
        (*str) += 1; return TOK_LT;
    }
    if ( tmp[0]=='>' )
    {
        if ( tmp[1]=='=' ) { (*str) += 2; return TOK_BE; }
        (*str) += 1; return TOK_BT;
    }
    if ( tmp[0]=='=' )
    {
        if ( tmp[1]=='=' ) { (*str) += 2; return TOK_EQ; }
        (*str) += 1; return TOK_EQ;
    }
    if ( tmp[0]=='(' ) { (*str) += 1; return TOK_LFT; }
    if ( tmp[0]==')' ) { (*str) += 1; return TOK_RGT; }
    if ( tmp[0]=='&' && tmp[1]=='&' ) { (*str) += 2; return TOK_AND; }
    if ( tmp[0]=='|' && tmp[1]=='|' ) { (*str) += 2; return TOK_OR; }
    if ( tmp[0]=='&' ) { (*str) += 1; return TOK_AND; }
    if ( tmp[0]=='|' ) { (*str) += 1; return TOK_OR; }
    if ( tmp[0]=='+' ) { (*str) += 1; return TOK_ADD; }
    if ( tmp[0]=='-' ) { (*str) += 1; return TOK_SUB; }
    if ( tmp[0]=='*' ) { (*str) += 1; return TOK_MULT; }
    if ( tmp[0]=='/' ) { (*str) += 1; return TOK_DIV; }

    while ( *tmp && !isspace(*tmp) )
    {
        if ( *tmp=='<' ) break;
        if ( *tmp=='>' ) break;
        if ( *tmp=='=' ) break;
        if ( *tmp=='&' ) break;
        if ( *tmp=='|' ) break;
        if ( *tmp=='(' ) break;
        if ( *tmp==')' ) break;
        if ( *tmp=='+' ) break;
        if ( *tmp=='-' ) break;
        if ( *tmp=='*' ) break;
        if ( *tmp=='/' ) break;
        tmp++;
    }
    *len = tmp - (*str);
    return TOK_VAL;
}

static void filters_set_qual(bcf1_t *line, token_t *tok)
{
    float *ptr = &line->qual;
    if ( bcf_float_is_missing(*ptr) )
        tok->missing_value = 1;
    else
        tok->num_value = line->qual;
}
static void filters_set_type(bcf1_t *line, token_t *tok)
{
    tok->num_value = bcf_get_variant_types(line);
}
static void filters_set_info(bcf1_t *line, token_t *tok)
{
    assert( tok->hdr_id >=0  );
    int i;
    for (i=0; i<line->n_info; i++)
        if ( line->d.info[i].key == tok->hdr_id ) break;

    if ( i==line->n_info ) 
        tok->missing_value = 1;
    else if ( line->d.info[i].type==BCF_BT_CHAR )
    {
        tok->str_value = (char*)line->d.info[i].vptr;       // string is typically not null-terminated
        tok->num_value = line->d.info[i].len;
    }
    else if ( line->d.info[i].type==BCF_BT_FLOAT )
    {
        tok->num_value = line->d.info[i].v1.f;
        tok->str_value = NULL;
    }
    else
    {
        tok->num_value = line->d.info[i].v1.i;
        tok->str_value = NULL;
    }
}
static void filters_set_filter(bcf1_t *line, token_t *tok) {}  // dummy
static int filters_cmp_filter(token_t *atok, token_t *btok, int op_type, bcf1_t *line)
{
    int i;
    if ( op_type==TOK_NE )  // AND logic: none of the filters can match
    {
        if ( !line->d.n_flt ) 
        {
            if ( atok->hdr_id==-1 ) return 0;   // missing value
            return 1; // no filter present, eval to true
        }
        for (i=0; i<line->d.n_flt; i++)
            if ( atok->hdr_id==line->d.flt[i] ) return 0;
        return 1;
    }
    // TOK_EQ with OR logic: at least one of the filters must match
    if ( !line->d.n_flt ) 
    {
        if ( atok->hdr_id==-1 ) return 1;
        return 0; // no filter present, eval to false
    }
    for (i=0; i<line->d.n_flt; i++)
        if ( atok->hdr_id==line->d.flt[i] ) return 1;
    return 0;
}

/**
 *  bcf_get_info_value() - get single INFO value, int or float
 *  @line:      BCF line
 *  @info_id:   tag ID, as returned by bcf_hdr_id2int
 *  @ivec:      0-based index to retrieve
 *  @vptr:      pointer to memory location of sufficient size to accomodate
 *              info_id's type
 *
 *  The returned value is -1 if tag is not present, 0 if present but
 *  values is missing or ivec is out of range, and 1 on success.
 */
int bcf_get_info_value(bcf1_t *line, int info_id, int ivec, void *value)
{
    int j;
    for (j=0; j<line->n_info; j++)
        if ( line->d.info[j].key == info_id ) break;
    if ( j==line->n_info ) return -1;

    bcf_info_t *info = &line->d.info[j];
    if ( info->len == 1 )
    {
        if ( info->type==BCF_BT_FLOAT ) *((float*)value) = info->v1.f;
        else if ( info->type==BCF_BT_INT8 || info->type==BCF_BT_INT16 || info->type==BCF_BT_INT32 ) *((int*)value) = info->v1.i;
        return 1;
    }

    #define BRANCH(type_t, is_missing, is_vector_end, out_type_t) { \
        type_t *p = (type_t *) info->vptr; \
        for (j=0; j<ivec && j<info->len; j++) \
        { \
            if ( is_vector_end ) return 0; \
        } \
        if ( is_missing ) return 0; \
        *((out_type_t*)value) = p[j]; \
        return 1; \
    }
    switch (info->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  p[j]==bcf_int8_missing,  p[j]==bcf_int8_vector_end,  int); break;
        case BCF_BT_INT16: BRANCH(int16_t, p[j]==bcf_int16_missing, p[j]==bcf_int16_vector_end, int); break;
        case BCF_BT_INT32: BRANCH(int32_t, p[j]==bcf_int32_missing, p[j]==bcf_int32_vector_end, int); break;
        case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(p[j]), bcf_float_is_vector_end(p[j]), float); break;
        default: fprintf(stderr,"todo: type %d\n", info->type); exit(1); break;
    }
    #undef BRANCH
    return -1;  // this shouldn't happen
}

static void filters_set_info_int(bcf1_t *line, token_t *tok)
{
    int value;
    if ( bcf_get_info_value(line,tok->hdr_id,tok->idx,&value) <= 0 )
        tok->missing_value = 1;
    else
        tok->num_value = value;
}

static void filters_set_info_float(bcf1_t *line, token_t *tok)
{
    float value;
    if ( bcf_get_info_value(line,tok->hdr_id,tok->idx,&value) <= 0 )
        tok->missing_value = 1;
    else
        tok->num_value = value;
}

static void filters_set_info_flag(bcf1_t *line, token_t *tok)
{
    int j;
    for (j=0; j<line->n_info; j++)
        if ( line->d.info[j].key == tok->hdr_id ) break;
    tok->num_value = j==line->n_info ? 0 : 1;
}

static int filters_init1(filter_t *filter, char *str, int len, token_t *tok)
{
    tok->tok_type = TOK_VAL;
    tok->hdr_id   = -1;
    tok->pass     = -1;

    // is this a string constant?
    if ( str[0]=='"' )
    {
        if ( str[len-1] != '"' ) error("TODO: [%s]\n", filter->str);
        tok->key = (char*) calloc(len-1,sizeof(char));
        tok->num_value = len-2;
        memcpy(tok->key,str+1,len-2);
        tok->key[len-2] = 0;
        return 0;
    }

    if ( !strncmp(str,"%QUAL",len) )
    {
        tok->setter = filters_set_qual;
        tok->tag = strdup("%QUAL");
        return 0;
    }
    if ( !strncmp(str,"%TYPE",len) )
    {
        tok->setter = filters_set_type;
        tok->tag = strdup("%TYPE");
        return 0;
    }
    if ( !strncmp(str,"%FILTER",len) )
    {
        tok->setter = filters_set_filter;
        tok->comparator = filters_cmp_filter;
        tok->tag = strdup("%FILTER");
        return 0;
    }

    // is this one of the VCF tags? For now do only INFO and QUAL, to be extended...
    kstring_t tmp = {0,0,0};
    kputsn(str, len, &tmp);

    tok->hdr_id = bcf_hdr_id2int(filter->hdr, BCF_DT_ID, tmp.s);
    if ( tok->hdr_id>=0 ) 
    {
        if ( bcf_hdr_id2type(filter->hdr,BCF_HL_INFO,tok->hdr_id) == BCF_HT_FLAG )
            tok->setter = filters_set_info_flag;
        else
            tok->setter = filters_set_info;
        tok->tag = strdup(tmp.s);
        if ( tmp.s ) free(tmp.s);
        return 0;
    }

    // is it a substrict VCF vector tag?
    if ( tmp.s[tmp.l-1] == ']' )
    {
        int i;
        for (i=0; i<tmp.l; i++)
            if ( tmp.s[i]=='[' ) { tmp.s[i] = 0; break; }

        tok->hdr_id = bcf_hdr_id2int(filter->hdr, BCF_DT_ID, tmp.s);
        if ( tok->hdr_id>=0 )
        {
            switch ( bcf_hdr_id2type(filter->hdr,BCF_HL_INFO,tok->hdr_id) ) 
            {
                case BCF_HT_INT:  tok->setter = &filters_set_info_int; break;
                case BCF_HT_REAL: tok->setter = &filters_set_info_float; break;
                default: error("FIXME: not ready for this, sorry\n");
            }
            tok->idx = atoi(&tmp.s[i+1]);
            if ( tmp.s ) free(tmp.s);
            return 0;
        }
    }

    // is it a value?
    char *end;
    errno = 0;
    tok->threshold = strtod(tmp.s, &end);
    if ( errno!=0 || end==tmp.s ) error("[%s:%d %s] Error: the tag \"INFO/%s\" is not defined in the VCF header\n", __FILE__,__LINE__,__FUNCTION__,tmp.s);

    if ( tmp.s ) free(tmp.s);
    return 0;
}


void filter_debug_print(token_t *toks, int ntoks)
{
    int i;
    for (i=0; i<ntoks; i++)
    {
        token_t *tok = &toks[i];
        if ( tok->tok_type==TOK_VAL )
        {
            if ( tok->key )
                fprintf(stderr,"%s", tok->key);
            else if ( tok->tag )
                fprintf(stderr,"%s", tok->tag);
            else
                fprintf(stderr,"%e", tok->threshold);
        }
        else
            fprintf(stderr,"%c", "x()[<=>]!|&+-*/"[tok->tok_type]);
        if ( tok->setter ) fprintf(stderr,"\t[setter %p]", tok->setter);
        fprintf(stderr,"\n");
    }
}

// Parse filter expression and convert to reverse polish notation. Dijkstra's shunting-yard algorithm
filter_t *filter_init(bcf_hdr_t *hdr, const char *str)
{
    filter_t *filter = (filter_t *) calloc(1,sizeof(filter_t));
    filter->str = strdup(str);
    filter->hdr = hdr;

    int nops = 0, mops = 0, *ops = NULL;    // operators stack
    int nout = 0, mout = 0;                 // filter tokens, RPN
    token_t *out = NULL;
    char *tmp = filter->str;
    int last_op = -1;
    while ( *tmp )
    {
        int len, ret;
        ret = filters_next_token(&tmp, &len);
        if ( ret==-1 ) error("Missing quotes in: %s\n", str);

        // fprintf(stderr,"token=[%c] .. [%s] %d\n", "x()[<=>]!|&+-*/"[ret], tmp, len);
        // int i; for (i=0; i<nops; i++) fprintf(stderr," .%c.", "x()[<=>]!|&+-*/"[ops[i]]); fprintf(stderr,"\n");

        if ( ret==TOK_LFT )         // left bracket
        {
            nops++;
            hts_expand(int, nops, mops, ops);
            ops[nops-1] = ret;
        }
        else if ( ret==TOK_RGT )    // right bracket
        {
            while ( nops>0 && ops[nops-1]!=TOK_LFT )
            {
                nout++;
                hts_expand0(token_t, nout, mout, out);
                out[nout-1].tok_type = ops[nops-1];
                nops--;
            }
            if ( nops<=0 ) error("Could not parse: %s\n", str);
            nops--;
        }
        else if ( ret!=TOK_VAL )    // one of the operators
        {
            // detect unary minus: replace -value with -1*(value)
            if ( ret==TOK_SUB && last_op!=TOK_VAL && last_op!=TOK_RGT )
            {
                nout++;
                hts_expand0(token_t, nout, mout, out);
                token_t *tok = &out[nout-1];
                tok->tok_type  = TOK_VAL;
                tok->hdr_id    = -1;
                tok->pass      = -1;
                tok->threshold = -1.0;
                ret = TOK_MULT;
            }
            else
            {
                while ( nops>0 && op_prec[ret] < op_prec[ops[nops-1]] )
                {
                    nout++;
                    hts_expand0(token_t, nout, mout, out);
                    out[nout-1].tok_type = ops[nops-1];
                    nops--;
                }
            }
            nops++;
            hts_expand(int, nops, mops, ops);
            ops[nops-1] = ret;
        }
        else if ( !len ) 
        {
            if ( *tmp && !isspace(*tmp) ) error("Could not parse the expression: [%s]\n", str);
            break;     // all tokens read
        }
        else           // annotation name or filtering value
        {
            nout++;
            hts_expand0(token_t, nout, mout, out);
            filters_init1(filter, tmp, len, &out[nout-1]);
            tmp += len;
        }
        last_op = ret;
    }
    while ( nops>0 )
    {
        if ( ops[nops-1]==TOK_LFT || ops[nops-1]==TOK_RGT ) error("Could not parse the expression: [%s]\n", filter->str);
        nout++;
        hts_expand0(token_t, nout, mout, out);
        out[nout-1].tok_type = ops[nops-1];
        nops--;
    }

    // In the special cases of %TYPE and %FILTER the BCF header IDs are yet unknown. Walk through the
    // list of operators and convert the strings (e.g. "PASS") to BCF ids. The string value token must be
    // just before or after the %FILTER token and they must be followed with a comparison operator.
    // This code is fragile: improve me.
    int i;
    for (i=0; i<nout; i++)
    {
        if ( out[i].tok_type!=TOK_VAL ) continue;
        if ( !out[i].tag ) continue;
        if ( !strcmp(out[i].tag,"%TYPE") )
        {
            if ( i+1==nout ) error("Could not parse the expression: %s\n", filter->str);
            int j = i+1;
            if ( out[j].tok_type==TOK_EQ || out[j].tok_type==TOK_NE ) j = i - 1;
            if ( out[j].tok_type!=TOK_VAL || !out[j].key ) error("[%s:%d %s] Could not parse the expression: %s\n",  __FILE__,__LINE__,__FUNCTION__, filter->str);
            if ( !strcasecmp(out[j].key,"snp") || !strcasecmp(out[j].key,"snps") ) out[j].threshold = VCF_SNP;
            else if ( !strcasecmp(out[j].key,"indel") || !strcasecmp(out[j].key,"indels") ) out[j].threshold = VCF_INDEL;
            else if ( !strcasecmp(out[j].key,"mnp") || !strcasecmp(out[j].key,"mnps") ) out[j].threshold = VCF_MNP;
            else if ( !strcasecmp(out[j].key,"other") ) out[j].threshold = VCF_OTHER;
            else error("The type \"%s\" not recognised: %s\n", out[j].key, filter->str);
            out[j].tag = out[j].key; out[j].key = NULL;
            i = j;
            continue;
        }
        if ( !strcmp(out[i].tag,"%FILTER") )
        {
            if ( i+1==nout ) error("Could not parse the expression: %s\n", filter->str);
            int j = i+1;
            if ( out[j].tok_type==TOK_EQ || out[j].tok_type==TOK_NE ) j = i - 1;
            if ( out[j].tok_type!=TOK_VAL || !out[j].key ) error("[%s:%d %s] Could not parse the expression: %s\n", __FILE__,__LINE__,__FUNCTION__, filter->str);
            if ( strcmp(".",out[j].key) )
            {
                out[j].hdr_id = bcf_hdr_id2int(filter->hdr, BCF_DT_ID, out[j].key);
                if ( !bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_FLT,out[j].hdr_id) )
                    error("The filter \"%s\" not present in the VCF header\n", out[j].key);
            }
            else
                out[j].hdr_id = -1;
            out[j].tag = out[j].key; out[j].key = NULL;
            out[i].hdr_id = out[j].hdr_id;
            i = j;
            continue;
        }
    }

    // filter_debug_print(out, nout);

    if ( mops ) free(ops);
    filter->filters   = out;
    filter->nfilters  = nout;
    filter->flt_stack = (token_t **)malloc(sizeof(token_t*)*nout);
    return filter;
}

void filter_destroy(filter_t *filter)
{
    int i;
    for (i=0; i<filter->nfilters; i++)
    {
        free(filter->filters[i].key);
        free(filter->filters[i].tag);
    }
    free(filter->filters);
    free(filter->flt_stack);
    free(filter->str);
    free(filter);
}

int filter_test(filter_t *filter, bcf1_t *line)
{
    bcf_unpack(line, BCF_UN_INFO);

    int i, nstack = 0;
    for (i=0; i<filter->nfilters; i++)
    {
        filter->filters[i].missing_value = 0;
        filter->filters[i].str_value = NULL;
        filter->filters[i].pass = -1;

        if ( filter->filters[i].tok_type == TOK_VAL )
        {
            if ( filter->filters[i].setter ) 
                filter->filters[i].setter(line, &filter->filters[i]);
            else if ( filter->filters[i].key )
            {
                filter->filters[i].str_value = filter->filters[i].key;
                filter->filters[i].num_value = filter->filters[i].num_value;
            }
            else
                filter->filters[i].num_value = filter->filters[i].threshold;
            filter->flt_stack[nstack++] = &filter->filters[i];
            continue;
        }
        if ( nstack<2 ) 
            error("Error occurred while processing the filter \"%s\": too few values left on stack (%d)\n", filter->str,nstack);

        int is_str  = (filter->flt_stack[nstack-1]->str_value ? 1 : 0) + (filter->flt_stack[nstack-2]->str_value ? 1 : 0 );

        if ( filter->filters[i].tok_type == TOK_OR )
        {
            if ( filter->flt_stack[nstack-1]->pass<0 || filter->flt_stack[nstack-2]->pass<0 ) 
                error("Error occurred while processing the filter \"%s\" (%d %d OR)\n", filter->str,filter->flt_stack[nstack-2]->pass,filter->flt_stack[nstack-1]->pass);
            filter->flt_stack[nstack-2]->pass = filter->flt_stack[nstack-1]->pass + filter->flt_stack[nstack-2]->pass;
            nstack--;
            continue;
        }
        if ( filter->filters[i].tok_type == TOK_AND )
        {
            if ( filter->flt_stack[nstack-1]->pass<0 || filter->flt_stack[nstack-2]->pass<0 ) 
                error("Error occurred while processing the filter \"%s\" (%d %d AND)\n", filter->str,filter->flt_stack[nstack-2]->pass,filter->flt_stack[nstack-1]->pass);
            filter->flt_stack[nstack-2]->pass = filter->flt_stack[nstack-1]->pass * filter->flt_stack[nstack-2]->pass;
            nstack--;
            continue;
        }

        if ( filter->filters[i].tok_type == TOK_ADD )
        {
            filter->flt_stack[nstack-2]->num_value += filter->flt_stack[nstack-1]->num_value;
            nstack--;
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_SUB )
        {
            filter->flt_stack[nstack-2]->num_value -= filter->flt_stack[nstack-1]->num_value;
            nstack--;
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_MULT )
        {
            filter->flt_stack[nstack-2]->num_value *= filter->flt_stack[nstack-1]->num_value;
            nstack--;
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_DIV )
        {
            filter->flt_stack[nstack-2]->num_value /= filter->flt_stack[nstack-1]->num_value;
            nstack--;
            continue;
        }

        int is_true = 0;
        if ( filter->flt_stack[nstack-1]->missing_value || filter->flt_stack[nstack-2]->missing_value )
            is_true = 0;
        else if ( filter->filters[i].tok_type == TOK_EQ )
        {
            if ( filter->flt_stack[nstack-1]->comparator )
                is_true = filter->flt_stack[nstack-1]->comparator(filter->flt_stack[nstack-1],filter->flt_stack[nstack-2],TOK_EQ,line);
            else if ( filter->flt_stack[nstack-2]->comparator )
                is_true = filter->flt_stack[nstack-2]->comparator(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],TOK_EQ,line);
            else if ( is_str==2 ) 
            {
                int ncmp = filter->flt_stack[nstack-1]->num_value > filter->flt_stack[nstack-2]->num_value ? filter->flt_stack[nstack-1]->num_value : filter->flt_stack[nstack-2]->num_value;
                is_true = strncmp(filter->flt_stack[nstack-1]->str_value,filter->flt_stack[nstack-2]->str_value, ncmp) ? 0 : 1;
            }
            else if ( is_str==1 ) 
                error("Comparing string to numeric value: %s\n", filter->str);
            else
                is_true = (filter->flt_stack[nstack-1]->num_value == filter->flt_stack[nstack-2]->num_value) ? 1 : 0;
        }
        else if ( filter->filters[i].tok_type == TOK_NE )
        {
            if ( filter->flt_stack[nstack-1]->comparator ) 
                is_true = filter->flt_stack[nstack-1]->comparator(filter->flt_stack[nstack-1],filter->flt_stack[nstack-2],TOK_NE,line);
            else if ( filter->flt_stack[nstack-2]->comparator )
                is_true = filter->flt_stack[nstack-2]->comparator(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],TOK_NE,line);
            else if ( is_str==2 )
            {
                int ncmp = filter->flt_stack[nstack-1]->num_value > filter->flt_stack[nstack-2]->num_value ? filter->flt_stack[nstack-1]->num_value : filter->flt_stack[nstack-2]->num_value;
                is_true = strncmp(filter->flt_stack[nstack-1]->str_value,filter->flt_stack[nstack-2]->str_value, ncmp) ? 1 : 0;
            }
            else if ( is_str==1 )
                error("Comparing string to numeric value: %s\n", filter->str);
            else
                is_true = ( filter->flt_stack[nstack-2]->num_value != filter->flt_stack[nstack-1]->num_value ) ? 1 : 0;
        }
        else if ( is_str>0 ) error("Wrong operator in string comparison: %s [%s,%s]\n", filter->str, filter->flt_stack[nstack-1]->str_value, filter->flt_stack[nstack-2]->str_value);
        else if ( filter->filters[i].tok_type == TOK_LE )
            is_true = ( filter->flt_stack[nstack-2]->num_value <= filter->flt_stack[nstack-1]->num_value ) ? 1 : 0;
        else if ( filter->filters[i].tok_type == TOK_LT )
            is_true = ( filter->flt_stack[nstack-2]->num_value <  filter->flt_stack[nstack-1]->num_value ) ? 1 : 0;
        else if ( filter->filters[i].tok_type == TOK_EQ )
            is_true = ( filter->flt_stack[nstack-2]->num_value == filter->flt_stack[nstack-1]->num_value ) ? 1 : 0;
        else if ( filter->filters[i].tok_type == TOK_BT )
            is_true = ( filter->flt_stack[nstack-2]->num_value >  filter->flt_stack[nstack-1]->num_value ) ? 1 : 0;
        else if ( filter->filters[i].tok_type == TOK_BE )
            is_true = ( filter->flt_stack[nstack-2]->num_value >= filter->flt_stack[nstack-1]->num_value ) ? 1 : 0;
        else
            error("FIXME: did not expect this .. tok_type %d = %d\n", i, filter->filters[i].tok_type);

        filter->flt_stack[nstack-2]->pass = is_true;
        nstack--;
    }
    if ( nstack>1 ) error("Error occurred while processing the filter \"%s\": too many values left on stack (%d)\n", filter->str,nstack);
    return filter->flt_stack[0]->pass;
}

