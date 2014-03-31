#include <ctype.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "filter.h"
#include "bcftools.h"

typedef struct _token_t
{
    // read-only values, same for all VCF lines
    int tok_type;       // one of the TOK_* keys below
    char *key;          // set only for string constants, otherwise NULL
    char *tag;          // for debugging and printout only, VCF tag name
    float threshold;    // filtering threshold
    int hdr_id;         // BCF header lookup ID
    int idx;            // 1-based index to VCF vector
    void (*setter)(filter_t *, bcf1_t *, struct _token_t *);
    int (*comparator)(struct _token_t *, struct _token_t *, int op_type, bcf1_t *);

    // modified on filter evaluation at each VCF line
    float *values;      // In case str_value is set, values[0] is one sample's string length,
    char *str_value;    //  values[0]*nvalues gives the total length;
    int is_str;
    int pass;               // -1 not applicable, 0 fails, >0 pass
    uint8_t *pass_samples;  // status of individual samples
    int missing_value;
    int nvalues, mvalues;   // dimension of values
}
token_t;

struct _filter_t 
{
    bcf_hdr_t *hdr;
    char *str;
    int nfilters;
    token_t *filters, **flt_stack;  // filtering input tokens (in RPN) and evaluation stack
    int32_t *tmpi;
    int max_unpack, mtmpi, npass_samples;
};


#define TOK_VAL     0
#define TOK_LFT     1       // (
#define TOK_RGT     2       // )
#define TOK_LE      3       // less or equal
#define TOK_LT      4       // less than
#define TOK_EQ      5       // equal
#define TOK_BT      6       // bigger than
#define TOK_BE      7       // bigger or equal
#define TOK_NE      8       // not equal
#define TOK_OR      9       // |, ||   
#define TOK_AND     10      // &, &&
#define TOK_ADD     11      // +
#define TOK_SUB     12      // -
#define TOK_MULT    13      // *
#define TOK_DIV     14      // /
#define TOK_MAX     15
#define TOK_MIN     16
#define TOK_AVG     17
#define TOK_FUNC    18

//                      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
//                        ( ) [ < = > ] ! | &  +  -  *  /  M  m  a
static int op_prec[] = {0,1,1,5,5,5,5,5,5,2,3, 6, 6, 7, 7, 8, 8, 8};
#define TOKEN_STRING "x()[<=>]!|&+-*/Mmaf"

static int filters_next_token(char **str, int *len)
{
    char *tmp = *str;
    while ( *tmp && isspace(*tmp) ) tmp++;
    *str = tmp;
    *len = 0;

    // test for doubles: d.ddde[+-]dd
    if ( isdigit(*str[0]) || *str[0]=='.' )   // strtod would eat +/-
    {
        double val_unused v = strtod(*str, &tmp);
        if ( *str!=tmp && (!tmp[0] || !isalnum(tmp[0])) )
        {
            *len = tmp - (*str);
            return TOK_VAL;
        }
        tmp = *str;
    }

    if ( !strncmp(tmp,"%MAX(",5) ) { (*str) += 4; return TOK_MAX; }
    if ( !strncmp(tmp,"%MIN(",5) ) { (*str) += 4; return TOK_MIN; }
    if ( !strncmp(tmp,"%AVG(",5) ) { (*str) += 4; return TOK_AVG; }
    if ( !strncmp(tmp,"INFO/",5) ) tmp += 5;
    if ( !strncmp(tmp,"FORMAT/",7) ) tmp += 7;
    if ( !strncmp(tmp,"FMT/",4) ) tmp += 4;

    while ( tmp[0] )
    {
        if ( tmp[0]=='"' ) break;
        if ( tmp[0]=='\'' ) break;
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
    if ( tmp[0]=='"' || tmp[0]=='\'' )
    {
        int quote = tmp[0];
        tmp++;
        while ( *tmp && tmp[0]!=quote ) tmp++;
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

static void filters_set_qual(filter_t *flt, bcf1_t *line, token_t *tok)
{
    float *ptr = &line->qual;
    if ( bcf_float_is_missing(*ptr) )
    {
        tok->missing_value = 1;
        tok->nvalues = 0;
    }
    else
    {
        tok->values[0] = line->qual;
        tok->nvalues = 1;
    }
}
static void filters_set_type(filter_t *flt, bcf1_t *line, token_t *tok)
{
    tok->values[0] = bcf_get_variant_types(line);
    tok->nvalues = 1;
}
static void filters_set_info(filter_t *flt, bcf1_t *line, token_t *tok)
{
    assert( tok->hdr_id >=0  );
    int i;
    for (i=0; i<line->n_info; i++)
        if ( line->d.info[i].key == tok->hdr_id ) break;

    if ( i==line->n_info ) 
    {
        tok->missing_value = 1;
        tok->nvalues = 0;
    }
    else if ( line->d.info[i].type==BCF_BT_CHAR )
    {
        tok->str_value = (char*)line->d.info[i].vptr;       // string is typically not null-terminated
        tok->values[0] = line->d.info[i].len;
        tok->nvalues   = 1;
    }
    else if ( line->d.info[i].type==BCF_BT_FLOAT )
    {
        tok->values[0] = line->d.info[i].v1.f;
        tok->str_value = NULL;
        tok->nvalues   = 1;
    }
    else
    {
        tok->values[0] = line->d.info[i].v1.i;
        tok->str_value = NULL;
        tok->nvalues   = 1;
    }
}
static void filters_set_filter(filter_t *flt, bcf1_t *line, token_t *tok) {}  // dummy
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
static int bcf_get_info_value(bcf1_t *line, int info_id, int ivec, void *value)
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

static void filters_set_info_int(filter_t *flt, bcf1_t *line, token_t *tok)
{
    int value;
    if ( bcf_get_info_value(line,tok->hdr_id,tok->idx,&value) <= 0 )
    {
        tok->missing_value = 1;
        tok->nvalues = 0;
    }
    else
    {
        tok->values[0] = value;
        tok->nvalues = 1;
    }
}

static void filters_set_info_float(filter_t *flt, bcf1_t *line, token_t *tok)
{
    float value;
    if ( bcf_get_info_value(line,tok->hdr_id,tok->idx,&value) <= 0 )
    {
        tok->missing_value = 1;
        tok->nvalues = 0;
    }
    else
    {
        tok->values[0] = value;
        tok->nvalues = 1;
    }
}

static void filters_set_info_flag(filter_t *flt, bcf1_t *line, token_t *tok)
{
    int j;
    for (j=0; j<line->n_info; j++)
        if ( line->d.info[j].key == tok->hdr_id ) break;
    tok->values[0] = j==line->n_info ? 0 : 1;
    tok->nvalues = 1;
}

static void filters_set_format_int(filter_t *flt, bcf1_t *line, token_t *tok)
{
    int i;
    if ( (tok->nvalues=bcf_get_format_int32(flt->hdr,line,tok->tag,&flt->tmpi,&flt->mtmpi))<0 )
    {
        tok->missing_value = 1;
        tok->nvalues = 0;
    }
    else
    {
        int is_missing = 1;
        hts_expand(float,tok->nvalues,tok->mvalues,tok->values);
        for (i=0; i<tok->nvalues; i++)
        {
            if ( flt->tmpi[i]==bcf_int32_missing || flt->tmpi[i]==bcf_int32_vector_end )
                bcf_float_set_missing(tok->values[i]);
            else
            {
                tok->values[i] = flt->tmpi[i];
                is_missing = 0;
            }
        }
        tok->missing_value = is_missing;
    }
}
static void filters_set_format_float(filter_t *flt, bcf1_t *line, token_t *tok)
{
    if ( (tok->nvalues=bcf_get_format_float(flt->hdr,line,tok->tag,&tok->values,&tok->mvalues))<=0 ) 
        tok->missing_value = 1;
}
static void filters_set_format_string(filter_t *flt, bcf1_t *line, token_t *tok)
{
    int ndim = 0;
    bcf_get_format_char(flt->hdr,line,tok->tag,&tok->str_value,&ndim);
    tok->nvalues = bcf_hdr_nsamples(flt->hdr);
    tok->values[0] = ndim;
    tok->missing_value = ndim ? 0 : 1;
}


static void set_max(filter_t *flt, bcf1_t *line, token_t *tok) 
{ 
    float val = -HUGE_VAL;
    int i;
    for (i=0; i<tok->nvalues; i++) 
    {
        if ( !bcf_float_is_missing(tok->values[i]) && val < tok->values[i] ) val = tok->values[i];
    }
    tok->values[0] = val;
    tok->nvalues = 1;
}
static void set_min(filter_t *flt, bcf1_t *line, token_t *tok) 
{ 
    float val = HUGE_VAL;
    int i;
    for (i=0; i<tok->nvalues; i++) 
        if ( !bcf_float_is_missing(tok->values[i]) && val > tok->values[i] ) val = tok->values[i];
    tok->values[0] = val;
    tok->nvalues = 1;
}
static void set_avg(filter_t *flt, bcf1_t *line, token_t *tok) 
{ 
    float val = 0;
    int i, n = 0;
    for (i=0; i<tok->nvalues; i++) 
        if ( !bcf_float_is_missing(tok->values[i]) ) val += tok->values[i];
    tok->values[0] = n ? val / n : 0;
    tok->nvalues = 1;
}
#define VECTOR_ARITHMETICS(atok,btok,AOP) \
{ \
    int i, is_missing = 1; \
    if ( !(atok)->missing_value && !(btok)->missing_value ) \
    { \
        if ( (atok)->nvalues==(btok)->nvalues ) \
        { \
            for (i=0; i<(atok)->nvalues; i++) \
            { \
                if ( bcf_float_is_missing((atok)->values[i]) ) continue; \
                if ( bcf_float_is_missing((btok)->values[i]) ) { bcf_float_set_missing((atok)->values[i]); continue; } \
                is_missing = 0; \
                (atok)->values[i] = (atok)->values[i] AOP (btok)->values[i]; \
            } \
        } \
        else if ( (atok)->nvalues==1 ) \
        { \
            hts_expand(float,(btok)->nvalues,(atok)->mvalues,(atok)->values); \
            for (i=0; i<(btok)->nvalues; i++) \
            { \
                if ( bcf_float_is_missing((atok)->values[0]) || bcf_float_is_missing((btok)->values[i]) ) \
                { \
                    bcf_float_set_missing((atok)->values[i]); \
                    continue; \
                } \
                is_missing = 0; \
                (atok)->values[i] = (atok)->values[0] AOP (btok)->values[i]; \
            } \
            (atok)->nvalues = (btok)->nvalues; \
        } \
        else if ( (btok)->nvalues==1 ) \
        { \
            for (i=0; i<(atok)->nvalues; i++) \
            { \
                if ( bcf_float_is_missing((atok)->values[i]) || bcf_float_is_missing((btok)->values[0]) ) \
                { \
                    bcf_float_set_missing((atok)->values[i]); \
                    continue; \
                } \
                is_missing = 0; \
                (atok)->values[i] = (atok)->values[i] AOP (btok)->values[0]; \
            } \
        } \
        else error("[%s:%d %s] fixme: Cannot operate on vectors of different length\n", __FILE__,__LINE__,__FUNCTION__); \
    } \
    (atok)->missing_value = is_missing; \
}

static int vector_logic_and(filter_t *filter, token_t *atok, token_t *btok, uint8_t *samples)
{
    int i, pass = 0;
    if ( atok->missing_value || btok->missing_value )
    {
        for (i=0; i<filter->npass_samples; i++) samples[i] = 0;
        return 0;
    }
    for (i=0; i<filter->npass_samples; i++)
    {
        samples[i] = atok->pass_samples[i] && btok->pass_samples[i];
        if ( !pass && samples[i] ) pass = 1;
    }
    return pass;
}
static int vector_logic_or(filter_t *filter, token_t *atok, token_t *btok, uint8_t *samples)
{
    int i, pass = 0;
    if ( atok->missing_value && btok->missing_value )
    {
        for (i=0; i<filter->npass_samples; i++) samples[i] = 0;
        return 0;
    }
    if ( atok->missing_value )
    {
        for (i=0; i<filter->npass_samples; i++) 
        {
            samples[i] = btok->pass_samples[i];
            if ( !pass && samples[i] ) pass = 1;
        }
    }
    else if ( btok->missing_value )
    {
        for (i=0; i<filter->npass_samples; i++) 
        {
            samples[i] = atok->pass_samples[i];
            if ( !pass && samples[i] ) pass = 1;
        }
    }
    else
    {
        for (i=0; i<filter->npass_samples; i++)
        {
            samples[i] = atok->pass_samples[i] || btok->pass_samples[i];
            if ( !pass && samples[i] ) pass = 1;
        }
    }
    return pass;
}

#define CMP_VECTORS(atok,btok,smpls_pass,CMP_OP,ret) \
{ \
    int i, pass = 0, is_missing = 1; \
    if ( (atok)->missing_value || (btok)->missing_value ) (ret) = 0; \
    else \
    { \
        if ( (atok)->nvalues==(btok)->nvalues ) \
        { \
            for (i=0; i<(atok)->nvalues; i++) \
            { \
                if ( bcf_float_is_missing((atok)->values[i]) ) { (smpls_pass)[i] = 0; continue; } \
                if ( bcf_float_is_missing((btok)->values[i]) ) { (smpls_pass)[i] = 0; continue; } \
                is_missing = 0; \
                if ( (atok)->values[i] CMP_OP (btok)->values[i] ) { (smpls_pass)[i] = 1; pass = 1; } \
                else (smpls_pass)[i] = 0; \
            } \
        } \
        else if ( (atok)->nvalues==1 ) \
        { \
            for (i=0; i<(btok)->nvalues; i++) \
            { \
                if ( bcf_float_is_missing((atok)->values[0]) ) { (smpls_pass)[i] = 0; continue; } \
                if ( bcf_float_is_missing((btok)->values[i]) ) { (smpls_pass)[i] = 0; continue; } \
                is_missing = 0; \
                if ( (atok)->values[0] CMP_OP (btok)->values[i] ) { (smpls_pass)[i] = 1; pass = 1; } \
                else (smpls_pass)[i] = 0; \
            } \
        } \
        else if ( (btok)->nvalues==1 ) \
        { \
            for (i=0; i<(atok)->nvalues; i++) \
            { \
                if ( bcf_float_is_missing((atok)->values[i]) ) { (smpls_pass)[i] = 0; continue; } \
                if ( bcf_float_is_missing((btok)->values[0]) ) { (smpls_pass)[i] = 0; continue; } \
                is_missing = 0; \
                if ( (atok)->values[i] CMP_OP (btok)->values[0] ) { (smpls_pass)[i] = 1; pass = 1; } \
                else (smpls_pass)[i] = 0; \
            } \
        } \
        else error("[%s:%d %s] fixme: Cannot compare vectors of different length\n", __FILE__,__LINE__,__FUNCTION__); \
        (ret) = pass; \
    } \
    (atok)->missing_value = is_missing; \
}
static int cmp_vector_strings(token_t *atok, token_t *btok, uint8_t *samples, int logic)    // logic: TOK_EQ or TOK_NE
{
    if ( atok->missing_value ) return 0;
    if ( btok->missing_value ) { atok->missing_value = 1; return 0; }
    int i, pass = 0;
    if ( atok->nvalues==btok->nvalues )
    {
        char *astr = atok->str_value, *bstr = btok->str_value;
        for (i=0; i<atok->nvalues; i++)
        {
            char *aend = astr + (int)atok->values[0], *a = astr;
            while ( a<aend && *a ) a++;
            char *bend = bstr + (int)btok->values[0], *b = bstr;
            while ( b<bend && *b ) b++;
            if ( a-astr != b-bstr ) samples[i] = logic==TOK_EQ ? 0 : 1;
            else samples[i] = logic==TOK_EQ ? !strncmp(astr,bstr,a-astr) : strncmp(astr,bstr,a-astr);
            if ( !pass && samples[i] ) pass = 1;
            while ( astr<aend && !*astr ) astr++;
            while ( bstr<bend && !*bstr ) bstr++;
        }
    }
    else if ( atok->nvalues==1 || btok->nvalues==1 )
    {
        token_t *xtok, *ytok;
        if ( atok->nvalues==1 ) { xtok = atok; ytok = btok; }
        else { xtok = btok; ytok = atok; }
        char *xstr = xtok->str_value, *ystr = ytok->str_value;
        char *xend = xstr + (int)xtok->values[0], *x = xstr;
        while ( x<xend && *x ) x++;
        for (i=0; i<ytok->nvalues; i++)
        {
            char *yend = ystr + (int)ytok->values[0], *y = ystr;
            while ( y<yend && *y ) y++;
            if ( x-xstr != y-ystr ) samples[i] = logic==TOK_EQ ? 0 : 1;
            else samples[i] = logic==TOK_EQ ? !strncmp(xstr,ystr,x-xstr) : strncmp(xstr,ystr,x-xstr);
            if ( !pass && samples[i] ) pass = 1;
            while ( ystr<yend && !*ystr ) ystr++;
        }
    }
    else error("[%s:%d %s] todo: Cannot compare vectors of different length\n", __FILE__,__LINE__,__FUNCTION__);
    return pass;
}

static int filters_init1(filter_t *filter, char *str, int len, int inside_func, token_t *tok)
{
    tok->tok_type = TOK_VAL;
    tok->hdr_id   = -1;
    tok->pass     = -1;

    // is this a string constant?
    if ( str[0]=='"' || str[0]=='\'' )
    {
        int quote = str[0];
        if ( str[len-1] != quote ) error("TODO: [%s]\n", filter->str);
        tok->key = (char*) calloc(len-1,sizeof(char));
        hts_expand(float,1,tok->mvalues,tok->values);
        tok->values[0] = len-2;
        memcpy(tok->key,str+1,len-2);
        tok->key[len-2] = 0;
        tok->is_str = 1;
        return 0;
    }

    int is_fmt = -1;
    if ( !strncmp(str,"FMT/",4) || !strncmp(str,"FORMAT/",4) ) { str += 4; len -= 4; is_fmt = 1; }
    else
    {
        if ( !strncmp(str,"INFO/",5) ) { is_fmt = 0; str += 5; len -= 5; }
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
    }
    if ( is_fmt==-1 ) is_fmt = inside_func ? 1 : 0;
    if ( is_fmt ) filter->max_unpack |= BCF_UN_FMT;

    // is this one of the VCF tags? For now do only INFO and QUAL
    kstring_t tmp = {0,0,0};
    kputsn(str, len, &tmp);

    tok->hdr_id = bcf_hdr_id2int(filter->hdr,BCF_DT_ID,tmp.s);
    if ( tok->hdr_id>=0 ) 
    {
        if ( is_fmt )
        {
            if ( !bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_FMT,tok->hdr_id) )
                error("No such FORMAT field: %s\n", tmp.s);
            if ( bcf_hdr_id2number(filter->hdr,BCF_HL_FMT,tok->hdr_id)!=1 )
                error("Error: Arrays must be subscripted, e.g. %s[0]\n", tmp.s);
            switch ( bcf_hdr_id2type(filter->hdr,BCF_HL_FMT,tok->hdr_id) )
            {
                case BCF_HT_INT:  tok->setter = &filters_set_format_int; break;
                case BCF_HT_REAL: tok->setter = &filters_set_format_float; break;
                case BCF_HT_STR:  tok->setter = &filters_set_format_string; tok->is_str = 1; break;
                default: error("[%s:%d %s] FIXME\n", __FILE__,__LINE__,__FUNCTION__);
            }
        }
        else if ( !bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_INFO,tok->hdr_id) )
            error("No such INFO field: %s\n", tmp.s);
        else
        {
            if ( bcf_hdr_id2type(filter->hdr,BCF_HL_INFO,tok->hdr_id) == BCF_HT_FLAG )
                tok->setter = filters_set_info_flag;
            else
            {
                if ( bcf_hdr_id2type(filter->hdr,BCF_HL_INFO,tok->hdr_id) == BCF_HT_STR ) tok->is_str = 1;
                if ( bcf_hdr_id2number(filter->hdr,BCF_HL_INFO,tok->hdr_id)!=1 )
                    error("Error: Arrays must be subscripted, e.g. %s[0]\n", tmp.s);
                tok->setter = filters_set_info;
            }
            filter->max_unpack |= BCF_UN_INFO;
        }
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
            if ( is_fmt )
            {
                if ( !bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_FMT,tok->hdr_id) )
                    error("No such FORMAT field: %s\n", tmp.s);
                if ( bcf_hdr_id2number(filter->hdr,BCF_HL_FMT,tok->hdr_id)!=1 )
                    error("Error: Arrays must be subscripted, e.g. %s[0]\n", tmp.s);
                switch ( bcf_hdr_id2type(filter->hdr,BCF_HL_FMT,tok->hdr_id) )
                {
                    case BCF_HT_INT:  tok->setter = &filters_set_format_int; break;
                    case BCF_HT_REAL: tok->setter = &filters_set_format_float; break;
                    case BCF_HT_STR:  tok->setter = &filters_set_format_string; tok->is_str = 1; break;
                    default: error("[%s:%d %s] FIXME\n", __FILE__,__LINE__,__FUNCTION__);
                }
            }
            else if ( !bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_INFO,tok->hdr_id) )
                error("No such INFO field: %s\n", tmp.s);
            else
            {
                switch ( bcf_hdr_id2type(filter->hdr,BCF_HL_INFO,tok->hdr_id) )
                {
                    case BCF_HT_INT:  tok->setter = &filters_set_info_int; break;
                    case BCF_HT_REAL: tok->setter = &filters_set_info_float; break;
                    case BCF_HT_STR:  error("fixme: String vectors not supported yet\n"); break;
                    default: error("[%s:%d %s] FIXME\n", __FILE__,__LINE__,__FUNCTION__);
                }
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
    if ( errno!=0 || end!=tmp.s+len ) error("[%s:%d %s] Error: the tag \"INFO/%s\" is not defined in the VCF header\n", __FILE__,__LINE__,__FUNCTION__,tmp.s);

    if ( tmp.s ) free(tmp.s);
    return 0;
}


static void filter_debug_print(token_t *toks, token_t **tok_ptrs, int ntoks)
{
    int i;
    for (i=0; i<ntoks; i++)
    {
        token_t *tok = toks ? &toks[i] : tok_ptrs[i];
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
            fprintf(stderr,"%c", TOKEN_STRING[tok->tok_type]);
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
    int nfunc = 0; // inside funcs the default of tags is FMT, otherwise INFO
    while ( *tmp )
    {
        int len, ret;
        ret = filters_next_token(&tmp, &len);
        if ( ret==-1 ) error("Missing quotes in: %s\n", str);

        // fprintf(stderr,"token=[%c] .. [%s] %d\n", TOKEN_STRING[ret], tmp, len);
        // int i; for (i=0; i<nops; i++) fprintf(stderr," .%c.", TOKEN_STRING[ops[i]]); fprintf(stderr,"\n");

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
                    if ( ops[nops-1]==TOK_MAX || ops[nops-1]==TOK_MIN || ops[nops-1]==TOK_AVG ) nfunc--;
                    nops--;
                }
            }
            nops++;
            hts_expand(int, nops, mops, ops);
            ops[nops-1] = ret;
            if ( ops[nops-1]==TOK_MAX || ops[nops-1]==TOK_MIN || ops[nops-1]==TOK_AVG ) nfunc++;
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
            filters_init1(filter, tmp, len, nfunc, &out[nout-1]);
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
            if ( out[j].tok_type!=TOK_VAL || !out[j].key )
                error("[%s:%d %s] Could not parse the expression, an unquoted string value perhaps? %s\n", __FILE__,__LINE__,__FUNCTION__, filter->str);
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
    filter->npass_samples = filter->max_unpack&BCF_UN_FMT ? bcf_hdr_nsamples(filter->hdr) : 1;
    for (i=0; i<nout; i++)
    {
        if ( out[i].tok_type==TOK_MAX )      { out[i].setter = set_max; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_MIN ) { out[i].setter = set_min; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_AVG ) { out[i].setter = set_avg; out[i].tok_type = TOK_FUNC; }
        hts_expand(float,1,out[i].mvalues,out[i].values);
        out[i].pass_samples = (uint8_t*)malloc(filter->npass_samples);
        int j;
        for (j=0; j<filter->npass_samples; j++) out[i].pass_samples[j] = 1;
    }

    if (0) filter_debug_print(out, NULL, nout);

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
        free(filter->filters[i].str_value);
        free(filter->filters[i].tag);
        free(filter->filters[i].values);
        free(filter->filters[i].pass_samples);
    }
    free(filter->filters);
    free(filter->flt_stack);
    free(filter->str);
    free(filter->tmpi);
    free(filter);
}

int filter_test(filter_t *filter, bcf1_t *line, const uint8_t **samples)
{
    bcf_unpack(line, filter->max_unpack);

    int i, nstack = 0;
    for (i=0; i<filter->nfilters; i++)
    {
        filter->filters[i].missing_value = 0;
        filter->filters[i].pass = -1;

        if ( filter->filters[i].tok_type == TOK_VAL )
        {
            if ( filter->filters[i].setter )    // variable, query the VCF line
                filter->filters[i].setter(filter, line, &filter->filters[i]);
            else if ( filter->filters[i].key )  // string constant
            {
                filter->filters[i].str_value = filter->filters[i].key;
                filter->filters[i].values[0] = filter->filters[i].values[0];
                filter->filters[i].nvalues   = 1;
            }
            else    // numeric constant
            {
                filter->filters[i].values[0] = filter->filters[i].threshold;
                filter->filters[i].nvalues   = 1;
            }

            filter->flt_stack[nstack++] = &filter->filters[i];
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_FUNC ) // all functions take only one argument
        {
            filter->filters[i].setter(filter, line, filter->flt_stack[nstack-1]);
            continue;
        }
        if ( nstack<2 ) 
            error("Error occurred while processing the filter \"%s\" (1:%d)\n", filter->str,nstack);  // too few values left on the stack

        int is_str  = filter->flt_stack[nstack-1]->is_str + filter->flt_stack[nstack-2]->is_str;

        if ( filter->filters[i].tok_type == TOK_OR )
        {
            if ( filter->flt_stack[nstack-1]->pass<0 || filter->flt_stack[nstack-2]->pass<0 ) 
                error("Error occurred while processing the filter \"%s\" (%d %d OR)\n", filter->str,filter->flt_stack[nstack-2]->pass,filter->flt_stack[nstack-1]->pass);
            filter->flt_stack[nstack-2]->pass = vector_logic_or(filter,filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],filter->flt_stack[nstack-2]->pass_samples);
            nstack--;
            continue;
        }
        if ( filter->filters[i].tok_type == TOK_AND )
        {
            if ( filter->flt_stack[nstack-1]->pass<0 || filter->flt_stack[nstack-2]->pass<0 ) 
                error("Error occurred while processing the filter \"%s\" (%d %d AND)\n", filter->str,filter->flt_stack[nstack-2]->pass,filter->flt_stack[nstack-1]->pass);
            filter->flt_stack[nstack-2]->pass = vector_logic_and(filter,filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],filter->flt_stack[nstack-2]->pass_samples);
            nstack--;
            continue;
        }

        if ( filter->filters[i].tok_type == TOK_ADD )
        {
            VECTOR_ARITHMETICS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],+);
            nstack--;
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_SUB )
        {
            VECTOR_ARITHMETICS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],-);
            nstack--;
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_MULT )
        {
            VECTOR_ARITHMETICS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],*);
            nstack--;
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_DIV )
        {
            VECTOR_ARITHMETICS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],/);
            nstack--;
            continue;
        }

        int is_true = 0;
        if ( filter->flt_stack[nstack-1]->missing_value || filter->flt_stack[nstack-2]->missing_value )
        {
            filter->flt_stack[nstack-2]->missing_value = 1;
            is_true = 0;
        }
        else if ( filter->filters[i].tok_type == TOK_EQ )
        {
            if ( filter->flt_stack[nstack-1]->comparator )
                is_true = filter->flt_stack[nstack-1]->comparator(filter->flt_stack[nstack-1],filter->flt_stack[nstack-2],TOK_EQ,line);
            else if ( filter->flt_stack[nstack-2]->comparator )
                is_true = filter->flt_stack[nstack-2]->comparator(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],TOK_EQ,line);
            else if ( is_str==2 )   // both are strings
                is_true = cmp_vector_strings(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],filter->flt_stack[nstack-2]->pass_samples,TOK_EQ);
            else if ( is_str==1 ) 
                error("Comparing string to numeric value: %s\n", filter->str);
            else
                CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],filter->flt_stack[nstack-2]->pass_samples,==,is_true);
        }
        else if ( filter->filters[i].tok_type == TOK_NE )
        {
            if ( filter->flt_stack[nstack-1]->comparator ) 
                is_true = filter->flt_stack[nstack-1]->comparator(filter->flt_stack[nstack-1],filter->flt_stack[nstack-2],TOK_NE,line);
            else if ( filter->flt_stack[nstack-2]->comparator )
                is_true = filter->flt_stack[nstack-2]->comparator(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],TOK_NE,line);
            else if ( is_str==2 )
                is_true = cmp_vector_strings(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],filter->flt_stack[nstack-2]->pass_samples,TOK_NE);
            else if ( is_str==1 )
                error("Comparing string to numeric value: %s\n", filter->str);
            else
                CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],filter->flt_stack[nstack-2]->pass_samples,!=,is_true);
        }
        else if ( is_str>0 ) error("Wrong operator in string comparison: %s [%s,%s]\n", filter->str, filter->flt_stack[nstack-1]->str_value, filter->flt_stack[nstack-2]->str_value);
        else if ( filter->filters[i].tok_type == TOK_LE )
            CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],filter->flt_stack[nstack-2]->pass_samples,<=,is_true)
        else if ( filter->filters[i].tok_type == TOK_LT )
            CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],filter->flt_stack[nstack-2]->pass_samples,<,is_true)
        else if ( filter->filters[i].tok_type == TOK_BT )
            CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],filter->flt_stack[nstack-2]->pass_samples,>,is_true)
        else if ( filter->filters[i].tok_type == TOK_BE )
            CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],filter->flt_stack[nstack-2]->pass_samples,>=,is_true)
        else
            error("FIXME: did not expect this .. tok_type %d = %d\n", i, filter->filters[i].tok_type);

        filter->flt_stack[nstack-2]->pass = is_true;
        nstack--;
    }
    if ( nstack>1 ) error("Error occurred while processing the filter \"%s\" (2:%d)\n", filter->str,nstack);    // too few values left on the stack
    if ( samples ) *samples = filter->max_unpack&BCF_UN_FMT ? filter->flt_stack[0]->pass_samples : NULL;
    return filter->flt_stack[0]->pass;
}

void filter_expression_info(FILE *fp)
{
    fprintf(fp, "Filter expressions may contain\n");
    fprintf(fp, "    - numerical constants and string constants\n");
    fprintf(fp, "        .. 1, 1.0, 1e-4\n");
    fprintf(fp, "        .. \"String\"\n");
    fprintf(fp, "    - arithmetic operators: +,*,-,/\n");
    fprintf(fp, "    - comparison operators: == (same as =), >, >=, <=, <, !=\n");
    fprintf(fp, "    - parentheses: (, )\n");
    fprintf(fp, "    - logical operators: && (same as &), || (same as |)\n");
    fprintf(fp, "    - INFO tags, FORMAT tags, column names\n");
    fprintf(fp, "        .. INFO/DP or DP\n");
    fprintf(fp, "        .. FORMAT/DV, FMT/DV, or DV\n");
    fprintf(fp, "        .. %%FILTER, %%QUAL\n");
    fprintf(fp, "    - 1 (or 0) to test the presence (or absence) of a flag\n");
    fprintf(fp, "        .. FlagA=1 && FlagB=0\n");
    fprintf(fp, "    - %%TYPE for variant type: indel,snp,mnp,other\n");
    fprintf(fp, "        .. %%TYPE=\"indel\" | %%TYPE=\"snp\"\n");
    fprintf(fp, "    - array subscripts\n");
    fprintf(fp, "        .. (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3\n");
    fprintf(fp, "    - operations on FORMAT fields: MAX, MIN, AVG\n");
    fprintf(fp, "        .. %%MIN(DV)>5\n");
    fprintf(fp, "        .. %%MIN(DV/DP)>0.3\n");
    fprintf(fp, "        .. %%MIN(DP)>10 & %%MIN(DV)>3\n");
    fprintf(fp, "        .. FMT/DP>10 & FMT/DV>3   .. the conditions must be valid in a single sample\n");
}

