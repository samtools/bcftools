#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kseq.h>
#include <dlfcn.h>
#include "bcftools.h"
#include "vcmp.h"

/** Plugin API */
typedef void (*dl_init_f) (bcf_hdr_t *);
typedef char* (*dl_about_f) (void);
typedef void (*dl_process_f) (bcf1_t *);
typedef void (*dl_destroy_f) (void);

typedef struct
{
    char *name;
    dl_init_f init;
    dl_about_f about;
    dl_process_f process;
    dl_destroy_f destroy;
    void *handle;
}
plugin_t;

struct _args_t;

typedef struct _rm_tag_t
{
    char *key;
    int hdr_id;
    void (*handler)(struct _args_t *, bcf1_t *, struct _rm_tag_t *);
}
rm_tag_t;

typedef struct
{
    char **cols;
    int ncols, mcols;
    char **als;
    int nals, mals;
    kstring_t line;
    int rid, start, end;
}
annot_line_t;

typedef struct _annot_col_t
{
    int icol;
    char *hdr_key;
    int (*setter)(struct _args_t *, bcf1_t *, struct _annot_col_t *, void*);
}
annot_col_t;

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr, *hdr_out;
    htsFile *out_fh;
    int output_type;
    bcf_sr_regions_t *tgts;

    plugin_t *plugins;      // user plugins
    int nplugins, nplugin_paths;
    char **plugin_paths;

    rm_tag_t *rm;           // tags scheduled for removal
    int nrm;

    vcmp_t *vcmp;           // for matching annotation and VCF lines by allele
    annot_line_t *alines;   // buffered annotation lines
    int nalines, malines;
    int ref_idx, alt_idx, chr_idx, from_idx, to_idx;   // -1 if not present
    annot_col_t *cols;      // column indexes and setters
    int ncols;

    int ntmpi, mtmpi, ntmpf, mtmpf, ntmps, mtmps;
    int32_t *tmpi;
    float *tmpf;
    char *tmps;

    char **argv, *targets_fname, *regions_fname, *header_fname;
    char *remove_annots, *columns;
    int argc;
}
args_t;

char *msprintf(const char *fmt, ...);

static void init_plugin_paths(args_t *args)
{
    if ( args->nplugin_paths!=-1 ) return;

    char *path = getenv("BCFTOOLS_PLUGINS");
    if ( path )
    {
        args->nplugin_paths = 1;
        args->plugin_paths  = (char**) malloc(sizeof(char*));
        char *ss = args->plugin_paths[0] = strdup(path);
        while ( *ss ) 
        { 
            if ( *ss==':' ) 
            {   
                *ss = 0;
                args->plugin_paths = (char**) realloc(args->plugin_paths,sizeof(char*)*(args->nplugin_paths+1));
                args->plugin_paths[args->nplugin_paths] = ss+1;
                args->nplugin_paths++;
            }
            ss++;
        }
    }
    else
        args->nplugin_paths = 0;
}

static void *dlopen_plugin(args_t *args, const char *fname)
{
    init_plugin_paths(args);

    char *tmp;
    void *handle;
    int i;
    for (i=0; i<args->nplugin_paths; i++)
    {
        tmp = msprintf("%s/%s.so", args->plugin_paths[i],fname);
        handle = dlopen(tmp, RTLD_NOW);
        free(tmp);
        if ( handle ) return handle;
    }

    handle = dlopen(fname, RTLD_NOW);
    if ( handle ) return handle;

    return NULL;
}

static int load_plugin(args_t *args, const char *fname, int exit_on_error)
{
    char *ss, *rmme;
    ss = rmme = strdup(fname);
    while ( ss )
    {
        char *se = strchr(ss,',');
        if ( se ) *se = 0;

        void *handle = dlopen_plugin(args, ss);
        if ( !handle )
        {
            if ( exit_on_error ) error("Could not load %s: %s\n", ss, dlerror());
            return -1;
        }

        args->nplugins++;
        args->plugins = (plugin_t*) realloc(args->plugins, sizeof(plugin_t)*args->nplugins);
        plugin_t *plugin = &args->plugins[args->nplugins-1];
        plugin->handle = handle;
        plugin->name   = strdup(ss);

        dlerror();
        plugin->init = (dl_init_f) dlsym(plugin->handle, "init");
        char *ret = dlerror();
        if ( ret ) 
        {
            if ( exit_on_error ) error("Could not initialize %s: %s\n", plugin->name, ret);
            return -1;
        }

        plugin->about = (dl_about_f) dlsym(plugin->handle, "about");
        ret = dlerror();
        if ( ret ) 
        {
            if ( exit_on_error ) error("Could not initialize %s: %s\n", plugin->name, ret);
            return -1;
        }

        plugin->process = (dl_process_f) dlsym(plugin->handle, "process");
        ret = dlerror();
        if ( ret ) 
        {
            if ( exit_on_error ) error("Could not initialize %s: %s\n", plugin->name, ret);
            return -1;
        }

        plugin->destroy = (dl_destroy_f) dlsym(plugin->handle, "destroy");
        ret = dlerror();
        if ( ret ) 
        {
            if ( exit_on_error ) error("Could not initialize %s: %s\n", plugin->name, ret);
            return -1;
        }

        if ( se ) { *se = ','; ss = se+1; }
        else ss = NULL;
    }
    free(rmme);
    return 0;
}

static void init_plugins(args_t *args)
{
    int i;
    for (i=0; i<args->nplugins; i++)
        args->plugins[i].init(args->hdr);
}

static int list_plugins(args_t *args)
{
    init_plugin_paths(args);

    int i, nprinted = 0;
    for (i=0; i<args->nplugin_paths; i++)
    {
        DIR *dp = opendir(args->plugin_paths[i]);
        if ( dp==NULL ) continue;

        struct dirent *ep;
        while ( (ep=readdir(dp)) )
        {
            char *tmp = msprintf("%s/%s", args->plugin_paths[i],ep->d_name);
            int ret = load_plugin(args, tmp, 0);
            free(tmp);
            if ( ret ) continue;
            printf("\n-- %s:\n%s", ep->d_name, args->plugins[args->nplugins-1].about());
            nprinted++;
        }
        closedir(dp);
    }
    if ( !nprinted ) 
        fprintf(stderr,
            "No functional bcftools plugins found:\n"
            " - is the environment variable BCFTOOLS_PLUGINS set?\n"
            " - are shared libraries accesible? (Check with `ldd your/plugin.so`.)\n"
            );
    else
        printf("\n");
    return nprinted ? 0 : 1;
}

void remove_id(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_id(args->hdr,line,NULL);
}
void remove_filter(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    if ( !tag->key ) bcf_update_filter(args->hdr, line, NULL, 0);
    else bcf_remove_filter(args->hdr, line, tag->hdr_id, 1);
}
void remove_qual(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_float_set_missing(line->qual);
}
void remove_info(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    // remove all INFO fields
    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);

    int i;
    for (i=0; i<line->n_info; i++)
        bcf_update_info(args->hdr, line, bcf_hdr_int2id(args->hdr,BCF_DT_ID,line->d.info[i].key), NULL, 0, BCF_HT_INT);  // the type does not matter with n=0
}
void remove_info_tag(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_info(args->hdr, line, tag->key, NULL, 0, BCF_HT_INT);  // the type does not matter with n=0
}
void remove_format_tag(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_format(args->hdr, line, tag->key, NULL, 0, BCF_HT_INT);  // the type does not matter with n=0
}

static void init_remove_annots(args_t *args)
{
    kstring_t str = {0,0,0};
    char *ss = args->remove_annots;
    while ( *ss )
    {
        args->nrm++;
        args->rm = (rm_tag_t*) realloc(args->rm,sizeof(rm_tag_t)*args->nrm);
        rm_tag_t *tag = &args->rm[args->nrm-1];
        tag->key = NULL;

        int type = BCF_HL_GEN;
        if ( !strncmp("INFO/",ss,5) ) { type = BCF_HL_INFO; ss += 5; }
        else if ( !strncmp("INF/",ss,4) ) { type = BCF_HL_INFO; ss += 4; }
        else if ( !strncmp("FORMAT/",ss,7) ) { type = BCF_HL_FMT; ss += 7; }
        else if ( !strncmp("FMT/",ss,4) ) { type = BCF_HL_FMT; ss += 4; }

        char *se = ss;
        while ( *se && *se!=',' ) se++;
        str.l = 0;
        kputsn(ss, se-ss, &str);

        if ( type!= BCF_HL_GEN )
        {
            int id = bcf_hdr_id2int(args->hdr,BCF_DT_ID,str.s);
            if ( !bcf_hdr_idinfo_exists(args->hdr,type,id) ) 
            {
                fprintf(stderr,"Warning: The tag \"%s\" not defined in the header\n", str.s);
                args->nrm--;
            }
            else
            {
                tag->key = strdup(str.s);
                if ( type==BCF_HL_INFO ) tag->handler = remove_info_tag;
                else if ( type==BCF_HL_FMT ) tag->handler = remove_format_tag;
                bcf_hdr_remove(args->hdr_out,type,tag->key);
            }
        }
        else if ( !strcmp("ID",str.s) ) tag->handler = remove_id;
        else if ( !strncmp("FILTER/",str.s,7) ) 
        { 
            tag->handler = remove_filter; 
            tag->key = strdup(str.s+7);
            tag->hdr_id = bcf_hdr_id2int(args->hdr, BCF_DT_ID, tag->key);
            if ( !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FLT,tag->hdr_id) ) error("Cannot remove %s, not defined in the header.\n", str.s);
            bcf_hdr_remove(args->hdr_out,BCF_HL_FLT,tag->key);
        }
        else if ( !strcmp("FILTER",str.s) ) tag->handler = remove_filter;
        else if ( !strcmp("QUAL",str.s) ) tag->handler = remove_qual;
        else if ( !strcmp("INFO",str.s) ) tag->handler = remove_info;
        else if ( str.l )
        {
            if ( str.s[0]=='#' && str.s[1]=='#' )
                bcf_hdr_remove(args->hdr_out,BCF_HL_GEN,str.s+2);
            else
                bcf_hdr_remove(args->hdr_out,BCF_HL_STR,str.s);
            args->nrm--;
        }

        ss = *se ? se+1 : se;
    }
    free(str.s);
    if ( !args->nrm ) error("No matching tag in -R %s\n", args->remove_annots);
}
static void init_header_lines(args_t *args)
{
    htsFile *file = hts_open(args->header_fname, "rb");
    if ( !file ) error("Error reading %s\n", args->header_fname);
    kstring_t str = {0,0,0};
    while ( hts_getline(file, KS_SEP_LINE, &str) > 0 ) 
    {
        if ( bcf_hdr_append(args->hdr_out,str.s) ) error("Could not parse %s: %s\n", args->header_fname, str.s);
    }
    hts_close(file);
    free(str.s);
}
static int setter_id(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    return bcf_update_id(args->hdr_out,line,tab->cols[col->icol]);
}
static int vcf_setter_id(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    return bcf_update_id(args->hdr_out,line,rec->d.id);
}
static int setter_qual(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) return 0;

    line->qual = strtod(str, &str);
    if ( str == tab->cols[col->icol] ) 
        error("Could not parse %s at %s:%d .. [%s]\n", col->hdr_key,bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);
    return 0;
}
static int vcf_setter_qual(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    line->qual = rec->qual;
    return 0;
}
static int setter_info_flag(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) return 0;

    if ( str[0]=='1' && str[1]==0 ) return bcf_update_info_flag(args->hdr_out,line,col->hdr_key,NULL,1);
    if ( str[0]=='0' && str[1]==0 ) return bcf_update_info_flag(args->hdr_out,line,col->hdr_key,NULL,0);
    error("Could not parse %s at %s:%d .. [%s]\n", bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);
    return -1;
}
static int vcf_setter_info_flag(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int flag = bcf_get_info_flag(args->files->readers[1].header,rec,col->hdr_key,NULL,NULL);
    bcf_update_info_flag(args->hdr_out,line,col->hdr_key,NULL,flag);
    return 0;
}
static int setter_info_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data) 
{ 
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol], *end = str;
    if ( str[0]=='.' && str[1]==0 ) return 0;

    args->ntmpi = 0;
    while ( *end )
    {
        int val = strtol(str, &end, 10);
        if ( end==str )
            error("Could not parse %s at %s:%d .. [%s]\n", bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);
        args->ntmpi++;
        hts_expand(int,args->ntmpi,args->mtmpi,args->tmpi);
        args->tmpi[args->ntmpi-1] = val;
        str = end+1;
    }
    return bcf_update_info_int32(args->hdr_out,line,col->hdr_key,args->tmpi,args->ntmpi);
}
static int vcf_setter_info_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    args->ntmpi = bcf_get_info_int32(args->files->readers[1].header,rec,col->hdr_key,&args->tmpi,&args->mtmpi);
    if ( args->ntmpi >=0 )
        bcf_update_info_int32(args->hdr_out,line,col->hdr_key,args->tmpi,args->ntmpi);
    return 0;
}
static int setter_info_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{ 
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol], *end = str;
    if ( str[0]=='.' && str[1]==0 ) return 0;

    args->ntmpf = 0;
    while ( *end )
    {
        double val = strtod(str, &end);
        if ( end==str )
            error("Could not parse %s at %s:%d .. [%s]\n", bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);
        args->ntmpf++;
        hts_expand(float,args->ntmpf,args->mtmpf,args->tmpf);
        args->tmpf[args->ntmpf-1] = val;
        str = end+1;
    }
    return bcf_update_info_float(args->hdr_out,line,col->hdr_key,args->tmpf,args->ntmpf);
}
static int vcf_setter_info_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    args->ntmpf = bcf_get_info_float(args->files->readers[1].header,rec,col->hdr_key,&args->tmpf,&args->mtmpf);
    if ( args->ntmpf >=0 )
        bcf_update_info_float(args->hdr_out,line,col->hdr_key,args->tmpf,args->ntmpf);
    return 0;
}
static int setter_info_str(args_t *args, bcf1_t *line, annot_col_t *col, void *data) 
{
    annot_line_t *tab = (annot_line_t*) data;
    return bcf_update_info_string(args->hdr_out,line,col->hdr_key,tab->cols[col->icol]);
}
static int vcf_setter_info_str(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    args->ntmps = bcf_get_info_string(args->files->readers[1].header,rec,col->hdr_key,&args->tmps,&args->mtmps);
    if ( args->ntmps >=0 )
        bcf_update_info_string(args->hdr_out,line,col->hdr_key,args->tmps);
    return 0;
}
static void init_columns(args_t *args)
{
    kstring_t str = {0,0,0}, tmp = {0,0,0};
    char *ss = args->columns, *se = ss;
    args->ncols = 0;
    int i = -1;
    se = ss;
    while ( *ss )
    {
        if ( *se && *se!=',' ) { se++; continue; }
        i++;
        str.l = 0;
        kputsn(ss, se-ss, &str);
        if ( !strcasecmp("-",str.s) ) ;
        else if ( !strcasecmp("CHROM",str.s) ) args->chr_idx = i;
        else if ( !strcasecmp("POS",str.s) ) args->from_idx = i;
        else if ( !strcasecmp("FROM",str.s) ) args->from_idx = i;
        else if ( !strcasecmp("TO",str.s) ) args->to_idx = i;
        else if ( !strcasecmp("REF",str.s) ) args->ref_idx = i;
        else if ( !strcasecmp("ALT",str.s) ) args->alt_idx = i;
        else if ( !strcasecmp("ID",str.s) )
        {
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = i;
            col->setter = args->tgts ? setter_id : vcf_setter_id;
            col->hdr_key = strdup(str.s);
        }
        else if ( !strcasecmp("QUAL",str.s) )
        {
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = i;
            col->setter = args->tgts ? setter_qual : vcf_setter_qual;
            col->hdr_key = strdup(str.s);
        }
        else 
        {
            if ( !strncasecmp("INFO/",str.s,5) ) { memmove(str.s,str.s+5,str.l-4); }
            int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, str.s);
            if ( !bcf_hdr_idinfo_exists(args->hdr_out,BCF_HL_INFO,hdr_id) )
            {
                if ( args->files->require_index ) // reading annotations from a VCF, add a new header line
                {
                    bcf_hrec_t *hrec = bcf_hdr_get_hrec(args->files->readers[1].header, BCF_HL_INFO, str.s);
                    if ( !hrec ) error("The tag \"%s\" is not defined in %s\n", str.s,args->files->readers[1].fname);
                    tmp.l = 0;
                    bcf_hrec_format(hrec, &tmp);
                    bcf_hdr_append(args->hdr_out, tmp.s);
                    hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, str.s);
                }
                if ( !bcf_hdr_idinfo_exists(args->hdr_out,BCF_HL_INFO,hdr_id) )
                    error("The column not recognised: [%s] .. %d\n", str.s, hdr_id);
            }

            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = i;
            col->hdr_key = strdup(str.s);
            switch ( bcf_hdr_id2type(args->hdr_out,BCF_HL_INFO,hdr_id) )
            {
                case BCF_HT_FLAG:   col->setter = args->tgts ? setter_info_flag : vcf_setter_info_flag; break;
                case BCF_HT_INT:    col->setter = args->tgts ? setter_info_int  : vcf_setter_info_int; break;
                case BCF_HT_REAL:   col->setter = args->tgts ? setter_info_real : vcf_setter_info_real; break;
                case BCF_HT_STR:    col->setter = args->tgts ? setter_info_str  : vcf_setter_info_str; break;
                default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(args->hdr_out,BCF_HL_INFO,hdr_id));
            }
        }
        if ( !*se ) break;
        ss = ++se;
    }
    free(str.s);
    free(tmp.s);
    if ( args->to_idx==-1 ) args->to_idx = args->from_idx;
}

static void init_data(args_t *args)
{
    args->hdr = args->files->readers[0].header;
    args->hdr_out = bcf_hdr_dup(args->hdr);

    if ( args->targets_fname )
    {
        if ( args->files->require_index )   // reading annots from a VCF
        {
            if ( !bcf_sr_add_reader(args->files, args->targets_fname) )
                error("Failed to open or the file not indexed: %s\n", args->targets_fname);
        }
        else
        {
            args->tgts = bcf_sr_regions_init(args->targets_fname,args->chr_idx,args->from_idx,args->to_idx);
            if ( !args->tgts ) error("Could not initialize the annotation file: %s\n", args->targets_fname);
            if ( !args->tgts->tbx ) error("Expected tabix-indexed annotation file: %s\n", args->targets_fname);
            args->vcmp = vcmp_init();
        }
    }
    if ( args->remove_annots ) init_remove_annots(args);
    if ( args->header_fname ) init_header_lines(args);
    if ( args->columns ) init_columns(args);
    init_plugins(args);

    bcf_hdr_append_version(args->hdr_out, args->argc, args->argv, "bcftools_annotate");
    args->out_fh = hts_open("-",hts_bcf_wmode(args->output_type));
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->nplugins; i++)
    {
        free(args->plugins[i].name);
        args->plugins[i].destroy();
        dlclose(args->plugins[i].handle);
    }
    free(args->plugins);
    for (i=0; i<args->nrm; i++) free(args->rm[i].key);
    free(args->rm);
    bcf_hdr_destroy(args->hdr_out);
    if ( args->nplugin_paths>0 )
    {   
        free(args->plugin_paths[0]);
        free(args->plugin_paths);
    }
    if (args->vcmp) vcmp_destroy(args->vcmp);
    for (i=0; i<args->ncols; i++)
        free(args->cols[i].hdr_key);
    free(args->cols);
    for (i=0; i<args->malines; i++)
    {
        free(args->alines[i].cols);
        free(args->alines[i].als);
        free(args->alines[i].line.s);
    }
    free(args->alines);
    if ( args->tgts ) bcf_sr_regions_destroy(args->tgts);
    free(args->tmpi);
    free(args->tmpf);
    free(args->tmps);
}

static void buffer_annot_lines(args_t *args, bcf1_t *line, int start_pos, int end_pos)
{
    if ( args->nalines && args->alines[0].rid != line->rid ) args->nalines = 0;

    int i = 0;
    while ( i<args->nalines )
    {
        if ( line->pos > args->alines[i].end )
        {
            args->nalines--;
            if ( args->nalines && i<args->nalines )
            {
                annot_line_t tmp = args->alines[i];
                memmove(&args->alines[i],&args->alines[i+1],(args->nalines-i)*sizeof(annot_line_t));
                args->alines[args->nalines] = tmp;
            }
        }
        else i++;
    }

    if ( args->ref_idx==-1 && args->nalines ) return;

    while ( !bcf_sr_regions_overlap(args->tgts, bcf_seqname(args->hdr,line), start_pos,end_pos) )
    {
        args->nalines++;
        hts_expand0(annot_line_t,args->nalines,args->malines,args->alines);
        annot_line_t *tmp = &args->alines[args->nalines-1];
        tmp->rid   = line->rid;
        tmp->start = args->tgts->start;
        tmp->end   = args->tgts->end;
        tmp->line.l = 0;
        kputs(args->tgts->line.s, &tmp->line);
        char *s = tmp->line.s;
        tmp->ncols = 1;
        hts_expand(char*,tmp->ncols,tmp->mcols,tmp->cols);
        tmp->cols[0] = s;
        while ( *s )
        {
            if ( *s=='\t' )
            {
                tmp->ncols++;
                hts_expand(char*,tmp->ncols,tmp->mcols,tmp->cols);
                tmp->cols[tmp->ncols-1] = s+1;
                *s = 0;
            }
            s++;
        }
        if ( args->ref_idx != -1 )
        {
            assert( args->ref_idx < tmp->ncols );
            assert( args->alt_idx < tmp->ncols );
            s = tmp->cols[args->alt_idx];
            tmp->nals = 1;
            hts_expand(char*,tmp->nals,tmp->mals,tmp->als);
            tmp->als[0] = s;
            while ( *s )
            {
                if ( *s==',' )
                {
                    tmp->nals++;
                    hts_expand(char*,tmp->nals,tmp->mals,tmp->als);
                    tmp->als[tmp->nals-1] = s;
                    *s = 0;
                }
                s++;
            }
            int iseq = args->tgts->iseq;
            if ( bcf_sr_regions_next(args->tgts)<0 || args->tgts->iseq!=iseq ) break;
        }
        else break;
    }
}

static void annotate(args_t *args, bcf1_t *line)
{
    int i, j;
    for (i=0; i<args->nrm; i++)
        args->rm[i].handler(args, line, &args->rm[i]);

    if ( args->tgts )
    {
        // Buffer annotation lines. When multiple ALT alleles are present in the
        // annotation file, at least one must match one of the VCF alleles.
        int len = 0;
        bcf_get_variant_types(line);
        for (i=1; i<line->n_allele; i++)
            if ( len > line->d.var[i].n ) len = line->d.var[i].n;
        int end_pos = len<0 ? line->pos - len : line->pos;
        buffer_annot_lines(args, line, line->pos, end_pos);
        for (i=0; i<args->nalines; i++)
        {
            if ( line->pos > args->alines[i].end || end_pos < args->alines[i].start ) continue;
            if ( args->ref_idx != -1 )
            {
                if ( vcmp_set_ref(args->vcmp, line->d.allele[0], args->alines[i].cols[args->ref_idx]) < 0 ) continue;   // refs not compatible
                for (j=0; j<args->alines[i].nals; j++)
                {
                    if ( line->n_allele==1 && args->alines[i].als[j][0]=='.' && args->alines[i].als[j][1]==0 ) break;   // no ALT allele in VCF and annot file has "."
                    if ( vcmp_find_allele(args->vcmp, line->d.allele+1, line->n_allele - 1, args->alines[i].als[j]) >= 0 ) break;
                }
                if ( j==args->alines[i].nals ) continue;    // none of the annot alleles present in VCF's ALT
            }
            break;
        }

        if ( i<args->nalines )
        {
            for (j=0; j<args->ncols; j++)
                if ( args->cols[j].setter(args,line,&args->cols[j],&args->alines[i]) ) 
                    error("fixme: Could not set %s at %s:%d\n", args->cols[j].hdr_key,bcf_seqname(args->hdr,line),line->pos+1);
        }
    }
    else if ( args->files->nreaders == 2 && bcf_sr_has_line(args->files,1) )
    {
        bcf1_t *aline = bcf_sr_get_line(args->files,1);
        for (j=0; j<args->ncols; j++)
            if ( args->cols[j].setter(args,line,&args->cols[j],aline) ) 
                error("fixme: Could not set %s at %s:%d\n", args->cols[j].hdr_key,bcf_seqname(args->hdr,line),line->pos+1);
    }

    for (i=0; i<args->nplugins; i++)
        args->plugins[i].process(line);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Annotate and edit VCF/BCF files.\n");
    fprintf(stderr, "Usage:   bcftools annotate [options] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -a, --annotations <file>       VCF file or tabix-indexed file with annotations: CHR\\tPOS[\\tVALUE]+\n");
    fprintf(stderr, "   -c, --columns <list>           list of columns in the annotation file, e.g. CHROM,POS,REF,ALT,-,INFO/TAG. See man page for details\n");
    fprintf(stderr, "   -h, --header-lines <file>      lines which should to appended to the VCF header\n");
	fprintf(stderr, "   -l, --list-plugins             list available plugins. See BCFTOOLS_PLUGINS environment variable and man page for details\n");
	fprintf(stderr, "   -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
	fprintf(stderr, "   -p, --plugins <name|...>       comma-separated list of dynamically loaded user-defined plugins. See man page for details\n");
    fprintf(stderr, "   -r, --regions <reg|file>       restrict to comma-separated list of regions or regions listed in a file, see man page for details\n");
    fprintf(stderr, "   -R, --remove <list>            list of annotations to remove (e.g. ID,INFO/DP,FORMAT/DP,FILTER). See man page for details\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfannotate(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->output_type = FT_VCF;
    args->nplugin_paths = -1;
    args->ref_idx = args->alt_idx = args->chr_idx = args->from_idx = args->to_idx = -1;
    int plist_only = 0;

    static struct option loptions[] = 
    {
        {"output-type",1,0,'O'},
        {"annotations",1,0,'a'},
        {"plugin",1,0,'p'},
        {"list-plugins",0,0,'l'},
        {"regions",1,0,'r'},
        {"remove",1,0,'R'},
        {"columns",1,0,'c'},
        {"header-lines",1,0,'h'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "h:?O:r:a:p:R:c:l",loptions,NULL)) >= 0) 
    {
        switch (c) {
    	    case 'c': args->columns = optarg; break;
    	    case 'O': 
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                };
                break;
            case 'R': args->remove_annots = optarg; break;
            case 'a': args->targets_fname = optarg; break;
            case 'r': args->regions_fname = optarg; break;
            case 'p': load_plugin(args, optarg, 1); break;
            case 'l': plist_only = 1; break;
            case 'h': args->header_fname = optarg; break;
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( plist_only ) return list_plugins(args);

    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else usage(args);
    }
    else fname = argv[optind];

    if ( args->regions_fname )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_fname)<0 )
            error("Failed to read the regions: %s\n", args->regions_fname);
    }
    if ( args->targets_fname && hts_file_type(args->targets_fname) & (FT_VCF|FT_BCF) ) args->files->require_index = 1;
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open or the file not indexed: %s\n", fname);
    
    init_data(args);
    bcf_hdr_write(args->out_fh, args->hdr_out);
    while ( bcf_sr_next_line(args->files) )
    {
        if ( !bcf_sr_has_line(args->files,0) ) continue;
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( line->errcode ) error("Encountered error, cannot proceed. Please check the error output above.\n");
        annotate(args, line);
        bcf_write1(args->out_fh, args->hdr_out, line);
    }
    hts_close(args->out_fh);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}
