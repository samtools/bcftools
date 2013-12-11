#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <dlfcn.h>
#include "bcftools.h"

/** Plugin API */
typedef void (*dl_init_f) (bcf_hdr_t *);
typedef void (*dl_process_f) (bcf1_t *);
typedef void (*dl_destroy_f) (void);

typedef struct
{
    char *name;
    dl_init_f init;
    dl_process_f process;
    dl_destroy_f destroy;
    void *handle;
}
plugin_t;

struct _args_t;

typedef struct _rm_tag_t
{
    char *key;
    void (*handler)(struct _args_t *, bcf1_t *, struct _rm_tag_t *);
}
rm_tag_t;

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr, *hdr_out;
    htsFile *out_fh;
    int output_type;
    bcf_sr_regions_t *tgts;

    plugin_t *plugins;
    int nplugins, nplugin_paths;
    char **plugin_paths;

    rm_tag_t *rm;
    int nrm;

    char **argv, *targets_fname, *regions_fname;
    char *remove_annots;
    int argc;
}
args_t;

char *msprintf(const char *fmt, ...);


static void *dlopen_plugin(args_t *args, const char *fname)
{
    if ( args->nplugin_paths==-1 )
    {
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

    tmp = msprintf("%s.so", fname);
    handle = dlopen(tmp, RTLD_NOW);
    free(tmp);
    if ( handle ) return handle;

    return NULL;
}

static void load_plugin(args_t *args, const char *fname)
{
    char *ss, *rmme;
    ss = rmme = strdup(fname);
    while ( ss )
    {
        char *se = strchr(ss,',');
        if ( se ) *se = 0;

        void *handle = dlopen_plugin(args, ss);
        if ( !handle ) error("Could not load %s: %s\n", ss, dlerror());

        args->nplugins++;
        args->plugins = (plugin_t*) realloc(args->plugins, sizeof(plugin_t)*args->nplugins);
        plugin_t *plugin = &args->plugins[args->nplugins-1];
        plugin->handle = handle;
        plugin->name   = strdup(ss);

        dlerror();
        plugin->init = (dl_init_f) dlsym(plugin->handle, "init");
        char *ret = dlerror();
        if ( ret ) error("Could not initialize %s: %s\n", plugin->name, ret);

        plugin->process = (dl_process_f) dlsym(plugin->handle, "process");
        ret = dlerror();
        if ( ret ) error("Could not initialize %s: %s\n", plugin->name, ret);

        plugin->destroy = (dl_destroy_f) dlsym(plugin->handle, "destroy");
        ret = dlerror();
        if ( ret ) error("Could not initialize %s: %s\n", plugin->name, ret);

        if ( se ) { *se = ','; ss = se+1; }
        else ss = NULL;
    }
    free(rmme);
}

static void init_plugins(args_t *args)
{
    int i;
    for (i=0; i<args->nplugins; i++)
        args->plugins[i].init(args->hdr);
}

void remove_id(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    error("todo: -r ID");
}
void remove_filter(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    error("todo: -r FILTER");
}
void remove_qual(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    error("todo: -r QUAL");
}
void remove_info(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    error("todo: -r INFO");
}
void remove_info_tag(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_info(args->hdr, line, tag->key, NULL, 0, BCF_HT_INT);  // the type does not matter with n=0
}
void remove_format_tag(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_format(args->hdr, line, tag->key, NULL, 0, BCF_HT_INT);  // the type does not matter with n=0
}

static void init_data(args_t *args)
{
    args->hdr = args->files->readers[0].header;
    args->hdr_out = bcf_hdr_dup(args->hdr);

    if ( args->remove_annots ) 
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

    bcf_hdr_append_version(args->hdr_out, args->argc, args->argv, "bcftools_annotate");
    args->out_fh = hts_open("-",hts_bcf_wmode(args->output_type));

    if ( args->targets_fname ) args->tgts = bcf_sr_regions_init(args->targets_fname);
    init_plugins(args);
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
}

static void annotate(args_t *args, bcf1_t *line)
{
    int i;
    for (i=0; i<args->nrm; i++)
        args->rm[i].handler(args, line, &args->rm[i]);

    for (i=0; i<args->nplugins; i++)
        args->plugins[i].process(line);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Annotate and edit VCF/BCF files.\n");
    fprintf(stderr, "Usage:   bcftools annotate [OPTIONS] <in.bcf>|<in.vcf>|<in.vcf.gz>|-\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -a, --annotations <file>       tabix-indexed file with annotations: CHR\\tPOS[\\tVALUE]+\n");
	fprintf(stderr, "   -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
	fprintf(stderr, "   -p, --plugins <name|...>       comma-separated list of dynamically loaded user-defined plugins\n");
    fprintf(stderr, "   -r, --regions <reg|file>       restrict to comma-separated list of regions or regions listed in a file, see man page for details\n");
    fprintf(stderr, "   -R, --remove <list>            list of annotations to remove (e.g. ID,INFO/DP,FORMAT/DP,FILTER)\n");
    fprintf(stderr, "\n");
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

    static struct option loptions[] = 
    {
        {"output-type",1,0,'O'},
        {"annotations",1,0,'a'},
        {"plugin",1,0,'p'},
        {"regions",1,0,'r'},
        {"remove",1,0,'R'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "h?O:r:a:p:R:",loptions,NULL)) >= 0) 
    {
        switch (c) {
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
            case 'p': load_plugin(args, optarg); break;
            case 'h': 
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( argc<optind+1 ) usage(args);
    if ( args->regions_fname )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_fname)<0 )
            error("Failed to read the regions: %s\n", args->regions_fname);
    }
    if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open or the file not indexed: %s\n", argv[optind]);
    
    init_data(args);
    bcf_hdr_write(args->out_fh, args->hdr_out);
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = args->files->readers[0].buffer[0];
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
