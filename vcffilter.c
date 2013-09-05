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
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "filter.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

// FILTER columns annotation: replace or add to existing FILTERs; set FILTER to PASS at good sites?
#define ANNOT_ADD   1
#define ANNOT_RESET 2

typedef struct _args_t
{
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE
    char *soft_filter;  // drop failed sites or annotate FILTER column?
    int annot_mode;     // add to existing FILTER annotation or replace? Otherwise reset FILTER to PASS or leave as it is?
    int flt_fail, flt_pass;     // BCF ids of fail and pass filters

    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    htsFile *out_fh;
    int output_type;

    char **argv, *targets_fname, *regions_fname;
    int argc;
}
args_t;

static void init_data(args_t *args)
{
    args->out_fh = hts_open("-",hts_bcf_wmode(args->output_type),0);

    args->hdr = args->files->readers[0].header;
    if ( args->soft_filter )
    {
        kstring_t flt_name = {0,0,0}, hdr_line = {0,0,0};
        if ( strcmp(args->soft_filter,"+") )
            kputs(args->soft_filter, &flt_name);
        else
        {
            // Make up a filter name
            int i = 0, id = -1;
            do
            {
                ksprintf(&flt_name,"Filter%d", ++i);
                id = bcf_id2int(args->hdr,BCF_DT_ID,flt_name.s);
            }
            while ( id>=0 && bcf_idinfo_exists(args->hdr,BCF_HL_FLT,id) );
        }
        ksprintf(&hdr_line,"##FILTER=<ID=%s,Description=\"Set if %s: %s\">", flt_name.s,args->filter_logic & FLT_INCLUDE ? "not true" : "true", args->filter_str);
        bcf_hdr_append(args->hdr, hdr_line.s);

        args->flt_pass = bcf_id2int(args->hdr,BCF_DT_ID,"PASS"); assert( !args->flt_pass );  // sanity check: required by BCF spec
        args->flt_fail = bcf_id2int(args->hdr,BCF_DT_ID,flt_name.s); assert( args->flt_fail>=0 );
        free(flt_name.s);
        free(hdr_line.s);
    }
    bcf_hdr_append_version(args->hdr, args->argc, args->argv, "bcftools_filter");
    bcf_hdr_fmt_text(args->hdr);

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);
}

static void destroy_data(args_t *args)
{
    if ( args->filter )
        filter_destroy(args->filter);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Apply fixed-threshold filters.\n");
    fprintf(stderr, "Usage:   bcftools filter [options] <in.bcf>|<in.vcf>|<in.vcf.gz> [region1 [...]]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -e, --exclude <expr>          exclude sites for which the expression is true (e.g. '%%TYPE=snp && %%QUAL>=10 && (DP4[2]+DP4[3] > 2')\n");
    fprintf(stderr, "    -i, --include <expr>          include only sites for which the expression is true\n");
    fprintf(stderr, "    -m, --mode <+|x>              \"+\": do not replace but add to existing FILTER; \"x\": reset filters at sites which pass (invokes -s)\n");
    fprintf(stderr, "    -o, --output-type <b|u|z|v>   b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "    -r, --regions <reg|file>      same as -t but index-jumps rather than streams to a region (requires indexed VCF/BCF)\n");
    fprintf(stderr, "    -s, --soft-filter <string>    annotate FILTER column with <string> or unique filter name (\"Filter%%d\") made up by the program (\"+\")\n");
    fprintf(stderr, "    -t, --targets <reg|file>      restrict to positions in tab-delimited tabix indexed file <chr,pos> or <chr,from,to>, 1-based, inclusive\n");
    fprintf(stderr, "Expressions:\n");
    fprintf(stderr, "    - arithmetic perators: +,*,-,/\n");
    fprintf(stderr, "    - logical operators: && (same as &), || (same as |)\n");
    fprintf(stderr, "    - conditional operators: == (same as =), >, >=, <=, <, !=\n");
    fprintf(stderr, "    - parentheses: (, )\n");
    fprintf(stderr, "    - array subscripts, such as AC[0]\n");
    fprintf(stderr, "    - double quotes for string values\n");
    fprintf(stderr, "    - %%QUAL or INFO tag names\n");
    fprintf(stderr, "    - %%TYPE for variant type, such as %%TYPE=indel|snp|mnp|other\n");
    exit(1);
}

int main_vcffilter(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->output_type = FT_VCF;

    static struct option loptions[] = 
    {
        {"mode",1,0,'m'},
        {"soft-filter",1,0,'s'},
        {"exclude",1,0,'e'},
        {"include",1,0,'i'},
        {"targets",1,0,'t'},
        {"regions",1,0,'r'},
        {"output-type",1,0,'o'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "e:i:t:r:h?s:m:o:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'o': 
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                }
                break;
            case 's': args->soft_filter = optarg; break;
            case 'm': 
                if ( strchr(optarg,'x') ) args->annot_mode |= ANNOT_RESET; 
                if ( strchr(optarg,'+') ) args->annot_mode |= ANNOT_ADD; 
                if ( !args->soft_filter ) args->soft_filter = "+";
                break;
            case 't': args->targets_fname = optarg; break;
            case 'r': args->regions_fname = optarg; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'h':
            case '?': usage(args);
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");
    if ( argc<optind+1 ) usage(args);
    // read in the regions from the command line
    if ( args->regions_fname )
    {
        args->files->require_index = 1;
        if ( bcf_sr_set_regions(args->files, args->regions_fname)<0 )
            error("Failed to read the regions: %s\n", args->regions_fname);
    }
    else if ( optind+1 < argc ) 
    {
        int i;
        kstring_t tmp = {0,0,0};
        kputs(argv[optind+1],&tmp);
        for (i=optind+2; i<argc; i++) { kputc(',',&tmp); kputs(argv[i],&tmp); }
        args->files->require_index = 1;
        if ( bcf_sr_set_regions(args->files, args->regions_fname)<0 )
            error("Failed to read the regions: %s\n", args->regions_fname);
    }
    if ( args->targets_fname )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_fname,0)<0 )
            error("Failed to read the targets: %s\n", args->targets_fname);
    }
    if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open or the file not indexed: %s\n", argv[optind]);
    
    init_data(args);
    vcf_hdr_write(args->out_fh, args->hdr);
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files, 0);
        int pass = filter_test(args->filter, line);
        if ( args->filter_logic & FLT_EXCLUDE ) 
            pass = pass ? 0 : 1;
        if ( args->soft_filter || pass )
        {
            if ( pass ) 
            {
                if ( (args->annot_mode & ANNOT_RESET) ) bcf1_add_filter(args->hdr, line, args->flt_pass);
            }
            else
            {
                if ( (args->annot_mode & ANNOT_ADD) ) bcf1_add_filter(args->hdr, line, args->flt_fail);
                else bcf1_update_filter(args->hdr, line, &args->flt_fail, 1);
            }
            if ( args->output_type & FT_BCF ) bcf1_sync(line);
            vcf_write1(args->out_fh, args->hdr, line);
        }
    }
    hts_close(args->out_fh);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}
