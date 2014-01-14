#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kseq.h>
#include "bcftools.h"

typedef struct _args_t
{
    bcf_srs_t *files;
    htsFile *out_fh;
    int output_type;
    bcf_hdr_t *hdr_out;

    char **argv, *file_list;
    int argc;
}
args_t;

static void init_data(args_t *args)
{
#if 0
    int i;
    bcf_sr_regions_t *regs = args->files->regions;
    args->hdr_out = bcf_hdr_dup(args->files->readers[0].header);
    for (i=0; i<regs->nseqs; i++)
    {
        
    }
    bcf_hdr_append_version(args->hdr_out, args->argc, args->argv, "bcftools_concat");
    args->out_fh = hts_open("-",hts_bcf_wmode(args->output_type));
#endif
}

static void destroy_data(args_t *args)
{
}

#if 0
static void concat(args_t *args, bcf1_t *line)
{
}
#endif

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Concatenate and combine VCF/BCF files split e.g. by chromosome or variant type.\n");
    fprintf(stderr, "         All source files must have the same samples in the same order.\n");
    fprintf(stderr, "Usage:   bcftools concat [OPTIONS] <file> [...]\n");
    fprintf(stderr, "Options:\n");
	fprintf(stderr, "   -a, --allow-overlaps           The files may overlap, require.\n");
	fprintf(stderr, "   -l, --file-list <file>         Read the list of files from a file.\n");
	fprintf(stderr, "   -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfconcat(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->output_type = FT_VCF;

    error("not ready yet, sorry\n");

    static struct option loptions[] = 
    {
        {"output-type",1,0,'O'},
        {"file-list",1,0,'l'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "h:?O:l:",loptions,NULL)) >= 0) 
    {
        switch (c) {
    	    case 'l': args->file_list = optarg; break;
    	    case 'O': 
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                };
                break;
            case 'h': 
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( args->file_list ) error("todo\n");
    else if ( optind==argc ) usage(args);

    while ( optind<argc )
    {
        if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open or the file not indexed: %s\n", argv[optind]);
        optind++;
    }
    
    init_data(args);
    bcf_hdr_write(args->out_fh, args->hdr_out);
    while ( bcf_sr_next_line(args->files) )
    {
    }
    hts_close(args->out_fh);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}
