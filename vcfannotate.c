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
#include "bcftools.h"

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    htsFile *out_fh;
    int output_type;
    bcf_sr_regions_t *tgts;

    char **argv, *targets_fname, *regions_fname;
    int argc;
}
args_t;

static void init_data(args_t *args)
{
    args->hdr = args->files->readers[0].header;
    
    // bcf_hdr_append(args->hdr,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
    bcf_hdr_append_version(args->hdr, args->argc, args->argv, "bcftools_annotate");
    bcf_hdr_fmt_text(args->hdr);

    args->out_fh = hts_open("-",hts_bcf_wmode(args->output_type));

    args->tgts = bcf_sr_regions_init(args->targets_fname);
}

static void destroy_data(args_t *args)
{
    bcf_sr_regions_destroy(args->tgts);
}

static void annotate(args_t *args, bcf1_t *line)
{
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Annotate and edit VCF/BCF files.\n");
    fprintf(stderr, "Usage:   bcftools annotate [OPTIONS] <in.bcf>|<in.vcf>|<in.vcf.gz>|-\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -a, --annotations <file>       tabix-indexed file with annotations: CHR\tPOS[\tVALUE]+\n");
	fprintf(stderr, "    -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "    -r, --regions <reg|file>       same as -t but index-jumps rather than streams to a region (requires indexed VCF/BCF)\n");
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

    error("not finished yet, sorry\n");

    static struct option loptions[] = 
    {
        {"output-type",1,0,'O'},
        {"annotations",1,0,'a'},
        {"regions",1,0,'r'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "h?O:r:a:",loptions,NULL)) >= 0) {
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
            case 'a': args->targets_fname = optarg; break;
            case 'r': args->regions_fname = optarg; break;
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
    if ( !args->targets_fname ) error("Missing the -a option\n");
    if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open or the file not indexed: %s\n", argv[optind]);
    
    init_data(args);
    bcf_hdr_write(args->out_fh, args->hdr);
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = args->files->readers[0].buffer[0];
        annotate(args, line);
        bcf_write1(args->out_fh, args->hdr, line);
    }
    hts_close(args->out_fh);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}
