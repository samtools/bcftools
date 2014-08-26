/*  vcfconvert.c -- convert between VCF/BCF and related formats.

    Copyright (C) 2013-2014 Genome Research Ltd.

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
THE SOFTWARE.  */

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "filter.h"
#include "convert.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct _args_t args_t;
struct _args_t
{
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE
    convert_t *convert;
    bcf_srs_t *files;
    bcf_hdr_t *header;
    void (*convert_func)(struct _args_t *);
    int nsamples, *samples, sample_is_file, targets_is_file, regions_is_file;
    char **argv, *sample_list, *targets_list, *regions_list, *tag;
    char *outfname, *infname;
    int argc;
};

static void destroy_data(args_t *args)
{
    if ( args->convert) convert_destroy(args->convert);
    if ( args->filter ) filter_destroy(args->filter);
    free(args->samples);
    if ( args->files ) bcf_sr_destroy(args->files);
}

static void open_vcf(args_t *args, const char *format_str)
{
    args->files = bcf_sr_init();
    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_list, args->regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_list, args->targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }
    if ( !bcf_sr_add_reader(args->files, args->infname) )
        error("Failed to open or the file not indexed: %s\n", args->infname);
    if ( args->filter_str )
        args->filter = filter_init(args->header, args->filter_str);

    args->header = args->files->readers[0].header;

    int i, nsamples = 0, *samples = NULL;
    if ( args->sample_list && strcmp("-",args->sample_list) )
    {
        for (i=0; i<args->files->nreaders; i++)
        {
            int ret = bcf_hdr_set_samples(args->files->readers[i].header,args->sample_list,args->sample_is_file);
            if ( ret<0 ) error("Error parsing the sample list\n");
            else if ( ret>0 ) error("Sample name mismatch: sample #%d not found in the header\n", ret);
        }

        if ( args->sample_list[0]!='^' )
        {
            // the sample ordering may be different if not negated
            int n;
            char **smpls = hts_readlist(args->sample_list, args->sample_is_file, &n);
            if ( !smpls ) error("Could not parse %s\n", args->sample_list);
            if ( n!=bcf_hdr_nsamples(args->files->readers[0].header) )
                error("The number of samples does not match, perhaps some are present multiple times?\n");
            nsamples = bcf_hdr_nsamples(args->files->readers[0].header);
            samples = (int*) malloc(sizeof(int)*nsamples);
            for (i=0; i<n; i++)
            {
                samples[i] = bcf_hdr_id2int(args->files->readers[0].header, BCF_DT_SAMPLE,smpls[i]);
                free(smpls[i]);
            }
            free(smpls);
        }
    }
    args->convert = convert_init(args->header, samples, nsamples, format_str);
    free(samples);

    if ( args->filter_str )
        args->filter = filter_init(args->header, args->filter_str);
}

static void vcf_to_gensample(args_t *args)
{
    kstring_t str = {0,0,0};
    kputs("%CHROM:%POS\\_%REF\\_%FIRST_ALT %_CHROM_POS_ID %POS %REF %FIRST_ALT", &str);
    if ( !args->tag || !strcmp(args->tag,"GT") ) kputs("%_GT_TO_PROB3",&str);
    else if ( !strcmp(args->tag,"PL") ) kputs("%_PL_TO_PROB3",&str);
    else error("todo: --tag %s\n", args->tag);
    kputs("\n", &str);
    open_vcf(args,str.s);

    int is_compressed = 1;
    char *gen_fname = NULL, *sample_fname = NULL;
    str.l = 0;
    kputs(args->outfname,&str);
    char *t = str.s;
    while ( *t && *t!=',' ) t++;
    if ( *t )
    {
        *t = 0;
        gen_fname = strdup(str.s);
        if ( strlen(gen_fname)<3 || strcasecmp(".gz",gen_fname+strlen(gen_fname)-3) ) is_compressed = 0;
        sample_fname = strdup(t+1);
    }
    else
    {
        int l = str.l;
        kputs(".samples",&str);
        sample_fname = strdup(str.s);
        str.l = l;
        kputs(".gen.gz",&str);
        gen_fname = strdup(str.s);
    }

    int i;
    FILE *fh = fopen(sample_fname,"w");
    if ( !fh ) error("Failed to write %s: %s\n", sample_fname,strerror(errno));
    fprintf(fh,"ID_1 ID_2 missing\n0 0 0\n");
    for (i=0; i<bcf_hdr_nsamples(args->header); i++)
        fprintf(fh,"%s %s 0\n", args->header->samples[i],args->header->samples[i]);
    fclose(fh);

    BGZF *out = bgzf_open(gen_fname, is_compressed ? "w" : "wu");
    free(sample_fname);

    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( line->n_allele==1 ) continue;  // alternate allele is required
        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) continue;
        }

        str.l = 0;
        convert_line(args->convert, line, &str);
        if ( str.l )
        {
            int ret = bgzf_write(out, str.s, str.l);
            if ( ret!= str.l ) error("Error writing %s: %s\n", gen_fname,strerror(errno));
        }
    }
    if ( str.m ) free(str.s);
    if ( bgzf_close(out)!=0 ) error("Error closing %s: %s\n", gen_fname,strerror(errno));
    free(gen_fname);
}

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Converts VCF/BCF to other formats and back. See man page for file\n");
    fprintf(stderr, "         formats details\n");
    fprintf(stderr, "Usage:   bcftools convert [OPTIONS] <input_file>\n");
    fprintf(stderr, "VCF input options:\n");
    fprintf(stderr, "   -e, --exclude <expr>        exclude sites for which the expression is true\n");
    fprintf(stderr, "   -i, --include <expr>        select sites for which the expression is true\n");
    fprintf(stderr, "   -r, --regions <region>      restrict to comma-separated list of regions\n");
    fprintf(stderr, "   -R, --regions-file <file>   restrict to regions listed in a file\n");
    fprintf(stderr, "   -s, --samples <list>        list of samples to include\n");
    fprintf(stderr, "   -S, --samples-file <file>   file of samples to include\n");
    fprintf(stderr, "   -t, --targets <region>      similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "   -T, --targets-file <file>   similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "gen/sample options:\n");
//  fprintf(stderr, "   -G, --gensample2vcf  <prefix> or <gen-file>,<sample-file>\n");
    fprintf(stderr, "   -g, --gensample      <prefix> or <gen-file>,<sample-file>\n");
    fprintf(stderr, "   --tag <string>       tag to take values for .gen file: GT,PL,GL,GP [GT]\n");
    fprintf(stderr, "\n");
//    fprintf(stderr, "hap/legend/sample options:\n");
//    fprintf(stderr, "   -h, --haplegend     <prefix> or <hap-file>,<legend-file>,<sample-file>\n");
//    fprintf(stderr, "\n");
//    fprintf(stderr, "plink options:\n");
//    fprintf(stderr, "   -p, --plink         <prefix> or <ped>,<map>,<fam>/<bed>,<bim>,<fam>/<tped>,<tfam>\n");
//    fprintf(stderr, "   --tped              make tped file instead\n");
//    fprintf(stderr, "   --bin               make binary bed/fam/bim files\n");
//    fprintf(stderr, "\n");
//    fprintf(stderr, "pbwt options:\n");
//    fprintf(stderr, "   -b, --pbwt          <prefix> or <pbwt>,<sites>,<sample>,<missing>\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfconvert(int argc, char *argv[])
{
    int c;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;

    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"include",1,0,'i'},
        {"exclude",1,0,'e'},
        {"output-file",1,0,'o'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"gensample",1,0,'g'},
        {"tag",1,0,1},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "?hr:R:s:S:t:T:i:e:g:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; args->regions_is_file = 1; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; args->targets_is_file = 1; break;
            case 's': args->sample_list = optarg; break;
            case 'S': args->sample_list = optarg; args->sample_is_file = 1; break;
            case 'g': args->convert_func = vcf_to_gensample; args->outfname = optarg; break;
            case  1 : args->tag = optarg; break;
            case 'h':
            case '?': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    args->infname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->infname = "-";
    }
    else args->infname = argv[optind];
    if ( !args->infname || !args->convert_func ) usage();

    args->convert_func(args);

    destroy_data(args);
    free(args);
    return 0;
}


