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
    bcf_hdr_t *out_hdr;
    int *seen_seq;

    char **argv, *file_list, **fnames;
    int argc, nfnames, allow_overlaps, phased_concat;
}
args_t;

static void init_data(args_t *args)
{
    kstring_t str = {0,0,0};
    int i;
    for (i=0; i<args->nfnames; i++)
    {
        htsFile *fp = hts_open(args->fnames[i], "r"); if ( !fp ) error("Failed to open: %s\n", args->fnames[i]);
        bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) error("Failed to parse header: %s\n", args->fnames[i]);
        hts_close(fp);
        if ( !args->out_hdr ) 
            args->out_hdr = bcf_hdr_dup(hdr);
        else
        {
            int nseqs, j;
            const char **seqs = bcf_hdr_seqnames(hdr,&nseqs);
            for (j=0; j<nseqs; j++)
            {
                int id;
                if ( (id=bcf_hdr_id2int(args->out_hdr,BCF_DT_CTG,seqs[j])) >=0  ) continue;
                id = bcf_hdr_id2int(hdr,BCF_DT_CTG,seqs[j]);
                bcf_hrec_t *hrec = bcf_hdr_id2hrec(hdr,BCF_DT_CTG,id);
                str.l = 0;
                bcf_hrec_format(hrec, &str);
                bcf_hdr_append(args->out_hdr, str.s);
            }
            free(seqs);

            if ( bcf_hdr_nsamples(hdr) != bcf_hdr_nsamples(args->out_hdr) ) 
                error("Different number of samples in %s. Perhaps \"bcftools merge\" is what you are looking for?\n", args->fnames[i]);

            for (j=0; j<bcf_hdr_nsamples(hdr); j++)
                if ( strcmp(args->out_hdr->samples[j],hdr->samples[j]) )
                    error("Different sample names in %s. Perhaps \"bcftools merge\" is what you are looking for?\n", args->fnames[i]);
        }
        bcf_hdr_destroy(hdr);
    }
    free(str.s);

    args->seen_seq = (int*) calloc(args->out_hdr->n[BCF_DT_CTG],sizeof(int));

    bcf_hdr_append_version(args->out_hdr, args->argc, args->argv, "bcftools_concat");
    args->out_fh = hts_open("-",hts_bcf_wmode(args->output_type));
    bcf_hdr_write(args->out_fh, args->out_hdr);

    if ( args->allow_overlaps ) 
    {
        args->files = bcf_sr_init();
        args->files->require_index = 1;
        for (i=0; i<args->nfnames; i++)
            if ( !bcf_sr_add_reader(args->files,args->fnames[i]) ) error("Failed to open, is the file indexed? %s\n", args->fnames[i]);
    }
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->nfnames; i++) free(args->fnames[i]);
    free(args->fnames);
    if ( args->files ) bcf_sr_destroy(args->files);
    if ( hts_close(args->out_fh)!=0 ) error("hts_close error\n");
    bcf_hdr_destroy(args->out_hdr);
    free(args->seen_seq);
}

int vcf_write_line(htsFile *fp, kstring_t *line);

static void concat(args_t *args)
{
    int i;
    if ( args->files )  // combining overlapping files, using synced reader 
    {
        while ( bcf_sr_next_line(args->files) )
        {
            for (i=0; i<args->files->nreaders; i++)
            {
                bcf1_t *line = bcf_sr_get_line(args->files,i);
                if ( !line ) continue;
                int id = bcf_hdr_name2id(args->out_hdr, bcf_seqname(args->files->readers[i].header, line));
                if ( id<0 ) error("FIXME: sequence name %s in %s\n", bcf_seqname(args->files->readers[i].header, line), args->fnames[i]);
                line->rid = id;
                bcf_write1(args->out_fh, args->out_hdr, line);
            }
        }
    }
    else    // concatenating
    {
        kstring_t tmp = {0,0,0};
        int prev_chr_id = -1, prev_pos;
        bcf1_t *line = bcf_init();
        for (i=0; i<args->nfnames; i++)
        {
            htsFile *fp = hts_open(args->fnames[i], "r"); if ( !fp ) error("Failed to open: %s\n", args->fnames[i]);
            bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) error("Failed to parse header: %s\n", args->fnames[i]);
            if ( !fp->is_bin && args->output_type&FT_VCF )
            {
                line->max_unpack = BCF_UN_STR;
                // if VCF is on both input and output, avoid VCF to BCF conversion
                while ( hts_getline(fp, KS_SEP_LINE, &fp->line) >=0 )
                {
                    char *str = fp->line.s;
                    while ( *str && *str!='\t' ) str++;
                    tmp.l = 0;
                    kputsn(fp->line.s,str-fp->line.s,&tmp);
                    int chr_id = bcf_hdr_name2id(args->out_hdr, tmp.s);
                    if ( chr_id<0 ) error("FIXME: sequence name %s in %s\n", tmp.s, args->fnames[i]);
                    if ( prev_chr_id!=chr_id ) 
                    {
                        prev_pos = -1;
                        if ( args->seen_seq[chr_id] )
                            error("\nThe chromosome block %s is not contiguous, consider running with -a.\n", tmp.s);
                    }
                    char *end;
                    int pos = strtol(str+1,&end,10) - 1;
                    if ( end==str+1 ) error("Could not parse line: %s\n", fp->line.s);
                    if ( prev_pos > pos )
                        error("The chromosome block %s is not sorted, consider running with -a.\n", tmp.s);
                    args->seen_seq[chr_id] = 1;
                    prev_chr_id = chr_id;

                    if ( vcf_write_line(args->out_fh, &fp->line)!=0 ) error("Failed to write %d bytes\n", fp->line.l);
                }
            }
            else
            {
                line->max_unpack = 0;
                // BCF conversion is required
                while ( bcf_read(fp, hdr, line)==0 )
                {
                    int chr_id = bcf_hdr_name2id(args->out_hdr, bcf_seqname(hdr, line));
                    if ( chr_id<0 ) error("FIXME: sequence name %s in %s\n", tmp.s, bcf_seqname(hdr, line));
                    if ( prev_chr_id!=chr_id ) 
                    {
                        prev_pos = -1;
                        if ( args->seen_seq[chr_id] )
                            error("\nThe chromosome block %s is not contiguous, consider running with -a.\n", bcf_seqname(hdr, line));
                    }
                    if ( prev_pos > line->pos )
                        error("The chromosome block %s is not sorted, consider running with -a.\n", bcf_seqname(hdr, line));
                    args->seen_seq[chr_id] = 1;
                    line->rid = chr_id;
                    prev_chr_id = chr_id;

                    if ( bcf_write(args->out_fh, args->out_hdr, line)!=0 ) error("Failed to write\n");
                }
            }
            bcf_hdr_destroy(hdr);
            hts_close(fp);
        }
        bcf_destroy(line);
        free(tmp.s);
    }
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Concatenate and/or combine overlapping VCF/BCF files split e.g. by chromosome\n");
    fprintf(stderr, "         or variant type. All source files must have the same sample columns appearing\n");
    fprintf(stderr, "         in the same order.\n");
    fprintf(stderr, "Usage:   bcftools concat [OPTIONS] <file> [...]\n");
    fprintf(stderr, "Options:\n");
	fprintf(stderr, "   -a, --allow-overlaps           Combine overlapping files. Requires indexed files.\n");
	fprintf(stderr, "   -l, --file-list <file>         Read the list of files from a file.\n");
	fprintf(stderr, "   -p, --phased-concat            Concatenate phased files.\n");
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
    args->output_type = FT_VCF;

    static struct option loptions[] = 
    {
        {"allow-overlap",1,0,'a'},
        {"phased-concat",1,0,'p'},
        {"output-type",1,0,'O'},
        {"file-list",1,0,'l'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "h:?O:l:ap",loptions,NULL)) >= 0) 
    {
        switch (c) {
    	    case 'a': args->allow_overlaps = 1; break;
    	    case 'p': args->phased_concat = 1; break;
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
    if ( args->phased_concat ) error("TODO: -p\n");
    while ( optind<argc ) 
    { 
        args->nfnames++;
        args->fnames = (char **)realloc(args->fnames,sizeof(char*)*args->nfnames);
        args->fnames[args->nfnames-1] = strdup(argv[optind]);
        optind++;
    }
    if ( args->file_list )
    {
        if ( args->nfnames ) error("Cannot combine -l with file names on command line.\n");
        args->fnames = hts_readlines(args->file_list, &args->nfnames);
    }
    if ( !args->nfnames ) usage(args);
    init_data(args);
    concat(args);
    destroy_data(args);
    free(args);
    return 0;
}
