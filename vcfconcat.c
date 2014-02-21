#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <math.h>
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

    // phasing
    int *start_pos, ifname;
    int *swap_phase, nswap, *nmatch, *nmism;
    bcf1_t **buf;
    int nbuf, mbuf, prev_chr, min_PQ, prev_pos_check;
    int32_t *GTa, *GTb, mGTa, mGTb, *phase_qual, *phase_set;

    char **argv, *file_list, **fnames;
    int argc, nfnames, allow_overlaps, phased_concat;
}
args_t;

static void init_data(args_t *args)
{
    bcf1_t *line = NULL;

    // With phased concat, the chunks overlap and come in the right order.  To
    // avoid opening all files at once, store start positions to recognise need
    // for the next one. This way we can keep only two open chunks at once.
    if ( args->phased_concat )
    {
        args->start_pos = (int*) malloc(sizeof(int)*args->nfnames);
        line = bcf_init();
    }

    kstring_t str = {0,0,0};
    int i, prev_chrid = -1;
    for (i=0; i<args->nfnames; i++)
    {
        htsFile *fp = hts_open(args->fnames[i], "r"); if ( !fp ) error("Failed to open: %s\n", args->fnames[i]);
        bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) error("Failed to parse header: %s\n", args->fnames[i]);
        if ( !args->out_hdr ) 
            args->out_hdr = bcf_hdr_dup(hdr);
        else
        {
            bcf_hdr_combine(args->out_hdr, hdr);

            if ( bcf_hdr_nsamples(hdr) != bcf_hdr_nsamples(args->out_hdr) ) 
                error("Different number of samples in %s. Perhaps \"bcftools merge\" is what you are looking for?\n", args->fnames[i]);

            int j;
            for (j=0; j<bcf_hdr_nsamples(hdr); j++)
                if ( strcmp(args->out_hdr->samples[j],hdr->samples[j]) )
                    error("Different sample names in %s. Perhaps \"bcftools merge\" is what you are looking for?\n", args->fnames[i]);
        }
        if ( args->phased_concat )
        {
            int ret = bcf_read(fp, hdr, line);
            if ( ret!=0 ) args->start_pos[i] = -2;
            else
            {
                int chrid = bcf_hdr_id2int(args->out_hdr,BCF_DT_CTG,bcf_seqname(hdr,line));
                args->start_pos[i] = chrid==prev_chrid ? line->pos : -1;
                prev_chrid = chrid;
            }
        }
        bcf_hdr_destroy(hdr);
        hts_close(fp);
    }
    free(str.s);
    if ( line ) bcf_destroy(line);

    args->seen_seq = (int*) calloc(args->out_hdr->n[BCF_DT_CTG],sizeof(int));

    // We do not output PS in order to save on VCF size. We can make a command line option if people need it.
    if ( args->phased_concat )
    {
        bcf_hdr_append(args->out_hdr,"##FORMAT=<ID=PQ,Number=1,Type=Integer,Description=\"Phasing Quality (bigger is better)\">");
        bcf_hdr_append(args->out_hdr,"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">");
    }
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
    else if ( args->phased_concat )
    {
        args->prev_chr = -1;
        args->swap_phase = (int*) calloc(bcf_hdr_nsamples(args->out_hdr),sizeof(int));
        args->nmatch = (int*) calloc(bcf_hdr_nsamples(args->out_hdr),sizeof(int));
        args->nmism  = (int*) calloc(bcf_hdr_nsamples(args->out_hdr),sizeof(int));
        args->phase_qual = (int32_t*) malloc(bcf_hdr_nsamples(args->out_hdr)*sizeof(int32_t));
        args->phase_set  = (int32_t*) malloc(bcf_hdr_nsamples(args->out_hdr)*sizeof(int32_t));
        args->files = bcf_sr_init();
        args->files->require_index = 1;
        args->ifname = 0;
        while ( args->ifname<args->nfnames && args->start_pos[args->ifname]==-2 ) args->ifname++;   // skip empty files
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
    free(args->start_pos);
    free(args->swap_phase);
    for (i=0; i<args->mbuf; i++) bcf_destroy(args->buf[i]);
    free(args->buf);
    free(args->GTa);
    free(args->GTb);
    free(args->nmatch);
    free(args->nmism);
    free(args->phase_qual);
    free(args->phase_set);
}

int vcf_write_line(htsFile *fp, kstring_t *line);

#define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }
static void phase_update(args_t *args, bcf_hdr_t *hdr, bcf1_t *rec)
{
    int i, nGTs = bcf_get_genotypes(hdr, rec, &args->GTa, &args->mGTa);
    for (i=0; i<bcf_hdr_nsamples(hdr); i++)
    {
        if ( !args->swap_phase[i] ) continue;
        int *gt = &args->GTa[i*2];
        if ( gt[0]==bcf_gt_missing || gt[1]==bcf_int32_vector_end ) continue;
        SWAP(int, gt[0], gt[1]);
        gt[1] |= 1;
    }
    bcf_update_genotypes(hdr,rec,args->GTa,nGTs);
}

static void phased_flush(args_t *args)
{
    if ( !args->nbuf ) return;

    bcf_hdr_t *ahdr = args->files->readers[0].header;
    bcf_hdr_t *bhdr = args->files->readers[1].header;

    int i, j, nsmpl = bcf_hdr_nsamples(args->out_hdr);

    for (i=0; i<args->nbuf; i+=2)
    {
        bcf1_t *arec = args->buf[i];
        bcf1_t *brec = args->buf[i+1];

        int nGTs = bcf_get_genotypes(ahdr, arec, &args->GTa, &args->mGTa);
        if ( nGTs < 0 ) error("GT is not present at %s:%d\n", bcf_seqname(ahdr,arec), arec->pos+1);
        if ( nGTs != 2*nsmpl ) error("TODO: Not diploid at %s:%d?\n", bcf_seqname(ahdr,arec), arec->pos+1);
        nGTs = bcf_get_genotypes(bhdr, brec, &args->GTb, &args->mGTb);
        if ( nGTs < 0 ) error("GT is not present at %s:%d\n", bcf_seqname(bhdr,brec), brec->pos+1);
        if ( nGTs != 2*nsmpl ) error("TODO: Not diploid at %s:%d?\n", bcf_seqname(bhdr,brec), brec->pos+1);

        for (j=0; j<nsmpl; j++)
        {
            int *gta = &args->GTa[j*2];
            int *gtb = &args->GTb[j*2];
            if ( gta[1]==bcf_int32_vector_end || gtb[1]==bcf_int32_vector_end ) continue;
            if ( gta[0]==bcf_gt_missing || gta[1]==bcf_gt_missing || gtb[0]==bcf_gt_missing || gtb[1]==bcf_gt_missing ) continue;
            if ( !bcf_gt_is_phased(gta[1]) || !bcf_gt_is_phased(gtb[1]) ) continue;
            if ( bcf_gt_allele(gta[0])==bcf_gt_allele(gta[1]) || bcf_gt_allele(gtb[0])==bcf_gt_allele(gtb[1]) ) continue;
            if ( bcf_gt_allele(gta[0])==bcf_gt_allele(gtb[0]) && bcf_gt_allele(gta[1])==bcf_gt_allele(gtb[1]) )
            {
                if ( args->swap_phase[j] ) args->nmism[j]++; else args->nmatch[j]++;
            }
            if ( bcf_gt_allele(gta[0])==bcf_gt_allele(gtb[1]) && bcf_gt_allele(gta[1])==bcf_gt_allele(gtb[0]) )
            {
                if ( args->swap_phase[j] ) args->nmatch[j]++; else args->nmism[j]++;
            }
        }
    }
    for (i=0; i<args->nbuf/2; i+=2)
    {
        bcf1_t *arec = args->buf[i];
        bcf_translate(args->out_hdr, args->files->readers[0].header, arec);
        if ( args->nswap )
            phase_update(args, args->out_hdr, arec);
        bcf_update_format_int32(args->out_hdr,arec,"PS",args->phase_set,nsmpl);
        bcf_write(args->out_fh, args->out_hdr, arec);

        if ( arec->pos < args->prev_pos_check ) error("FIXME, disorder: %s:%d vs %d  [1]\n", bcf_seqname(args->files->readers[0].header,arec),arec->pos+1,args->prev_pos_check+1);
        args->prev_pos_check = arec->pos;
    }
    args->nswap = 0;
    for (j=0; j<nsmpl; j++)
    {
        if ( args->nmatch[j] >= args->nmism[j] )
            args->swap_phase[j] = 0;
        else
        {
            args->swap_phase[j] = 1;
            args->nswap++;
        }
        if ( args->nmatch[j] && args->nmism[j] )
        {
            // Entropy-inspired quality. The factor 0.7 shifts and scales to (0,1)
            double f = (double)args->nmatch[j]/(args->nmatch[j]+args->nmism[j]);
            args->phase_qual[j] = 99*(0.7 + f*log(f) + (1-f)*log(1-f))/0.7;
        }
        else
            args->phase_qual[j] = 99;
        args->nmatch[j] = 0;
        args->nmism[j]  = 0;
    }
    int PQ_printed = 0;
    for (; i<args->nbuf; i+=2)
    {
        bcf1_t *brec = args->buf[i+1];
        bcf_translate(args->out_hdr, args->files->readers[1].header, brec);
        if ( !PQ_printed )
        {
            bcf_update_format_int32(args->out_hdr,brec,"PQ",args->phase_qual,nsmpl);
            PQ_printed = 1;
            for (j=0; j<nsmpl; j++)
                if ( args->phase_qual[j] < args->min_PQ ) args->phase_set[j] = brec->pos+1;
        }
        if ( args->nswap )
            phase_update(args, args->out_hdr, brec);
        bcf_update_format_int32(args->out_hdr,brec,"PS",args->phase_set,nsmpl);
        bcf_write(args->out_fh, args->out_hdr, brec);

        if ( brec->pos < args->prev_pos_check ) error("FIXME, disorder: %s:%d vs %d  [2]\n", bcf_seqname(args->files->readers[1].header,brec),brec->pos+1,args->prev_pos_check+1);
        args->prev_pos_check = brec->pos;
    }
    args->nbuf = 0;
}

static void phased_push(args_t *args, bcf1_t *arec, bcf1_t *brec)
{
    int i, nsmpl = bcf_hdr_nsamples(args->out_hdr);
    int chr_id = bcf_hdr_name2id(args->out_hdr, bcf_seqname(args->files->readers[0].header,arec));
    if ( args->prev_chr<0 || args->prev_chr!=chr_id )
    {
        if ( args->prev_chr>=0 ) phased_flush(args);

        for (i=0; i<nsmpl; i++) 
            args->phase_set[i] = arec->pos+1;

        args->prev_chr = chr_id;
        args->prev_pos_check = -1;
    }

    if ( !brec )
    {
        bcf_translate(args->out_hdr, args->files->readers[0].header, arec);
        if ( args->nswap )
            phase_update(args, args->out_hdr, arec);
        bcf_update_format_int32(args->out_hdr,arec,"PS",args->phase_set,nsmpl);
        bcf_write(args->out_fh, args->out_hdr, arec);

        if ( arec->pos < args->prev_pos_check ) error("FIXME, disorder: %s:%d vs %d  [3]\n", bcf_seqname(args->files->readers[0].header,arec), arec->pos+1,args->prev_pos_check+1);
        args->prev_pos_check = arec->pos;
        return;
    }

    int m = args->mbuf;
    args->nbuf += 2;
    hts_expand(bcf1_t*,args->nbuf,args->mbuf,args->buf);
    for (i=m; i<args->mbuf; i++)
        args->buf[i] = bcf_init1();

    SWAP(bcf1_t*, args->files->readers[0].buffer[0], args->buf[args->nbuf-2]);
    SWAP(bcf1_t*, args->files->readers[1].buffer[0], args->buf[args->nbuf-1]);
}

static void concat(args_t *args)
{
    int i;
    if ( args->phased_concat )  // phased concat
    {
        // keep only two open files at a time
        while ( args->ifname < args->nfnames )
        {
            if ( !bcf_sr_add_reader(args->files,args->fnames[args->ifname]) ) error("Failed to open %s\n", args->fnames[args->ifname]);
            args->ifname++;
            while ( args->ifname<args->nfnames && args->start_pos[args->ifname]==-2 ) args->ifname++;   // skip empty files
            int seek_pos = -1;
            int seek_chr = -1;
            if ( args->files->nreaders < 2 && args->ifname < args->nfnames )
            {
                if ( !bcf_sr_add_reader(args->files,args->fnames[args->ifname]) ) error("Failed to open %s\n", args->fnames[args->ifname]);
                args->ifname++;
                while ( args->ifname<args->nfnames && args->start_pos[args->ifname]==-2 ) args->ifname++;   // skip empty files
                if ( bcf_sr_has_line(args->files,0) )
                {
                    bcf1_t *line = bcf_sr_get_line(args->files,0);
                    bcf_sr_seek(args->files, bcf_seqname(args->files->readers[0].header,line), line->pos);
                    seek_pos = line->pos;
                    seek_chr = bcf_hdr_name2id(args->out_hdr, bcf_seqname(args->files->readers[0].header,line));
                }
            }

            int nret;
            while ( (nret = bcf_sr_next_line(args->files)) )
            {
                if ( !bcf_sr_has_line(args->files,0) )  // no input from the first reader
                {
                    // We are assuming that there is a perfect overlap, sites which are not present in both files are dropped
                    if ( ! bcf_sr_region_done(args->files,0) ) continue;

                    phased_flush(args);
                    bcf_sr_remove_reader(args->files, 0);
                }

                // Set chromosome, the IDs may be different across the files
                for (i=0; i<args->files->nreaders; i++)
                    if ( bcf_sr_has_line(args->files,i) ) break;
                bcf1_t *line = bcf_sr_get_line(args->files,i);

                // This can happen after bcf_sr_seek: indel may start before the coordinate which we seek to.
                if ( seek_chr>=0 && seek_pos>line->pos && seek_chr==bcf_hdr_name2id(args->out_hdr, bcf_seqname(args->files->readers[i].header,line)) ) continue;
                seek_pos = seek_chr = -1;

                int must_seek = 0;
                while ( args->ifname < args->nfnames && line->pos >= args->start_pos[args->ifname] )   // overlaps with the next reader
                {
                    must_seek = 1;
                    bcf_sr_add_reader(args->files,args->fnames[args->ifname]);
                    args->ifname++;
                    while ( args->ifname<args->nfnames && args->start_pos[args->ifname]==-2 ) args->ifname++;   // skip empty files
                }
                if ( must_seek )
                {
                    bcf_sr_seek(args->files, bcf_seqname(args->files->readers[i].header,line), line->pos);
                    seek_pos = line->pos;
                    seek_chr = bcf_hdr_name2id(args->out_hdr, bcf_seqname(args->files->readers[i].header,line));
                    continue;
                }

                // We are assuming that there is a perfect overlap, sites which are not present in both files are dropped
                if ( args->files->nreaders>1 && !bcf_sr_has_line(args->files,1) && !bcf_sr_region_done(args->files,1) ) continue;

                phased_push(args, bcf_sr_get_line(args->files,0), args->files->nreaders>1 ? bcf_sr_get_line(args->files,1) : NULL);
            }
            args->ifname++;
        }
    }
    else if ( args->files )  // combining overlapping files, using synced reader 
    {
        while ( bcf_sr_next_line(args->files) )
        {
            for (i=0; i<args->files->nreaders; i++)
            {
                bcf1_t *line = bcf_sr_get_line(args->files,i);
                if ( !line ) continue;
                bcf_translate(args->out_hdr, args->files->readers[i].header, line);
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
                // BCF conversion is required
                line->max_unpack = 0;
                while ( bcf_read(fp, hdr, line)==0 )
                {
                    bcf_translate(args->out_hdr, hdr, line);

                    if ( prev_chr_id!=line->rid ) 
                    {
                        prev_pos = -1;
                        if ( args->seen_seq[line->rid] )
                            error("\nThe chromosome block %s is not contiguous, consider running with -a.\n", bcf_seqname(args->out_hdr, line));
                    }
                    if ( prev_pos > line->pos )
                        error("The chromosome block %s is not sorted, consider running with -a.\n", bcf_seqname(args->out_hdr, line));
                    args->seen_seq[line->rid] = 1;
                    prev_chr_id = line->rid;

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
    fprintf(stderr, "Usage:   bcftools concat [options] <A.bcf>|<A.vcf>|<A.vcf.gz> [<B.bcf>|<B.vcf>|<B.vcf.gz> ...]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
	fprintf(stderr, "   -a, --allow-overlaps           Combine overlapping files. Requires indexed files.\n");
	fprintf(stderr, "   -l, --file-list <file>         Read the list of files from a file.\n");
	fprintf(stderr, "   -p, --phased-concat            Concatenate phased files.\n");
	fprintf(stderr, "   -q, --min-PQ <int>             Break phase set if phasing quality is lower than <int> [30]\n");
	fprintf(stderr, "   -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfconcat(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->output_type = FT_VCF;
    args->min_PQ  = 30;

    static struct option loptions[] = 
    {
        {"allow-overlap",1,0,'a'},
        {"phased-concat",1,0,'p'},
        {"output-type",1,0,'O'},
        {"file-list",1,0,'l'},
        {"min-PQ",1,0,'q'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "h:?O:l:apq:",loptions,NULL)) >= 0) 
    {
        switch (c) {
    	    case 'q': args->min_PQ = atoi(optarg); break;
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
