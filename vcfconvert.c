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
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <htslib/kseq.h>
#include "bcftools.h"
#include "filter.h"
#include "convert.h"
#include "tsv2vcf.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct _args_t args_t;
struct _args_t
{
    faidx_t *ref;
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE
    convert_t *convert;
    bcf_srs_t *files;
    bcf_hdr_t *header;
    void (*convert_func)(struct _args_t *);
    struct {
        int total, skipped, hom_rr, het_ra, hom_aa, het_aa, missing; 
    } n;
    kstring_t str;
    int32_t *gts;
    int nsamples, *samples, sample_is_file, targets_is_file, regions_is_file, output_type;
    char **argv, *sample_list, *targets_list, *regions_list, *tag, *columns;
    char *outfname, *infname, *ref_fname;
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

static void hapsample_to_vcf(args_t *args)
{
    kstring_t str = {0,0,0}, line = {0,0,0};

    char *hap_fname = NULL, *sample_fname = NULL;
    sample_fname = strchr(args->infname,',');
    if ( !sample_fname )
    {
        ksprintf(&str,"%s.hap.gz", args->infname);
        hap_fname = strdup(str.s);
        str.l = 0;
        ksprintf(&str,"%s.samples", args->infname);
        sample_fname = strdup(str.s);
    }
    else
    {
        *sample_fname = 0;
        hap_fname = strdup(args->infname);
        sample_fname = strdup(sample_fname+1);
    }
    htsFile *hap_fh = hts_open(hap_fname, "r");
    if ( !hap_fh ) error("Could not read: %s\n", hap_fname);
    if ( hts_getline(hap_fh, KS_SEP_LINE, &line) <= 0 ) error("Empty file: %s\n", hap_fname);

    // Find out the chromosome name, sample names, init and print the VCF header
    str.l = 0;
    char *ss, *se = strchr(line.s,':');
    if ( !se ) error("Expected CHROM:POS_REF_ALT in first column of %s\n", hap_fname);
    kputsn(line.s, se-line.s, &str);

    args->header = bcf_hdr_init("w");
    bcf_hdr_append(args->header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_printf(args->header, "##contig=<ID=%s,length=%d>", str.s,0x7fffffff);   // MAX_CSI_COOR

    int i, nsamples;
    char **samples = hts_readlist(sample_fname, 1, &nsamples);
    for (i=2; i<nsamples; i++)
    {
        se = samples[i]; while ( *se && !isspace(*se) ) se++;
        *se = 0;
        bcf_hdr_add_sample(args->header,samples[i]);
    }
    bcf_hdr_add_sample(args->header,NULL);
    for (i=0; i<nsamples; i++) free(samples[i]);
    free(samples);

    htsFile *out_fh = hts_open(args->outfname,hts_bcf_wmode(args->output_type));
    bcf_hdr_write(out_fh,args->header);
    bcf1_t *rec = bcf_init();

    nsamples -= 2;
    args->gts = (int32_t *) malloc(sizeof(int32_t)*nsamples*2);

    do
    {
        bcf_clear(rec);

        // CHROM:POS_REF_ALT
        ss = se = line.s;
        while ( *se && *se!=':' ) se++;
        if ( !*se ) error("Could not parse: %s\n", line.s);
        char tmp = *se; *se = 0;
        rec->rid = bcf_hdr_name2id(args->header,ss); *se = tmp;

        // POS
        rec->pos = strtol(se+1,&ss,10);
        if ( ss==se+1 ) error("Could not parse the POS part in CHROM:POS_REF_ALT: %s\n", line.s);
        rec->pos--;

        // REF,ALT
        str.l = 0;
        se = ++ss;
        while ( *se && *se!='_' ) se++; 
        if ( !*se ) error("Could not parse: %s\n", line.s);
        kputsn(ss,se-ss,&str);
        ss = ++se;
        while ( *se && !isspace(*se) ) se++;
        if ( !*se ) error("Could not parse: %s\n", line.s);
        kputc(',',&str);
        kputsn(ss,se-ss,&str);
        bcf_update_alleles_str(args->header, rec, str.s);

        // ID,POS: skip
        ss = ++se;
        while ( *se && !isspace(*se) ) se++;
        if ( !*se ) error("Could not parse: %s\n", line.s);
        ss = ++se;
        while ( *se && !isspace(*se) ) se++;
        if ( !*se ) error("Could not parse: %s\n", line.s);

        // REF,ALT
        int rev_als = 0;
        ss = ++se;
        while ( *se && !isspace(*se) ) se++;
        if ( !*se ) error("Could not parse: %s\n", line.s);
        tmp = *se; *se = 0;
        if ( strcmp(ss,rec->d.allele[0]) )
        {
            if ( strcmp(ss,rec->d.allele[1]) ) { *se = tmp; error("REF/ALT mismatch: %s\n", line.s); }
            rev_als = 1;
        }
        *se = tmp;
        ss = ++se;
        while ( *se && !isspace(*se) ) se++;
        if ( !*se ) error("Could not parse: %s\n", line.s);
        tmp = *se; *se = 0;
        if ( !rev_als && strcmp(ss,rec->d.allele[1]) ) { *se = tmp; error("REF/ALT mismatch: %s\n", line.s); }
        else if ( rev_als && strcmp(ss,rec->d.allele[0]) ) { *se = tmp; error("REF/ALT mismatch: %s\n", line.s); }
        *se = tmp;
        ss = se + 1;

        // Haplotypes
        int32_t a0, a1;
        if ( rev_als ) { a0 = bcf_gt_phased(1); a1 = bcf_gt_phased(0); }
        else { a0 = bcf_gt_phased(0); a1 = bcf_gt_phased(1); }
        for (i=0; i<nsamples; i++)
        {
            args->gts[2*i]   = ss[4*i]  =='0' ? a0 : a1;
            args->gts[2*i+1] = ss[4*i+2]=='0' ? a0 : a1;
        }

        if ( bcf_update_genotypes(args->header,rec,args->gts,nsamples*2) ) error("Could not update GT field\n");
        bcf_write(out_fh, args->header,rec);
    }
    while ( hts_getline(hap_fh, KS_SEP_LINE, &line)>0 );


    if ( hts_close(out_fh) ) error("Close failed: %s\n", args->outfname);
    if ( hts_close(hap_fh) ) error("Close failed: %s\n", hap_fname);
    bcf_hdr_destroy(args->header);
    bcf_destroy(rec);
    free(sample_fname);
    free(hap_fname);
    free(str.s);
    free(line.s);
    free(args->gts);
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

    int ret, gen_compressed = 1, sample_compressed = 0;
    char *gen_fname = NULL, *sample_fname = NULL;
    str.l = 0;
    kputs(args->outfname,&str);
    int n_files, i;
    char **files = hts_readlist(str.s, 0, &n_files);
    if ( n_files==1 )
    {
        int l = str.l;
        kputs(".samples",&str);
        sample_fname = strdup(str.s);
        str.l = l;
        kputs(".gen.gz",&str);
        gen_fname = strdup(str.s);
    }
    else if ( n_files==2 )
    {
        if (strlen(files[0]) && strcmp(files[0],".")!=0) gen_fname = strdup(files[0]);
        if (strlen(files[1]) && strcmp(files[1],".")!=0) sample_fname = strdup(files[1]);
    }
    else
    {
        error("Error parsing --gensample filenames: %s\n", args->outfname);
    }
    for (i=0; i<n_files; i++) free(files[i]);
    free(files);

    if ( gen_fname && (strlen(gen_fname)<3 || strcasecmp(".gz",gen_fname+strlen(gen_fname)-3)) ) gen_compressed = 0;
    if ( sample_fname && strlen(sample_fname)>3 && strcasecmp(".gz",sample_fname+strlen(sample_fname)-3)==0 ) sample_compressed = 0;

    if (gen_fname) fprintf(stderr, "Gen file: %s\n", gen_fname);
    if (sample_fname) fprintf(stderr, "Sample file: %s\n", sample_fname);

    // write samples file
    if (sample_fname) {
        int i;
        BGZF *sout = bgzf_open(sample_fname, sample_compressed ? "w" : "wu");
        if ( sample_compressed ) sout->is_gzip = 1;     // impute2 can't read bgzf
        str.l = 0;
        kputs("ID_1 ID_2 missing\n0 0 0\n", &str);
        ret = bgzf_write(sout, str.s, str.l);
        if ( ret != str.l ) error("Error writing %s: %s\n", sample_fname, strerror(errno));
        for (i=0; i<bcf_hdr_nsamples(args->header); i++)
        {
            str.l = 0;
            ksprintf(&str, "%s %s 0\n", args->header->samples[i],args->header->samples[i]);
            ret = bgzf_write(sout, str.s, str.l);
            if ( ret != str.l ) error("Error writing %s: %s\n", sample_fname, strerror(errno));
        }
        if ( bgzf_close(sout)!=0 ) error("Error closing %s: %s\n", sample_fname, strerror(errno));
        free(sample_fname);
    }
    if (!gen_fname) {
        if ( str.m ) free(str.s);
        return;
    }

    int no_alt = 0, non_biallelic = 0, filtered = 0;
    BGZF *gout = bgzf_open(gen_fname, gen_compressed ? "w" : "wu");
    if ( gen_compressed ) gout->is_gzip = 1;
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) { filtered++; continue; }
        }

        // ALT allele is required
        if ( line->n_allele<2 ) { no_alt++; continue; }
        // biallelic required
        if ( line->n_allele>2 ) {
            if (!non_biallelic)
                fprintf(stderr, "Warning: non-biallelic records are skipped. Consider splitting multi-allelic records into biallelic records using 'bcftools norm -m-'.\n");
            non_biallelic++;
            continue;
        }

        str.l = 0;
        convert_line(args->convert, line, &str);
        if ( str.l )
        {
            int ret = bgzf_write(gout, str.s, str.l);
            if ( ret!= str.l ) error("Error writing %s: %s\n", gen_fname,strerror(errno));
        }
    }
    fprintf(stderr, "%d records skipped: %d/%d/%d no-ALT/non-biallelic/filtered\n", no_alt+non_biallelic+filtered, no_alt, non_biallelic, filtered);
    if ( str.m ) free(str.s);
    if ( bgzf_close(gout)!=0 ) error("Error closing %s: %s\n", gen_fname,strerror(errno));
    free(gen_fname);
}

static void vcf_to_haplegendsample(args_t *args)
{
    kstring_t str = {0,0,0};
    kputs("%_GT_TO_HAP\n", &str);
    open_vcf(args,str.s);

    int ret, hap_compressed = 1, legend_compressed = 1, sample_compressed = 0;
    char *hap_fname = NULL, *legend_fname = NULL, *sample_fname = NULL;
    str.l = 0;
    kputs(args->outfname,&str);
    int n_files, i;
    char **files = hts_readlist(str.s, 0, &n_files);
    if ( n_files==1 )
    {
        int l = str.l;
        kputs(".samples",&str);
        sample_fname = strdup(str.s);
        str.l = l;
        kputs(".legend.gz",&str);
        legend_fname = strdup(str.s);
        str.l = l;
        kputs(".hap.gz",&str);
        hap_fname = strdup(str.s);
    }
    else if ( n_files==3 )
    {
        if (strlen(files[0]) && strcmp(files[0],".")!=0) hap_fname = strdup(files[0]);
        if (strlen(files[1]) && strcmp(files[1],".")!=0) legend_fname = strdup(files[1]);
        if (strlen(files[2]) && strcmp(files[2],".")!=0) sample_fname = strdup(files[2]);
    }
    else
    {
        error("Error parsing --hapslegendsample filenames: %s\n", args->outfname);
    }
    for (i=0; i<n_files; i++) free(files[i]);
    free(files);

    if ( hap_fname && (strlen(hap_fname)<3 || strcasecmp(".gz",hap_fname+strlen(hap_fname)-3)) ) hap_compressed = 0;
    if ( legend_fname && (strlen(legend_fname)<3 || strcasecmp(".gz",legend_fname+strlen(legend_fname)-3)) ) legend_compressed = 0;
    if ( sample_fname && strlen(sample_fname)>3 && strcasecmp(".gz",sample_fname+strlen(sample_fname)-3)==0 ) sample_compressed = 0;

    if (hap_fname) fprintf(stderr, "Haps file: %s\n", hap_fname);
    if (legend_fname) fprintf(stderr, "Legend file: %s\n", legend_fname);
    if (sample_fname) fprintf(stderr, "Sample file: %s\n", sample_fname);

    // write samples file
    if (sample_fname) {
        int i;
        BGZF *sout = bgzf_open(sample_fname, sample_compressed ? "w" : "wu");
        if ( sample_compressed ) sout->is_gzip = 1;
        str.l = 0;
        kputs("sample population group sex\n", &str);
        ret = bgzf_write(sout, str.s, str.l);
        if ( ret != str.l ) error("Error writing %s: %s\n", sample_fname, strerror(errno));
        for (i=0; i<bcf_hdr_nsamples(args->header); i++)
        {
            str.l = 0;
            ksprintf(&str, "%s %s %s 2\n", args->header->samples[i], args->header->samples[i], args->header->samples[i]);
            ret = bgzf_write(sout, str.s, str.l);
            if ( ret != str.l ) error("Error writing %s: %s\n", sample_fname, strerror(errno));
        }
        if ( bgzf_close(sout)!=0 ) error("Error closing %s: %s\n", sample_fname, strerror(errno));
        free(sample_fname);
    }
    if (!hap_fname && !legend_fname) {
        if ( str.m ) free(str.s);
        return;
    }

    // open haps and legend outputs
    BGZF *hout = hap_fname ? bgzf_open(hap_fname, hap_compressed ? "w" : "wu") : NULL;
    BGZF *lout = legend_fname ? bgzf_open(legend_fname, legend_compressed ? "w" : "wu") : NULL;
    if ( hout && hap_compressed ) hout->is_gzip = 1;    // impute2 can't read BGZF
    if ( lout && legend_compressed ) lout->is_gzip = 1;
    if (legend_fname) {
        str.l = 0;
        kputs("id position a0 a1\n", &str);
        ret = bgzf_write(lout, str.s, str.l);
        if ( ret != str.l ) error("Error writing %s: %s\n", legend_fname, strerror(errno));
    }

    int no_alt = 0, non_biallelic = 0, filtered = 0;
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) { filtered++; continue; }
        }

        // ALT allele is required
        if ( line->n_allele<2 ) { no_alt++; continue; }
        // biallelic required
        if ( line->n_allele>2 ) {
            if (!non_biallelic)
                fprintf(stderr, "Warning: non-biallelic records are skipped. Consider splitting multi-allelic records into biallelic records using 'bcftools norm -m-'.\n");
            non_biallelic++;
            continue;
        }

        str.l = 0;
        convert_line(args->convert, line, &str);
        if ( str.l )
        {
            // write haps file
            if (hap_fname) {
                ret = bgzf_write(hout, str.s, str.l); // write hap file
                if ( ret != str.l ) error("Error writing %s: %s\n", hap_fname, strerror(errno));
            }
            if (legend_fname) {
                str.l = 0;
                if (line->d.id[0]!='.' || line->d.id[1]!=0) {
                    ksprintf(&str, "%s %d %s %s\n", line->d.id, line->pos+1, line->d.allele[0], line->d.allele[1]);
                }
                else {
                    int vartype = bcf_get_variant_types(line);
                    if (vartype&VCF_SNP) {
                        ksprintf(&str, "%s:%d_%s_%s %d %s %s\n", bcf_seqname(args->header, line), line->pos+1, line->d.allele[0], line->d.allele[1], line->pos+1, line->d.allele[0], line->d.allele[1]);
                    }
                    else {
                        char type;
                        if (vartype&VCF_INDEL) {
                            if (strlen(line->d.allele[0]) > strlen(line->d.allele[1]))
                                type = 'D';
                            else
                                type = 'I';
                        }
                        else if (vartype&VCF_MNP)
                            type = 'M';
                        else
                            type = 'O';
                        ksprintf(&str, "%s:%d_%s_%s:%c %d %s %s\n", bcf_seqname(args->header, line), line->pos+1, line->d.allele[0], line->d.allele[1], type, line->pos+1, line->d.allele[0], line->d.allele[1]);
                    }
                }
                // write legend file
                ret = bgzf_write(lout, str.s, str.l);
                if ( ret != str.l ) error("Error writing %s: %s\n", legend_fname, strerror(errno));
            }
        }
    }
    fprintf(stderr, "%d records skipped: %d/%d/%d no-ALT/non-biallelic/filtered\n", no_alt+non_biallelic+filtered, no_alt, non_biallelic, filtered);
    if ( str.m ) free(str.s);
    if ( hout && bgzf_close(hout)!=0 ) error("Error closing %s: %s\n", hap_fname, strerror(errno));
    if ( lout && bgzf_close(lout)!=0 ) error("Error closing %s: %s\n", legend_fname, strerror(errno));
    if (hap_fname) free(hap_fname);
    if (legend_fname) free(legend_fname);
}

static void bcf_hdr_set_chrs(bcf_hdr_t *hdr, faidx_t *fai)
{
    int i, n = faidx_nseq(fai);
    for (i=0; i<n; i++)
    {
        const char *seq = faidx_iseq(fai,i);
        int len = faidx_seq_len(fai, seq);
        bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%d>", seq,len);
    }
}
static inline int acgt_to_5(char base)
{
    if ( base=='A' ) return 0;
    if ( base=='C' ) return 1;
    if ( base=='G' ) return 2;
    if ( base=='T' ) return 3;
    return 4;
}
static inline int tsv_setter_aa1(args_t *args, char *ss, char *se, int alleles[], int *nals, int ref, int32_t *gts)
{
    if ( se - ss > 2 ) return -1;   // currently only SNPs

    if ( ss[0]=='-' )
    {
        // missing GT
        gts[0] = bcf_gt_missing;
        gts[1] = bcf_int32_vector_end;
        args->n.missing++;
        return 0;
    }
    if ( ss[0]=='I' ) return -2;    // skip insertions/deletions for now
    if ( ss[0]=='D' ) return -2;

    int a0 = acgt_to_5(toupper(ss[0]));
    int a1 = ss[1] ? acgt_to_5(toupper(ss[1])) : a0;
    if ( alleles[a0]<0 ) alleles[a0] = (*nals)++;
    if ( alleles[a1]<0 ) alleles[a1] = (*nals)++;

    gts[0] = bcf_gt_unphased(alleles[a0]); 
    gts[1] = ss[1] ? bcf_gt_unphased(alleles[a1]) : bcf_int32_vector_end;

    if ( ref==a0 && ref==a1  ) args->n.hom_rr++;    // hom ref: RR
    else if ( ref==a0 ) args->n.het_ra++;           // het: RA
    else if ( ref==a1 ) args->n.het_ra++;           // het: AR
    else if ( a0==a1 ) args->n.hom_aa++;            // hom-alt: AA
    else args->n.het_aa++;                          // non-ref het: AA

    return 0;
}
static int tsv_setter_aa(tsv_t *tsv, bcf1_t *rec, void *usr)
{
    args_t *args = (args_t*) usr;

    int len;
    char *ref = faidx_fetch_seq(args->ref, (char*)bcf_hdr_id2name(args->header,rec->rid), rec->pos, rec->pos, &len);
    if ( !ref ) error("faidx_fetch_seq failed at %s:%d\n", bcf_hdr_id2name(args->header,rec->rid), rec->pos+1);

    int nals = 1, alleles[5] = { -1, -1, -1, -1, -1 };    // a,c,g,t,n
    ref[0] = toupper(ref[0]);
    int iref = acgt_to_5(ref[0]);
    alleles[iref] = 0;

    rec->n_sample = bcf_hdr_nsamples(args->header);

    int i, ret;
    for (i=0; i<rec->n_sample; i++)
    {
        if ( i>0 )
        {
            ret = tsv_next(tsv);
            if ( ret==-1 ) error("Too few columns for %d samples at %s:%d\n", rec->n_sample,bcf_hdr_id2name(args->header,rec->rid), rec->pos+1);
        }
        ret = tsv_setter_aa1(args, tsv->ss, tsv->se, alleles, &nals, iref, args->gts+i*2);
        if ( ret==-1 ) error("Error parsing the site %s:%d, expected two characters\n", bcf_hdr_id2name(args->header,rec->rid), rec->pos+1);
        if ( ret==-2 ) 
        {
            // something else than a SNP
            free(ref);
            return -1;
        }
    }

    args->str.l = 0;
    kputc(ref[0], &args->str);
    for (i=0; i<5; i++) 
    {
        if ( alleles[i]>0 )
        {
            kputc(',', &args->str);
            kputc("ACGTN"[i], &args->str);
        }
    }
    bcf_update_alleles_str(args->header, rec, args->str.s);
    if ( bcf_update_genotypes(args->header,rec,args->gts,rec->n_sample*2) ) error("Could not update the GT field\n");

    free(ref);
    return 0;
}

static void tsv_to_vcf(args_t *args)
{
    if ( !args->ref_fname ) error("Missing the --ref option\n");
    if ( !args->sample_list ) error("Missing the --samples option\n");

    args->ref = fai_load(args->ref_fname);
    if ( !args->ref ) error("Could not load the reference %s\n", args->ref_fname);

    args->header = bcf_hdr_init("w");
    bcf_hdr_set_chrs(args->header, args->ref);
    bcf_hdr_append(args->header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

    int i, n;
    char **smpls = hts_readlist(args->sample_list, args->sample_is_file, &n);
    if ( !smpls ) error("Could not parse %s\n", args->sample_list);
    for (i=0; i<n; i++)
    {
        bcf_hdr_add_sample(args->header, smpls[i]);
        free(smpls[i]);
    }
    free(smpls);
    bcf_hdr_add_sample(args->header, NULL);
    args->gts = (int32_t *) malloc(sizeof(int32_t)*n*2);

    htsFile *out_fh = hts_open(args->outfname,hts_bcf_wmode(args->output_type));
    bcf_hdr_write(out_fh,args->header);

    tsv_t *tsv = tsv_init(args->columns ? args->columns : "ID,CHROM,POS,AA");
    if ( tsv_register(tsv, "CHROM", tsv_setter_chrom, args->header) < 0 ) error("Expected CHROM column\n");
    if ( tsv_register(tsv, "POS", tsv_setter_pos, NULL) < 0 ) error("Expected POS column\n");
    if ( tsv_register(tsv, "ID", tsv_setter_id, args->header) < 0 ) error("Expected ID column\n");
    if ( tsv_register(tsv, "AA", tsv_setter_aa, args) < 0 ) error("Expected AA column\n");

    bcf1_t *rec = bcf_init();
    bcf_float_set_missing(rec->qual);

    kstring_t line = {0,0,0};
    htsFile *in_fh = hts_open(args->infname, "r");
    if ( !in_fh ) error("Could not read: %s\n", args->infname);
    while ( hts_getline(in_fh, KS_SEP_LINE, &line) > 0 )
    {
        if ( line.s[0]=='#' ) continue;     // skip comments
        bcf_clear(rec);

        args->n.total++;
        if ( !tsv_parse(tsv, rec, line.s) )
            bcf_write(out_fh, args->header, rec);
        else
            args->n.skipped++;
    }
    if ( hts_close(in_fh) ) error("Close failed: %s\n", args->infname);
    free(line.s);

    fai_destroy(args->ref);
    bcf_hdr_destroy(args->header);
    hts_close(out_fh);
    tsv_destroy(tsv);
    bcf_destroy(rec);
    free(args->str.s);
    free(args->gts);

    fprintf(stderr,"Rows total: \t%d\n", args->n.total);
    fprintf(stderr,"Rows skipped: \t%d\n", args->n.skipped);
    fprintf(stderr,"Missing GTs: \t%d\n", args->n.missing);
    fprintf(stderr,"Hom RR: \t%d\n", args->n.hom_rr);
    fprintf(stderr,"Het RA: \t%d\n", args->n.het_ra);
    fprintf(stderr,"Hom AA: \t%d\n", args->n.hom_aa);
    fprintf(stderr,"Het AA: \t%d\n", args->n.het_aa);
}

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Converts VCF/BCF to other formats and back. See man page for file\n");
    fprintf(stderr, "         formats details. When specifying output files explicitly instead\n");
    fprintf(stderr, "         of with <prefix>, one can use '-' for stdout and '.' to suppress.\n");
    fprintf(stderr, "Usage:   bcftools convert [OPTIONS] <input_file>\n");
    fprintf(stderr, "\n");
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
    fprintf(stderr, "VCF output options:\n");
    fprintf(stderr, "   -o, --output <file>            output file name [stdout]\n");
    fprintf(stderr, "   -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "GEN/SAMPLE options:\n");
    fprintf(stderr, "   -g, --gensample             <prefix>|<gen-file>,<sample-file>\n");
    fprintf(stderr, "       --tag <string>          tag to take values for .gen file: GT,PL,GL,GP [GT]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "HAP/SAMPLE options (output from SHAPEIT):\n");
    fprintf(stderr, "       --hapsample2vcf         <prefix>|<haps-file>,<sample-file>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "HAP/LEGEND/SAMPLE options:\n");
    fprintf(stderr, "   -h, --haplegendsample       <prefix>|<hap-file>,<legend-file>,<sample-file>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "TSV options:\n");
    fprintf(stderr, "       --tsv2vcf <file>        \n");
    fprintf(stderr, "   -c, --columns <string>      columns of the input tsv file [CHROM,POS,ID,AA]\n");
    fprintf(stderr, "   -f, --fasta-ref <file>      reference sequence in fasta format\n");
    fprintf(stderr, "   -s, --samples <list>        list of sample names\n");
    fprintf(stderr, "   -S, --samples-file <file>   file of sample names\n");
    fprintf(stderr, "\n");
    // fprintf(stderr, "PLINK options:\n");
    // fprintf(stderr, "   -p, --plink <prefix>|<ped>,<map>,<fam>|<bed>,<bim>,<fam>|<tped>,<tfam>\n");
    // fprintf(stderr, "       --tped              make tped file instead\n");
    // fprintf(stderr, "       --bin               make binary bed/fam/bim files\n");
    // fprintf(stderr, "\n");
    // fprintf(stderr, "PBWT options:\n");
    // fprintf(stderr, "   -b, --pbwt          <prefix> or <pbwt>,<sites>,<sample>,<missing>\n");
    // fprintf(stderr, "\n");
    exit(1);
}

int main_vcfconvert(int argc, char *argv[])
{
    int c;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->outfname = "-";
    args->output_type = FT_VCF;

    static struct option loptions[] =
    {
        {"include",required_argument,NULL,'i'},
        {"exclude",required_argument,NULL,'e'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"targets",required_argument,NULL,'t'},
        {"targets-file",required_argument,NULL,'T'},
        {"samples",required_argument,NULL,'s'},
        {"samples-file",required_argument,NULL,'S'},
        {"gensample",required_argument,NULL,'g'},
        {"hapsample2vcf",required_argument,NULL,3},
        {"tag",required_argument,NULL,1},
        {"haplegendsample",required_argument,NULL,'h'},
        {"tsv2vcf",required_argument,NULL,2},
        {"columns",required_argument,NULL,'c'},
        {"fasta-ref",required_argument,NULL,'f'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "?h:r:R:s:S:t:T:i:e:g:o:O:c:f:",loptions,NULL)) >= 0) {
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
            case  2 : args->convert_func = tsv_to_vcf; args->infname = optarg; break;
            case  3 : args->convert_func = hapsample_to_vcf; args->infname = optarg; break;
            case 'f': args->ref_fname = optarg; break;
            case 'c': args->columns = optarg; break;
            case 'o': args->outfname = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                }
                break;
            case 'h': args->convert_func = vcf_to_haplegendsample; args->outfname = optarg; break;
            case '?': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( !args->infname )
    {
        if ( optind>=argc )
        {
            if ( !isatty(fileno((FILE *)stdin)) ) args->infname = "-";
        }
        else args->infname = argv[optind];
    }
    if ( !args->infname || !args->convert_func ) usage();

    args->convert_func(args);

    destroy_data(args);
    free(args);
    return 0;
}
