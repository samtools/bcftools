/* The MIT License

   Copyright (c) 2014 Genome Research Ltd.

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
   THE SOFTWARE.

 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <ctype.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kseq.h>
#include <htslib/regidx.h>
#include "bcftools.h"
#include "rbuf.h"

typedef struct
{
    kstring_t fa_buf;   // buffered reference sequence
    int fa_ori_pos;     // start position of the fa_buffer (wrt original sequence)
    int fa_frz_pos;     // protected position to avoid conflicting variants (last pos for SNPs/ins)
    int fa_mod_off;     // position difference of fa_frz_pos in the ori and modified sequence (ins positive)
    int fa_end_pos;     // region's end position in the original sequence
    int fa_case;        // output upper case or lower case?

    rbuf_t vcf_rbuf;
    bcf1_t **vcf_buf;
    int nvcf_buf, rid;

    regidx_t *mask;

    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    FILE *fp_out;
    char **argv;
    int argc, output_iupac, haplotype, isample;
    char *fname, *ref_fname, *sample, *output_fname, *mask_fname;
}
args_t;

static void init_data(args_t *args)
{
    args->files = bcf_sr_init();
    args->files->require_index = 1;
    if ( !bcf_sr_add_reader(args->files,args->fname) ) error("Failed to open %s: %s\n", args->fname, bcf_sr_strerror(args->files->errnum));
    args->hdr = args->files->readers[0].header;
    args->isample = -1;
    if ( args->sample )
    {
        args->isample = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,args->sample);
        if ( args->isample<0 ) error("No such sample: %s\n", args->sample);
    }
    if ( args->haplotype && args->isample<0 )
    {
        if ( bcf_hdr_nsamples(args->hdr) > 1 ) error("The --sample option is expected with --haplotype\n");
        args->isample = 0;
    }
    if ( args->mask_fname )
    {
        args->mask = regidx_init(args->mask_fname,NULL,NULL,0,NULL);
        if ( !args->mask ) error("Failed to initialize mask regions\n");
    }
    rbuf_init(&args->vcf_rbuf, 100);
    args->vcf_buf = (bcf1_t**) calloc(args->vcf_rbuf.m, sizeof(bcf1_t*));
    args->fp_out = args->output_fname ? fopen(args->output_fname,"w") : stdout;
}

static void destroy_data(args_t *args)
{
    bcf_sr_destroy(args->files);
    int i;
    for (i=0; i<args->vcf_rbuf.m; i++)
        if ( args->vcf_buf[i] ) bcf_destroy1(args->vcf_buf[i]);
    free(args->vcf_buf);
    free(args->fa_buf.s);
    if ( args->mask ) regidx_destroy(args->mask);
    if ( fclose(args->fp_out) ) error("Close failed: %s\n", args->output_fname);
}

static void init_region(args_t *args, char *line)
{
    char *ss, *se = line;
    while ( *se && !isspace(*se) && *se!=':' ) se++;
    int from = 0, to = 0x7fffffff;
    char tmp, *tmp_ptr = NULL;
    if ( *se )
    {
        tmp = *se; *se = 0; tmp_ptr = se;
        ss = ++se;
        from = strtol(ss,&se,10);
        if ( ss==se || !*se || *se!='-' ) from = 0;
        else
        {
            from--;
            ss = ++se;
            to = strtol(ss,&se,10);
            if ( ss==se || (*se && !isspace(*se)) ) { from = 0; to = 0x7fffffff; }
            else to--;
        }
    }
    args->rid = bcf_hdr_name2id(args->hdr,line);
    if ( args->rid<0 ) error("Sequence \"%s\" not in %s\n", line,args->fname);
    args->fa_buf.l = 0;
    args->fa_end_pos = to;
    args->fa_ori_pos = from;
    args->fa_mod_off = 0;
    args->fa_frz_pos = -1;
    args->fa_case    = -1;
    args->vcf_rbuf.n = 0;
    bcf_sr_seek(args->files,line,args->fa_ori_pos);
    if ( tmp_ptr ) *tmp_ptr = tmp;
    fprintf(args->fp_out,">%s\n",line);
}

static bcf1_t **next_vcf_line(args_t *args)
{
    if ( args->vcf_rbuf.n )
    {
        int i = rbuf_shift(&args->vcf_rbuf);
        return &args->vcf_buf[i];
    }
    else if ( bcf_sr_next_line(args->files) )
        return &args->files->readers[0].buffer[0];

    return NULL;
}
static void unread_vcf_line(args_t *args, bcf1_t **rec_ptr)
{
    bcf1_t *rec = *rec_ptr;
    if ( args->vcf_rbuf.n >= args->vcf_rbuf.m )
        error("FIXME: too many overlapping records near %s:%d\n", bcf_seqname(args->hdr,rec),rec->pos+1);

    // Insert the new record in the buffer. The line would be overwritten in
    // the next bcf_sr_next_line call, therefore we need to swap it with an
    // unused one
    int i = rbuf_append(&args->vcf_rbuf);
    if ( !args->vcf_buf[i] ) args->vcf_buf[i] = bcf_init1();
    bcf1_t *tmp = rec; *rec_ptr = args->vcf_buf[i]; args->vcf_buf[i] = tmp;
}
static void flush_fa_buffer(args_t *args, int len)
{
    if ( !args->fa_buf.l ) return;

    int nwr = 0;
    while ( nwr + 60 <= args->fa_buf.l )
    {
        if ( fwrite(args->fa_buf.s+nwr,1,60,args->fp_out) != 60 ) error("Could not write: %s\n", args->output_fname);
        if ( fwrite("\n",1,1,args->fp_out) != 1 ) error("Could not write: %s\n", args->output_fname);
        nwr += 60;
    }
    if ( nwr )
        args->fa_ori_pos += nwr;

    if ( len )
    {
        // not finished on this chr yet and the buffer cannot be emptied completely
        if ( nwr && nwr < args->fa_buf.l )
            memmove(args->fa_buf.s,args->fa_buf.s+nwr,args->fa_buf.l-nwr);
        args->fa_buf.l -= nwr;
        return;
    }

    // empty the whole buffer
    if ( nwr == args->fa_buf.l ) { args->fa_buf.l = 0; return; }

    if ( fwrite(args->fa_buf.s+nwr,1,args->fa_buf.l - nwr,args->fp_out) != args->fa_buf.l - nwr ) error("Could not write: %s\n", args->output_fname);
    if ( fwrite("\n",1,1,args->fp_out) != 1 ) error("Could not write: %s\n", args->output_fname);

    args->fa_ori_pos += args->fa_buf.l - nwr - args->fa_mod_off;
    args->fa_mod_off = 0;
    args->fa_buf.l = 0;
}
static void apply_variant(args_t *args, bcf1_t *rec)
{
    if ( rec->n_allele==1 ) return;

    if ( rec->pos <= args->fa_frz_pos )
    {
        if ( !args->mask )
            fprintf(stderr,"The site %s:%d overlaps with another variant, skipping...\n", bcf_seqname(args->hdr,rec),rec->pos+1);
        return;
    }

    int i, ialt = 1;
    if ( args->isample >= 0 )
    {
        bcf_fmt_t *fmt = bcf_get_fmt(args->hdr, rec, "GT");
        if ( !fmt ) return;
        if ( args->haplotype )
        {
            if ( args->haplotype > fmt->n ) error("Can't apply %d-th haplotype at %s:%d\n", args->haplotype,bcf_seqname(args->hdr,rec),rec->pos+1);
            uint8_t *ignore, *ptr = fmt->p + fmt->size*args->isample + args->haplotype - 1;
            ialt = bcf_dec_int1(ptr, fmt->type, &ignore);
            if ( bcf_gt_is_missing(ialt) || ialt==bcf_int32_vector_end ) return;
            ialt = bcf_gt_allele(ialt);
        }
        else if ( args->output_iupac ) 
        {
            uint8_t *ignore, *ptr = fmt->p + fmt->size*args->isample;
            ialt = bcf_dec_int1(ptr, fmt->type, &ignore);
            if ( bcf_gt_is_missing(ialt) || ialt==bcf_int32_vector_end ) return;
            ialt = bcf_gt_allele(ialt);

            int jalt;
            if ( fmt->n>1 )
            {
                ptr = fmt->p + fmt->size*args->isample + 1;
                jalt = bcf_dec_int1(ptr, fmt->type, &ignore);
                if ( bcf_gt_is_missing(jalt) || jalt==bcf_int32_vector_end ) jalt = ialt;
                else jalt = bcf_gt_allele(jalt);
            }
            else jalt = ialt;
            if ( rec->n_allele <= ialt || rec->n_allele <= jalt ) error("Broken VCF, too few alts at %s:%d\n", bcf_seqname(args->hdr,rec),rec->pos+1);
            if ( ialt!=jalt && !rec->d.allele[ialt][1] && !rec->d.allele[jalt][1] ) // is this a het snp?
            {
                char ial = rec->d.allele[ialt][0];
                char jal = rec->d.allele[jalt][0];
                rec->d.allele[ialt][0] = gt2iupac(ial,jal);
            }
        }
        else
        {
            for (i=0; i<fmt->n; i++)
            {
                uint8_t *ignore, *ptr = fmt->p + fmt->size*args->isample + i;
                ialt = bcf_dec_int1(ptr, fmt->type, &ignore);
                if ( bcf_gt_is_missing(ialt) || ialt==bcf_int32_vector_end ) return;
                ialt = bcf_gt_allele(ialt);
                if ( ialt ) break;
            }
        }
        if ( !ialt ) return;  // ref allele
        if ( rec->n_allele <= ialt ) error("Broken VCF, too few alts at %s:%d\n", bcf_seqname(args->hdr,rec),rec->pos+1);
    }
    else if ( args->output_iupac && !rec->d.allele[0][1] && !rec->d.allele[1][1] )
    {
        char ial = rec->d.allele[0][0];
        char jal = rec->d.allele[1][0];
        rec->d.allele[1][0] = gt2iupac(ial,jal);
    }

    int idx = rec->pos - args->fa_ori_pos + args->fa_mod_off;
    if ( idx<0 || idx>=args->fa_buf.l ) 
        error("FIXME: %s:%d .. idx=%d, ori_pos=%d, len=%d, off=%d\n",bcf_seqname(args->hdr,rec),rec->pos+1,idx,args->fa_ori_pos,args->fa_buf.l,args->fa_mod_off);

    // sanity check the reference base
    if ( strncasecmp(rec->d.allele[0],args->fa_buf.s+idx,rec->rlen) && strcasecmp(rec->d.allele[ialt], "<DEL>") )
    {
        // fprintf(stderr,"%d .. [%s], idx=%d ori=%d off=%d\n",args->fa_ori_pos,args->fa_buf.s,idx,args->fa_ori_pos,args->fa_mod_off);
        char tmp = 0;
        if ( args->fa_buf.l - idx > rec->rlen ) 
        { 
            tmp = args->fa_buf.s[idx+rec->rlen];
            args->fa_buf.s[idx+rec->rlen] = 0;
        }
        error(
            "The fasta sequence does not match the REF allele at %s:%d:\n"
            "   .vcf: [%s]\n" 
            "   .fa:  [%s]%c%s\n",
            bcf_seqname(args->hdr,rec),rec->pos+1, rec->d.allele[0], args->fa_buf.s+idx, 
            tmp?tmp:' ',tmp?args->fa_buf.s+idx+rec->rlen+1:""
            );
    }
    if ( !strcasecmp(rec->d.allele[ialt], "<DEL>") ) {
        rec->d.allele[ialt] = "\0";
    }

    int alen = strlen(rec->d.allele[ialt]);
    if ( args->fa_case )
        for (i=0; i<alen; i++) rec->d.allele[ialt][i] = toupper(rec->d.allele[ialt][i]);
    else
        for (i=0; i<alen; i++) rec->d.allele[ialt][i] = tolower(rec->d.allele[ialt][i]);

    int len_diff = alen - rec->rlen;
    if ( len_diff <= 0 )
    {
        // deletion or same size event
        for (i=0; i<alen; i++)
            args->fa_buf.s[idx+i] = rec->d.allele[ialt][i];
        if ( len_diff )
            memmove(args->fa_buf.s+idx+alen,args->fa_buf.s+idx+rec->rlen,args->fa_buf.l-idx-rec->rlen);
    }
    else
    {
        // insertion
        ks_resize(&args->fa_buf, args->fa_buf.l + len_diff);
        memmove(args->fa_buf.s + idx + rec->rlen + len_diff, args->fa_buf.s + idx + rec->rlen, args->fa_buf.l - idx - rec->rlen);
        for (i=0; i<alen; i++)
            args->fa_buf.s[idx+i] = rec->d.allele[ialt][i];
    }
    args->fa_buf.l += len_diff;
    args->fa_mod_off += len_diff;
    args->fa_frz_pos  = rec->pos + rec->rlen - 1;
}
static void mask_region(args_t *args, char *seq, int len)
{
    char *chr = (char*)bcf_hdr_id2name(args->hdr,args->rid);
    int start = args->fa_ori_pos + args->fa_buf.l;
    int end   = start + len - 1;

    regitr_t itr;
    if ( !regidx_overlap(args->mask, chr,start,end, &itr) ) return;

    int idx_start, idx_end, i;
    while ( REGITR_OVERLAP(itr,start,end) )
    {
        if ( args->fa_frz_pos<0 || args->fa_frz_pos < REGITR_END(itr) ) args->fa_frz_pos = REGITR_END(itr);
        idx_start = REGITR_START(itr) - args->fa_ori_pos;
        idx_end   = REGITR_END(itr) - args->fa_ori_pos;
        if ( idx_start < 0 ) idx_start = 0;
        if ( idx_end >= len ) idx_end = len - 1;
        for (i=idx_start; i<=idx_end; i++) seq[i] = 'N';
        itr.i++;
    }
}

static void consensus(args_t *args)
{
    htsFile *fasta = hts_open(args->ref_fname, "rb");
    if ( !fasta ) error("Error reading %s\n", args->ref_fname);
    kstring_t str = {0,0,0};
    while ( hts_getline(fasta, KS_SEP_LINE, &str) > 0 )
    {
        if ( str.s[0]=='>' )
        {
            while ( args->vcf_rbuf.n )
            {
                bcf1_t *rec = args->vcf_buf[args->vcf_rbuf.f];
                if ( rec->rid!=args->rid || rec->pos > args->fa_end_pos ) break;
                int i = rbuf_shift(&args->vcf_rbuf);
                apply_variant(args, args->vcf_buf[i]);
            }
            flush_fa_buffer(args, 0);
            init_region(args, str.s+1);
            continue;
        }

        // determine if uppercase or lowercase is used in this fasta file
        if ( args->fa_case==-1 ) args->fa_case = toupper(str.s[0])==str.s[0] ? 1 : 0;

        if ( args->mask ) mask_region(args, str.s, str.l);
        kputs(str.s, &args->fa_buf);

        bcf1_t **rec_ptr = NULL;
        while ( (rec_ptr = next_vcf_line(args)) )
        {
            bcf1_t *rec = *rec_ptr;

            // still the same chr and the same region? if not, fasta buf can be flushed
            if ( rec->rid!=args->rid || rec->pos > args->fa_end_pos )
            {
                // save the vcf record until next time and flush
                unread_vcf_line(args, rec_ptr);
                rec_ptr = NULL;
                break;
            }

            // is the vcf record well beyond cached fasta buffer? if yes, the buf can be flushed
            if ( args->fa_ori_pos + args->fa_buf.l - args->fa_mod_off <= rec->pos )
            {
                unread_vcf_line(args, rec_ptr);
                rec_ptr = NULL;
                break;
            }

            // is the cached fasta buffer full enough? if not, read more fasta, no flushing
            if ( args->fa_ori_pos + args->fa_buf.l - args->fa_mod_off <= rec->pos + rec->rlen )
            {
                unread_vcf_line(args, rec_ptr);
                break;
            }
            apply_variant(args, rec);
        }
        if ( !rec_ptr ) flush_fa_buffer(args, 60);
    }
    flush_fa_buffer(args, 0);
    hts_close(fasta);
    free(str.s);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Create consensus sequence by applying VCF variants to a reference\n");
    fprintf(stderr, "         fasta file.\n");
    fprintf(stderr, "Usage:   bcftools consensus [OPTIONS] <file.vcf>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -f, --fasta-ref <file>     reference sequence in fasta format\n");
    fprintf(stderr, "    -H, --haplotype <1|2>      apply variants for the given haplotype\n");
    fprintf(stderr, "    -i, --iupac-codes          output variants in the form of IUPAC ambiguity codes\n");
    fprintf(stderr, "    -m, --mask <file>          replace regions with N\n");
    fprintf(stderr, "    -o, --output <file>        write output to a file [standard output]\n");
    fprintf(stderr, "    -s, --sample <name>        apply variants of the given sample\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "   # Get the consensus for one region. The fasta header lines are then expected\n");
    fprintf(stderr, "   # in the form \">chr:from-to\".\n");
    fprintf(stderr, "   samtools faidx ref.fa 8:11870-11890 | bcftools consensus in.vcf.gz > out.fa\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_consensus(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;

    static struct option loptions[] = 
    {
        {"sample",1,0,'s'},
        {"iupac-codes",0,0,'i'},
        {"haplotype",1,0,'H'},
        {"output",1,0,'o'},
        {"fasta-ref",1,0,'f'},
        {"mask",1,0,'m'},
        {0,0,0,0}
    };
    char c;
    while ((c = getopt_long(argc, argv, "h?s:1iH:f:o:m:",loptions,NULL)) >= 0) 
    {
        switch (c) 
        {
            case 's': args->sample = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'i': args->output_iupac = 1; break;
            case 'f': args->ref_fname = optarg; break;
            case 'm': args->mask_fname = optarg; break;
            case 'H': 
                args->haplotype = optarg[0] - '0'; 
                if ( args->haplotype <=0 ) error("Expected positive integer with --haplotype\n");
                break;
            default: usage(args); break;
        }
    }
    if ( optind>=argc ) usage(args);
    args->fname = argv[optind];

    if ( !args->ref_fname && !isatty(fileno((FILE *)stdin)) ) args->ref_fname = "-";
    if ( !args->ref_fname ) usage(args);

    init_data(args);
    consensus(args);
    destroy_data(args);
    free(args);

    return 0;
}


