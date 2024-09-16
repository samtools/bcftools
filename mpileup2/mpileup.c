/*  mileup2.c -- mpileup v2.0 API

   Copyright (C) 2022-2023 Genome Research Ltd.

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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
//#include <htslib/regidx.h>    // fixme: itr_loop and itr->seq is not working for some reason
#include "regidx.h"
#include <sys/stat.h>
#include "mpileup.h"
#include "bam_sample.h"

typedef struct
{
    int icig;           // position in the cigar string
    int iseq;           // the cigar[icig] operation refers to seq[iseq]
    hts_pos_t ref_pos;  // ref coordinate of the first base (if the current icig op consumes ref)
}
cigar_state_t;

typedef struct
{
    bam1_t *bam;
    cigar_state_t cstate;
    uint32_t
        next:31,    // index of the next node, rbuf->tail's next points to itself
        is_set:1;   // if is_set==0 && next==0, the node is unitialized in a newly allocated block
}
rb_node_t;

// Reads from a small window around the current position, per sample, not per BAM
typedef struct
{
    int tid;
    hts_pos_t pos;  // current position
    int len;        // keep reads overlapping pos+len
    rb_node_t *buf; // the reads, unordered and with gaps, linked through 'next' indexes
    int n, m;       // the number of buffer elements used and allocated
    int head, tail; // the index of the first and last used record
    int empty;      // the index of the first unused record
    int nins, ndel, nbase;
}
read_buffer_t;

// Data associated with one bam file
typedef struct
{
    samFile *fp;
    bam_hdr_t *hdr;
    hts_idx_t *idx;
    hts_itr_t *iter;
    char *fname;
    int bam_id;
    bam1_t *cached_rec;
    int cached_smpl;        // short-term memory to remember read's sample id across function calls (unused)
    int32_t tid;
    hts_pos_t pos;
    mpileup_t *mplp;        // make the core struct accessible in legacy API
}
bam_aux_t;

typedef struct
{
    int active;
    bam_mplp_t iter;
    bam_pileup1_t **plp;
    int *nplp;
}
legacy_t;

typedef struct
{
    char *seq;
    int tid, len;
}
refseq_t;

struct mpileup_t_
{
    char *fasta_ref, *targets, *targets_fname, *samples, *samples_fname, *read_groups_fname;
    int smart_overlaps, smpl_ignore_rg, max_dp_per_sample, adjust_mq, min_mq, min_bq, max_bq, delta_bq;
    int min_realn_dp, max_realn_dp, max_realn_len, skip_any_unset, skip_all_unset, skip_any_set, skip_all_set;
    float min_realn_frac;
    int nbam, nbam_names;
    char **bam_names;
    regidx_t *regions_idx, *targets_idx;
    regitr_t *regions_itr;
    int nregions;
    bam_aux_t *bam;
    bam_hdr_t *hdr;
    bam_smpl_t *bsmpl;
    read_buffer_t *read_buf;    // read buffer for each sample, nsmpl elements
    int nsmpl;
    int32_t tid;
    char *chr;
    hts_pos_t pos;
    kstring_t tmp_str;
    legacy_t legacy;
    refseq_t ref[2];
    faidx_t *fai;
};

static int mplp_add_bam(mpileup_t *mplp, char *bam, int is_file);
static int mplp_add_region(mpileup_t *mplp, char *regions, int is_file);
static int mplp_add_target(mpileup_t *mplp, char *regions, int is_file);
static void read_buffer_clear(read_buffer_t *rbuf);
static void read_buffer_reset(read_buffer_t *rbuf);
static char *mplp_get_ref(mpileup_t *mplp, int *len);

static int legacy_mpileup_init(mpileup_t *mplp);
static int legacy_mpileup_destroy(mpileup_t *mplp);
static int legacy_mpileup_next(mpileup_t *mplp);

mpileup_t *mpileup_alloc(void)
{
    mpileup_t *mplp = (mpileup_t*) calloc(1,sizeof(mpileup_t));
    mpileup_set(mplp,MAX_DP_PER_SAMPLE,250);
    mpileup_set(mplp,MIN_BQ,1);
    mpileup_set(mplp,MAX_BQ,60);
    mpileup_set(mplp,DELTA_BQ,30);
    mpileup_set(mplp,MIN_REALN_FRAC,0.05);
    mpileup_set(mplp,MIN_REALN_DP,2);
    mpileup_set(mplp,MAX_REALN_DP,250);
    mpileup_set(mplp,MAX_REALN_LEN,500);
    mpileup_set(mplp,SKIP_ANY_SET,BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP);
    mplp->pos = -1;
    mplp->tid = -1;
    mplp->ref[0].tid = -1;
    mplp->ref[1].tid = -1;
    return mplp;
}

int mpileup_set(mpileup_t *mplp, mpileup_opt_t key, ...)
{
    va_list args;
    switch (key)
    {
        case LEGACY_MODE:
            va_start(args, key);
            mplp->legacy.active = va_arg(args,int);
            va_end(args);
            return 0;

        case FASTA_REF:
            free(mplp->fasta_ref);
            va_start(args, key);
            mplp->fasta_ref = strdup(va_arg(args,char*));
            va_end(args);
            if ( !mplp->fasta_ref ) return -1;
            if ( !(mplp->fai=fai_load(mplp->fasta_ref)) )
            {
                free(mplp->fasta_ref);
                mplp->fasta_ref = NULL;
                return -1;
            }
            return 0;

        case BAM:
        {
            va_start(args, key);
            char *bam = va_arg(args,char*);
            va_end(args);
            return mplp_add_bam(mplp,bam,0);
        }

        case BAM_FNAME:
        {
            va_start(args, key);
            char *bam = va_arg(args,char*);
            va_end(args);
            return mplp_add_bam(mplp,bam,1);
        }

        case REGIONS:
        {
            va_start(args, key);
            char *regions = va_arg(args,char*);
            va_end(args);
            return mplp_add_region(mplp,regions,0);
        }

        case REGIONS_FNAME:
        {
            va_start(args, key);
            char *regions = va_arg(args,char*);
            va_end(args);
            return mplp_add_region(mplp,regions,1);
        }

        case TARGETS:
        {
            va_start(args, key);
            char *targets = va_arg(args,char*);
            va_end(args);
            return mplp_add_target(mplp,targets,0);
        }

        case TARGETS_FNAME:
        {
            va_start(args, key);
            char *targets = va_arg(args,char*);
            va_end(args);
            return mplp_add_target(mplp,targets,1);
        }

        case SAMPLES:
            free(mplp->samples);
            va_start(args, key);
            mplp->samples = strdup(va_arg(args,char*));
            va_end(args);
            return mplp->samples ? 0 : -1;

        case SAMPLES_FNAME:
            free(mplp->samples_fname);
            va_start(args, key);
            mplp->samples_fname = strdup(va_arg(args,char*));
            va_end(args);
            return mplp->samples_fname ? 0 : -1;

        case READ_GROUPS_FNAME:
            free(mplp->read_groups_fname);
            va_start(args, key);
            mplp->read_groups_fname = strdup(va_arg(args,char*));
            va_end(args);
            return mplp->read_groups_fname ? 0 : -1;

        case SMART_OVERLAPS:
            va_start(args, key);
            mplp->smart_overlaps = va_arg(args,int);
            va_end(args);
            return 0;

        case SMPL_IGNORE_RG:
            va_start(args, key);
            mplp->smpl_ignore_rg = va_arg(args,int);
            va_end(args);
            return 0;

        case SKIP_ANY_UNSET:
            va_start(args, key);
            mplp->skip_any_unset = va_arg(args,int);
            va_end(args);
            return 0;

        case SKIP_ALL_UNSET:
            va_start(args, key);
            mplp->skip_all_unset = va_arg(args,int);
            va_end(args);
            return 0;

        case SKIP_ANY_SET:
            va_start(args, key);
            mplp->skip_any_set = va_arg(args,int);
            va_end(args);
            return 0;

        case SKIP_ALL_SET:
            va_start(args, key);
            mplp->skip_all_set = va_arg(args,int);
            va_end(args);
            return 0;

        case MAX_DP_PER_SAMPLE:
            va_start(args, key);
            mplp->max_dp_per_sample = va_arg(args,int);
            va_end(args);
            return 0;

        case ADJUST_MQ:
            va_start(args, key);
            mplp->adjust_mq = va_arg(args,int);
            va_end(args);
            return 0;

        case MIN_MQ:
            va_start(args, key);
            mplp->min_mq = va_arg(args,int);
            va_end(args);
            return 0;

        case MIN_BQ:
            va_start(args, key);
            mplp->min_bq = va_arg(args,int);
            va_end(args);
            return 0;

        case MAX_BQ:
            va_start(args, key);
            mplp->max_bq = va_arg(args,int);
            va_end(args);
            return 0;

        case DELTA_BQ:
            va_start(args, key);
            mplp->delta_bq = va_arg(args,int);
            va_end(args);
            return 0;

        case MIN_REALN_FRAC:
            va_start(args, key);
            mplp->min_realn_frac = va_arg(args,double);
            va_end(args);
            return 0;

        case MIN_REALN_DP:
            va_start(args, key);
            mplp->min_realn_dp = va_arg(args,int);
            va_end(args);
            return 0;

        case MAX_REALN_DP:
            va_start(args, key);
            mplp->max_realn_dp = va_arg(args,int);
            va_end(args);
            return 0;

        case MAX_REALN_LEN:
            va_start(args, key);
            mplp->max_realn_len = va_arg(args,int);
            va_end(args);
            return 0;

        default:
            hts_log_error("Todo: mpileup_set key=%d",(int)key);
            return -1;
            break;
    }
    return 0;
}

void *mpileup_get(mpileup_t *mplp, mpileup_opt_t key, ...)
{
    va_list args;
    switch (key)
    {
        case SMART_OVERLAPS: return &mplp->smart_overlaps; break;
        case SMPL_IGNORE_RG: return &mplp->smpl_ignore_rg; break;
        case MAX_DP_PER_SAMPLE: return &mplp->max_dp_per_sample; break;
        case ADJUST_MQ: return &mplp->adjust_mq; break;
        case MIN_MQ: return &mplp->min_mq; break;
        case MIN_BQ: return &mplp->min_bq; break;
        case MAX_BQ: return &mplp->max_bq; break;
        case DELTA_BQ: return &mplp->delta_bq; break;
        case MIN_REALN_FRAC: return &mplp->min_realn_frac; break;
        case MIN_REALN_DP: return &mplp->min_realn_dp; break;
        case MAX_REALN_DP: return &mplp->max_realn_dp; break;
        case MAX_REALN_LEN: return &mplp->max_realn_len; break;
        case SKIP_ANY_UNSET: return &mplp->skip_any_unset; break;
        case SKIP_ALL_UNSET: return &mplp->skip_all_unset; break;
        case SKIP_ANY_SET: return &mplp->skip_any_set; break;
        case SKIP_ALL_SET: return &mplp->skip_all_set; break;
        case N_SAMPLES: return &mplp->nsmpl; break;
        case N_READS: return &mplp->legacy.nplp; break;
        case CHROM: return &mplp->chr; break;
        case POS: return &mplp->pos; break;
        case LEGACY_PILEUP: return &mplp->legacy.plp; break;
        case REF:
            va_start(args, key);
            int *len_ptr = va_arg(args,int*);
            va_end(args);
            return mplp_get_ref(mplp, len_ptr);
            break;
        default: hts_log_error("Todo: mpileup_get key=%d",(int)key); return NULL; break;
    }
    return NULL;
}

void mpileup_destroy(mpileup_t *mplp)
{
    legacy_mpileup_destroy(mplp);
    int i;
    for (i=0; i<mplp->nbam; i++)
    {
        if ( mplp->bam[i].iter ) hts_itr_destroy(mplp->bam[i].iter);
        sam_close(mplp->bam[i].fp);
        if ( mplp->bam[i].cached_rec ) bam_destroy1(mplp->bam[i].cached_rec);
    }
    for (i=0; i<mplp->nbam_names; i++) free(mplp->bam_names[i]);
    free(mplp->bam);
    free(mplp->bam_names);
    for (i=0; i<mplp->nsmpl; i++) read_buffer_clear(&mplp->read_buf[i]);
    free(mplp->read_buf);
    bam_hdr_destroy(mplp->hdr);
    if ( mplp->bsmpl ) bam_smpl_destroy(mplp->bsmpl);
    free(mplp->fasta_ref);
    free(mplp->targets);
    free(mplp->targets_fname);
    if ( mplp->targets_idx ) regidx_destroy(mplp->targets_idx);
    if ( mplp->regions_idx ) regidx_destroy(mplp->regions_idx);
    if ( mplp->regions_itr ) regitr_destroy(mplp->regions_itr);
    free(mplp->samples);
    free(mplp->samples_fname);
    free(mplp->read_groups_fname);
    free(mplp->tmp_str.s);
    free(mplp->ref[0].seq);
    free(mplp->ref[1].seq);
    if ( mplp->fai ) fai_destroy(mplp->fai);
    free(mplp);
}

static int is_url(const char *str)
{
    static const char uri_scheme_chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+.-";
    return str[strspn(str, uri_scheme_chars)] == ':';
}

static int mplp_add_target(mpileup_t *mplp, char *regions, int is_file)
{
    if ( mplp->targets_idx ) return -1;
    if ( is_file )
        mplp->targets_idx = regidx_init(regions,NULL,NULL,0,NULL);
    else
    {
        mplp->targets_idx = regidx_init(NULL,regidx_parse_reg,NULL,sizeof(char*),NULL);
        if ( !mplp->targets_idx || regidx_insert_list(mplp->targets_idx,regions,',')!=0  )
        {
            hts_log_error("Could not initialize the targets: %s",regions);
            return -1;
        }
    }
    if ( !mplp->targets_idx )
    {
        hts_log_error("Could not initialize the targets: %s",regions);
        return -1;
    }
    return 0;
}
static int mplp_add_region(mpileup_t *mplp, char *regions, int is_file)
{
    if ( mplp->regions_idx ) return -1;     // regions can be added only once
    if ( is_file )
        mplp->regions_idx = regidx_init(regions,NULL,NULL,0,NULL);
    else
    {
        mplp->regions_idx = regidx_init(NULL,regidx_parse_reg,NULL,sizeof(char*),NULL);
        if ( !mplp->regions_idx || regidx_insert_list(mplp->regions_idx,regions,',')!=0  )
        {
            hts_log_error("Could not initialize the regions: %s",regions);
            return -1;
        }
    }
    if ( !mplp->regions_idx )
    {
        hts_log_error("Could not initialize the regions: %s",regions);
        return -1;
    }
    mplp->nregions = regidx_nregs(mplp->regions_idx);
    mplp->regions_itr = regitr_init(mplp->regions_idx);
    regitr_loop(mplp->regions_itr);   // region iterator now positioned at the first region
    return 0;
}
static int mplp_add_bam(mpileup_t *mplp, char *bam, int is_file)
{
    struct stat sbuf;
    if ( is_file )
    {
        int i, nnames = 0;
        char **names = hts_readlist(bam, 1, &nnames);
        for (i=0; i<nnames; i++)
        {
            if ( !(is_url(names[i]) || stat(names[i],&sbuf)==0) )
            {
                // no such file, check if it is safe to print its name
                char *tmp = names[i];
                while ( *tmp )
                {
                    if ( !isprint(*tmp) ) break;
                    tmp++;
                }
                if ( !*tmp )
                    hts_log_error("The file list \"%s\" appears broken, could not locate file: \"%s\"", bam,names[i]);
                else
                    hts_log_error("Does the file \"%s\" really contain a list of files and do all exist?", bam);
                break;
            }
            char **tmp = (char**)realloc(mplp->bam_names,(mplp->nbam_names+1)*sizeof(*mplp->bam_names));
            if ( !tmp )
            {
                hts_log_error("Could not allocate %zu bytes",(mplp->nbam_names+1)*sizeof(*mplp->bam_names));
                break;
            }
            mplp->bam_names = tmp;
            mplp->bam_names[mplp->nbam_names++] = strdup(names[i]);
        }
        for (i=0; i<nnames; i++) free(names[i]);
        free(names);
        return i==nnames ? 0 : -1;
    }
    char **tmp = (char**)realloc(mplp->bam_names,(mplp->nbam_names+1)*sizeof(*mplp->bam_names));
    if ( !tmp )
    {
        hts_log_error("Could not allocate %zu bytes",(mplp->nbam_names+1)*sizeof(*mplp->bam_names));
        return -1;
    }
    mplp->bam_names = tmp;
    mplp->bam_names[mplp->nbam_names++] = strdup(bam);
    return 0;
}

static int mplp_check_header_contigs(mpileup_t *mplp, bam_hdr_t *hdr)
{
    static int warned = 0;
    if ( !warned ) hts_log_error("todo: mplp_check_header_contigs");
    warned = 1;
    return 0;
}
int mpileup_init(mpileup_t *mplp)
{
    mplp->bsmpl = bam_smpl_init();
    if ( mplp->smpl_ignore_rg ) bam_smpl_ignore_readgroups(mplp->bsmpl);
    if ( mplp->read_groups_fname ) bam_smpl_add_readgroups(mplp->bsmpl, mplp->read_groups_fname, 1);
    if ( mplp->samples && !bam_smpl_add_samples(mplp->bsmpl,mplp->samples,0) )
    {
        hts_log_error("Could not parse samples: %s", mplp->samples);
        return -1;
    }
    if ( mplp->samples_fname && !bam_smpl_add_samples(mplp->bsmpl,mplp->samples_fname,0) )
    {
        hts_log_error("Could not read samples: %s", mplp->samples_fname);
        return -1;
    }

    int nregs = 0;
    char **regs = NULL;
    if ( mplp->regions_idx )
    {
        nregs = regidx_nregs(mplp->regions_idx);
        regs = calloc(nregs,sizeof(*regs));
        int i = 0;
        regitr_reset(mplp->regions_idx,mplp->regions_itr);
        while ( regitr_loop(mplp->regions_itr) )
        {
            assert(i < nregs);
            mplp->tmp_str.l = 0;
            ksprintf(&mplp->tmp_str,"%s:%u-%u",mplp->regions_itr->seq,mplp->regions_itr->beg+1,mplp->regions_itr->end+1);
            regs[i++] = strdup(mplp->tmp_str.s);
        }
    }

    mplp->bam = (bam_aux_t*)calloc(mplp->nbam_names,sizeof(*mplp->bam));
    if ( !mplp->bam )
    {
        hts_log_error("Could not allocate %zu bytes",mplp->nbam_names*sizeof(*mplp->bam));
        return -1;
    }
    int i;
    for (i=0; i<mplp->nbam_names; i++)
    {
        char *fname = mplp->bam_names[i];
        bam_aux_t *baux = &mplp->bam[mplp->nbam];
        baux->mplp = mplp;
        baux->cached_smpl = -1;
        baux->fname = fname;
        baux->fp = sam_open(fname,"rb");
        if ( !baux->fp )
        {
            hts_log_error("Failed to open %s: %s",fname,strerror(errno));
            return -1;
        }
        if ( hts_set_opt(baux->fp, CRAM_OPT_DECODE_MD, 1) )
        {
            hts_log_error("Failed to set CRAM_OPT_DECODE_MD value: %s\n",fname);
            return -1;
        }
        if ( mplp->fasta_ref && hts_set_fai_filename(baux->fp,mplp->fasta_ref)!=0 )
        {
            hts_log_error("Failed to process %s: %s\n",mplp->fasta_ref, strerror(errno));
            return -1;
        }
  // baux->conf = conf;
  // baux->ref = &mp_ref;
        bam_hdr_t *hdr = sam_hdr_read(baux->fp);
        if ( !hdr )
        {
            hts_log_error("Failed to read the header of %s\n",fname);
            return -1;
        }
        baux->bam_id = bam_smpl_add_bam(mplp->bsmpl,hdr->text,fname);
        if ( baux->bam_id < 0 )
        {
            // no usable read groups in this bam, can be skipped
            sam_close(baux->fp);
            bam_hdr_destroy(hdr);
        }
        mplp->nbam++;

        // iterators
        if ( nregs )
        {
            hts_idx_t *idx = sam_index_load(baux->fp,baux->fname);
            if ( !idx )
            {
                hts_log_error("Failed to load index: %s\n",fname);
                return 1;
            }
            baux->iter = sam_itr_regarray(idx,hdr,regs,nregs);
            if ( !baux->iter )
            {
                hts_log_error("Failed to initialize %d regions",nregs);
                return 1;
            }
            if ( mplp->nregions==1 ) // no need to keep the index in memory
                hts_idx_destroy(idx);
            else
                baux->idx = idx;
        }

        if ( !mplp->hdr ) mplp->hdr = hdr;
        else
        {
            int ret = mplp_check_header_contigs(mplp,hdr);
            bam_hdr_destroy(hdr);
            if ( ret!=0 ) return -1;
        }
        baux->hdr = mplp->hdr;
        baux->tid = -1;
    }
    for (i=0; i<nregs; i++) free(regs[i]);
    free(regs);
    bam_smpl_get_samples(mplp->bsmpl,&mplp->nsmpl);
    if ( !mplp->nsmpl )
    {
        hts_log_error("Failed to initialize samples");
        return -1;
    }
    mplp->read_buf = (read_buffer_t*) calloc(mplp->nsmpl,sizeof(*mplp->read_buf));

    if ( mplp->legacy.active ) return legacy_mpileup_init(mplp);

    return 0;
}


// Clear the per-sample buffer, destroying all its data
static void read_buffer_clear(read_buffer_t *rbuf)
{
    int i;
    for (i=0; i<rbuf->m; i++)
    {
        if ( rbuf->buf[i].bam ) bam_destroy1(rbuf->buf[i].bam);
    }
    free(rbuf->buf);
    memset(rbuf,0,sizeof(*rbuf));
}

// Reset the buffer to be ready for the next region
static void read_buffer_reset(read_buffer_t *rbuf)
{
    int i;
    for (i=0; i<rbuf->m; i++) rbuf->buf[i].is_set = 0;
    rbuf->n = 0;
    rbuf->tid = 0;
    rbuf->pos = 0;
    rbuf->empty = 0;
}

/*
 * Returns
 *   0 .. on success (pushed onto the buffer, swapped rec for NULL or an empty record)
 *   1 .. cannot push (different tid or not ready for the current position yet)
 *   2 .. cannot push (past the current position)
 *  -1 .. error
 */
static int read_buffer_push(read_buffer_t *rbuf, hts_pos_t pos, bam1_t **rec_ptr)
{
static int warned = 0;
if ( !warned ) hts_log_error("read_buffer_push - left softclip vs read start");
warned = 1;

    bam1_t *rec = *rec_ptr;
    if ( rec->core.tid != rbuf->tid ) return 1;
    if ( rbuf->len < rec->core.l_qseq ) rbuf->len = rec->core.l_qseq;
    if ( rec->core.pos + rbuf->len < pos ) return 2;
    if ( rec->core.pos > pos + rbuf->len ) return 2;
//fprintf(stderr,"pushing %d:  %d+%d  in %d-%d\n",rbuf->n,(int)rec->core.pos,(int)rec->core.l_qseq,(int)pos,rbuf->len);
    if ( rbuf->n==rbuf->m )
    {
        // Need to increase the buffer size
        if ( hts_resize(rb_node_t, rbuf->n + 1, &rbuf->m, &rbuf->buf, HTS_RESIZE_CLEAR)<0 )
        {
            hts_log_error("Failed to allocate memory");
            return -1;
        }
        rbuf->empty = rbuf->n;  // the first newly allocated block is the new empty
    }
    // claim an empty node and populate it
    rb_node_t *node = &rbuf->buf[rbuf->empty];
    *rec_ptr   = node->bam;
    node->bam  = rec;
    node->cstate.ref_pos = rec->core.pos;
    node->next = rbuf->empty;   // the tail's next points to itself
    if ( rbuf->n )
        rbuf->buf[rbuf->tail].next = rbuf->empty;
    else
        rbuf->head = rbuf->empty;
    rbuf->tail = rbuf->empty;

    // Assign the head of empty nodes. Note that the latter can go beyond the allocated buffer or nonsense value,
    // but that will be set above to a correct value later when the buffer size is increased
    rbuf->empty  = node->is_set ? node->next : rbuf->empty + 1;
    node->is_set = 1;
    rbuf->n++;

    if ( pos==0 && rbuf->n==1 )
    {
        rbuf->tid = rbuf->buf[rbuf->head].bam->core.tid;
        rbuf->pos = rbuf->buf[rbuf->head].bam->core.pos;
    }

    return 0;
}

/*
 * Read from a bam and fill a read buffer (separate for each sample) to allow read sub-sampling.
 *
 * Returns 0 on success (buffers full), 1 when done (i.e. no new reads), or a negative value on error
 */
static int mplp_read_bam(mpileup_t *mplp, bam_aux_t *baux)
{
    // if new mplp position is requested (different tid or mplp pos outside of the last region)
    //    update the sam iterator baux->iter = sam_itr_querys(conf->mplp_data[i]->idx, conf->mplp_data[i]->h, conf->buf.s);

    if ( baux->cached_smpl >= 0 )
    {
        read_buffer_t *rbuf = &mplp->read_buf[baux->cached_smpl];
        int ret = read_buffer_push(rbuf, mplp->pos, &baux->cached_rec);
        if ( ret<0 ) return ret;      // an error (malloc)
        if ( ret==1 ) return 0;       // cannot push, full buffer, not ready for the current position
        baux->cached_smpl = -1;
    }
    while (1)
    {
        if ( !baux->cached_rec ) baux->cached_rec = bam_init1();
        int ret = baux->iter ? sam_itr_next(baux->fp, baux->iter, baux->cached_rec) : sam_read1(baux->fp, baux->hdr, baux->cached_rec);
        if ( ret<0 )
        {
            if ( ret==-1 ) return 1;    // no more data
            return -1;                  // error
        }

        // Check if the records have tids in ascending order. This is not required by the SAM specification but
        // mpileup2 relies on it (for now)
        bam1_t *rec = baux->cached_rec;
        if ( baux->tid==-1 ) baux->tid = rec->core.tid, baux->pos = rec->core.pos;
        else if ( baux->tid > rec->core.tid )
        {
            hts_log_error("Failed assumption, sequence tids out of order in %s",baux->fname);
            return -1;
        }
        else if ( baux->pos > rec->core.pos )
        {
            hts_log_error("Unsorted alignments in %s",baux->fname);
            return -1;
        }
        baux->pos = rec->core.pos;

        // exclude unmapped reads and any filters provided by the user
        if ( rec->core.tid < 0 || (rec->core.flag&BAM_FUNMAP) ) continue;
        if ( mplp->skip_any_unset && (mplp->skip_any_unset & rec->core.flag)!=mplp->skip_any_unset ) continue;
        if ( mplp->skip_all_set && (mplp->skip_all_set & rec->core.flag)==mplp->skip_all_set ) continue;
        if ( mplp->skip_all_unset && !(mplp->skip_all_unset & rec->core.flag) ) continue;
        if ( mplp->skip_any_set && (mplp->skip_any_set & rec->core.flag) ) continue;

        // exclude if sample is not wanted
        baux->cached_smpl = bam_smpl_get_sample_id(mplp->bsmpl, baux->bam_id, rec);
        if ( baux->cached_smpl<0 ) continue;

        read_buffer_t *rbuf = &mplp->read_buf[baux->cached_smpl];
        ret = read_buffer_push(rbuf, mplp->pos, &baux->cached_rec);
        if ( ret==1 ) return 0;     // cannot push, full buffer, not ready for the current position
        if ( ret<0 ) return ret;    // an error (malloc)
        baux->cached_smpl = -1;
    }
    return -1;  // can never happen
}

/*
 *  Returns
 *      BAM_CMATCH when read overlaps pos and has a base match or mismatch. Sets ipos
 *      BAM_CINS when the next base is an insertion. Sets ipos to the current base
 *      BAM_CDEL when the next or overlapping base is a deletion. Sets ipos to the current
 *              base when the next base is a deletion and to -1 when inside the deletion
 *      -1 when there is no base to show (e.g. soft clip)
 *      -2 when the read ends before pos
 */
static int cigar_get_pos(cigar_state_t *cs, bam1_t *bam, hts_pos_t pos, int32_t *ipos)
{
    int ncig = bam->core.n_cigar;
    uint32_t *cigar = bam_get_cigar(bam);
    while ( cs->ref_pos <= pos )
    {
        if ( cs->icig >= ncig ) return -2;    // the read ends before pos
        int op  = cigar[cs->icig] &  BAM_CIGAR_MASK;
        int len = cigar[cs->icig] >> BAM_CIGAR_SHIFT;
        if ( op==BAM_CMATCH || op==BAM_CEQUAL || op==BAM_CDIFF )
        {
            hts_pos_t end_pos = cs->ref_pos + len - 1;
            if ( end_pos < pos )    // pos not covered by this cigar block
            {
                cs->ref_pos += len;
                cs->iseq += len;
                cs->icig++;
                continue;
            }
            *ipos = pos - cs->ref_pos + cs->iseq;

            // if the following base is an insertion, return BAM_CINS
            if ( end_pos==pos && cs->icig+1 < ncig )
            {
                int next_op = cigar[cs->icig+1] & BAM_CIGAR_MASK;
                if ( next_op==BAM_CINS ) return BAM_CINS;
                if ( next_op==BAM_CDEL ) return BAM_CDEL;
            }
            return BAM_CMATCH;
        }
        if ( op==BAM_CINS || op==BAM_CSOFT_CLIP )
        {
            cs->iseq += len;
            cs->icig++;
            continue;
        }
        if ( op==BAM_CDEL )
        {
            hts_pos_t end_pos = cs->ref_pos + len - 1;
            if ( end_pos < pos )    // pos not covered by this cigar block
            {
                cs->ref_pos += len;
                cs->icig++;
                continue;
            }
            *ipos = -1;
            return BAM_CDEL;
        }
        if ( op==BAM_CREF_SKIP )
        {
            hts_pos_t end_pos = cs->ref_pos + len - 1;
            if ( end_pos < pos )    // pos not covered by this cigar block
            {
                cs->ref_pos += len;
                cs->icig++;
                continue;
            }
            return -1;
        }
        hts_log_error("todo: CIGAR operator %d", op);
        assert(0);
    }
    return -1;
}

// Go through all reads in the buffer and make the pileup structures ready for output
static int mplp_set_pileup(mpileup_t *mplp, read_buffer_t *rbuf)
{
    if ( !rbuf->n ) return 0;
    rbuf->ndel = rbuf->nins = rbuf->nbase = 0;
    int iprev = -1, i = rbuf->head;
    while ( 1 )
    {
        rb_node_t *node = &rbuf->buf[i];    // rb_node_t keeps info about a single read
        int32_t ipos;
        int ret = cigar_get_pos(&node->cstate, node->bam, mplp->pos, &ipos);
        if ( ret==BAM_CMATCH )   // there is a base, match or mismatch
        {
            uint8_t *seq = bam_get_seq(node->bam);
            fprintf(stderr," %c","ACGTN"[seq_nt16_int[bam_seqi(seq,ipos)]]);
            rbuf->nbase++;
        }
        else if ( ret==BAM_CINS ) rbuf->nins++;
        else if ( ret==BAM_CDEL && ipos!=-1 ) rbuf->ndel++;
        else if ( ret==-2 )     // this read is done, remove the node
        {
            rbuf->n--;
            node->is_set = 0;
            int next = node->next;
            node->next = rbuf->empty;
            rbuf->empty = i;
            if ( next==i ) // node that points to itself is the last node in the chain
            {
                assert( rbuf->tail==i );
                rbuf->tail = iprev;
                if ( iprev==-1 ) rbuf->head = iprev;
                else rbuf->buf[iprev].next = iprev;
                break;
            }
            else if ( iprev==-1 )   // removing the first node
                rbuf->head = next;
            else                    // removing a middle node
                rbuf->buf[iprev].next = next;
            i = next;
            continue;
        }

        if ( rbuf->buf[i].next == i ) break;
        iprev = i;
        i = rbuf->buf[i].next;
    }
    fprintf(stderr,"\n");
    return 0;
}

static int mplp_region_next(mpileup_t *mplp)
{
    int i;

    mplp->pos++;
    if ( mplp->nregions && mplp->pos > mplp->regions_itr->end ) return 0;

    // fill buffers
    for (i=0; i<mplp->nbam; i++)
    {
        int ret = mplp_read_bam(mplp, &mplp->bam[i]);
        if ( ret<0 ) return ret;    // an error occurred
    }

    // this would be a good place for realignment

    // find the minimum position if first time here
    if ( mplp->tid==-1 )
    {
        mplp->tid = INT32_MAX;
        for (i=0; i<mplp->nsmpl; i++)
        {
            if ( mplp->tid > mplp->read_buf[i].tid ) mplp->tid = mplp->read_buf[i].tid, mplp->pos = HTS_POS_MAX;
            if ( mplp->pos > mplp->read_buf[i].pos ) mplp->pos = mplp->read_buf[i].pos;
        }
        if ( mplp->nregions && mplp->pos < mplp->regions_itr->beg ) mplp->pos = mplp->regions_itr->beg;
    }

    // prepare the pileup
    int nreads = 0;
    for (i=0; i<mplp->nsmpl; i++)
    {
        int ret = mplp_set_pileup(mplp, &mplp->read_buf[i]);
        if ( ret!=0 ) return -1;
        nreads += mplp->read_buf[i].n;
    }
    return nreads>0 /* || nregions */ ? 1 : 0;
}

int mpileup_next(mpileup_t *mplp)
{
    if ( mplp->legacy.active ) return legacy_mpileup_next(mplp);

    // if new region is needed
    //  o advance regitr_loop
    //  o reset read_buffers
    //  o reset mplp tid and pos

    if ( !mplp->nregions )
        return mplp_region_next(mplp);

    int i;
    while (1)
    {
        int ret = mplp_region_next(mplp);
        if ( ret>0 )
        {
            if ( mplp->pos < mplp->regions_itr->beg ) continue;
            if ( mplp->pos <= mplp->regions_itr->end ) return ret;
        }
        if ( !regitr_loop(mplp->regions_itr) ) return 0;

        for (i=0; i<mplp->nsmpl; i++) read_buffer_reset(&mplp->read_buf[i]);

        mplp->tmp_str.l = 0;
        //ksprintf(&mplp->tmp_str,"%s:%"PRIhts_pos"-%"PRIhts_pos,mplp->regions_itr->seq,mplp->regions_itr->beg+1,mplp->regions_itr->end+1);
        ksprintf(&mplp->tmp_str,"%s:%u-%u",mplp->regions_itr->seq,mplp->regions_itr->beg+1,mplp->regions_itr->end+1);
        for (i=0; i<mplp->nbam; i++)
        {
            bam_aux_t *baux = &mplp->bam[i];
            hts_itr_destroy(baux->iter);
            baux->iter = sam_itr_querys(baux->idx, baux->hdr, mplp->tmp_str.s);
            if ( !baux->iter )
            {
                baux->iter = sam_itr_querys(baux->idx, baux->hdr, mplp->regions_itr->seq);
                if ( baux->iter )
                {
                    hts_log_error("[Failed to parse the region: %s",mplp->tmp_str.s);
                    return 1;
                }
                hts_log_error("The sequence \"%s\" not found: %s",mplp->regions_itr->seq,baux->fname);
                return -1;
            }
            baux->tid = baux->pos = -1;
        }
        mplp->tid = mplp->pos = -1;
    }
    return -1;  // never happens
}


// keeps cache of the last two chromosomes
static int mplp_cache_refseq(mpileup_t *mplp, int tid, char **ref, int *ref_len)
{
    if  ( !mplp->fai )    // ref not available at all
    {
        *ref = NULL;
        *ref_len = 0;
        return 0;
    }

    if ( tid == mplp->ref[0].tid )  // the one used last time
    {
        *ref = mplp->ref[0].seq;
        *ref_len = mplp->ref[0].len;
        return 1;
    }
    if ( tid == mplp->ref[1].tid )
    {
        refseq_t tmp = mplp->ref[0]; mplp->ref[0] = mplp->ref[1]; mplp->ref[1] = tmp;
        *ref = mplp->ref[0].seq;
        *ref_len = mplp->ref[0].len;
        return 1;
    }

    // New, so migrate to old and load new
    free(mplp->ref[1].seq);
    mplp->ref[1] = mplp->ref[0];
    mplp->ref[0].seq = faidx_fetch_seq(mplp->fai, mplp->chr, 0, INT_MAX, &mplp->ref[0].len);

    if ( !mplp->ref[0].seq )    // sequence is not available
    {
        *ref = NULL;
        *ref_len = 0;
        return 0;
    }
    mplp->ref[0].tid = tid;

    *ref = mplp->ref[0].seq;
    *ref_len = mplp->ref[0].len;
    return 1;
}
static char *mplp_get_ref(mpileup_t *mplp, int *len)
{
    char *ref;
    if ( !mplp_cache_refseq(mplp, mplp->tid, &ref, len) ) return NULL;
    if ( mplp->pos >= *len ) return NULL;
    *len = 1;
    return mplp->ref[0].seq + mplp->pos;
}


// ---------------------------------------
// Legacy mpileup
// ---------------------------------------

static int legacy_mplp_func(void *data, bam1_t *rec)
{
    bam_aux_t *baux = (bam_aux_t*)data;
    mpileup_t *mplp = baux->mplp;
    int ret = 0;
    while (1)
    {
        ret = baux->iter ? sam_itr_next(baux->fp, baux->iter, rec) : sam_read1(baux->fp, baux->hdr, rec);
        if (ret < 0) break;
        if ( rec->core.tid < 0 || (rec->core.flag&BAM_FUNMAP) ) continue;
        if ( mplp->skip_any_unset && (mplp->skip_any_unset & rec->core.flag)!=mplp->skip_any_unset ) continue;
        if ( mplp->skip_all_set && (mplp->skip_all_set & rec->core.flag)==mplp->skip_all_set ) continue;
        if ( mplp->skip_all_unset && !(mplp->skip_all_unset & rec->core.flag) ) continue;
        if ( mplp->skip_any_set && (mplp->skip_any_set & rec->core.flag) ) continue;
// fai .. get_ref
// realign
        if ( rec->core.qual < mplp->min_mq ) continue;
        break;
    }
    return ret;
}
static int legacy_mpileup_init(mpileup_t *mplp)
{
    mplp->legacy.plp  = (bam_pileup1_t**) calloc(mplp->nsmpl,sizeof(*mplp->legacy.plp));
    mplp->legacy.nplp = (int*) calloc(mplp->nsmpl,sizeof(*mplp->legacy.nplp));
    mplp->legacy.iter = bam_mplp_init(mplp->nbam,legacy_mplp_func,(void**)&mplp->bam);
    return 0;
}
static int legacy_mpileup_destroy(mpileup_t *mplp)
{
    if ( !mplp->legacy.active ) return 0;
    free(mplp->legacy.plp);
    free(mplp->legacy.nplp);
    bam_mplp_destroy(mplp->legacy.iter);
    return 0;
}

static int legacy_mpileup_next(mpileup_t *mplp)
{
    while (1)
    {
        int pos;
        int ret = bam_mplp_auto(mplp->legacy.iter, &mplp->tid, &pos, mplp->legacy.nplp, (const bam_pileup1_t**)mplp->legacy.plp);
        mplp->pos = pos;
        if ( ret < 0 ) return ret;
        mplp->chr = (char*)sam_hdr_tid2name(mplp->bam[0].hdr,mplp->tid);
        if ( !ret ) return ret;
        if ( mplp->targets_idx && !regidx_overlap(mplp->targets_idx,mplp->chr,mplp->pos,mplp->pos,0) ) continue;
        if ( mplp->regions_idx && !regidx_overlap(mplp->regions_idx,mplp->chr,mplp->pos,mplp->pos,0) ) continue;
        return ret;
    }
}

