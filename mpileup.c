/*  mpileup.c -- mpileup subcommand. Previously bam_plcmd.c from samtools

    Copyright (C) 2008-2016 Genome Research Ltd.
    Portions copyright (C) 2009-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash_str2int.h>
#include <assert.h>
#include "regidx.h"
#include "bcftools.h"
#include "bam2bcf.h"
#include "sample.h"
#include "gvcf.h"

#define MPLP_BCF        1
#define MPLP_VCF        (1<<1)
#define MPLP_NO_COMP    (1<<2)
#define MPLP_NO_ORPHAN  (1<<3)
#define MPLP_REALN      (1<<4)
#define MPLP_NO_INDEL   (1<<5)
#define MPLP_REDO_BAQ   (1<<6)
#define MPLP_ILLUMINA13 (1<<7)
#define MPLP_IGNORE_RG  (1<<8)
#define MPLP_PRINT_POS  (1<<9)
#define MPLP_PRINT_MAPQ (1<<10)
#define MPLP_PER_SAMPLE (1<<11)
#define MPLP_SMART_OVERLAPS (1<<12)

typedef struct _mplp_aux_t mplp_aux_t;
typedef struct _mplp_pileup_t mplp_pileup_t;

// Data shared by all bam files
typedef struct {
    int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag;
    int rflag_require, rflag_filter, output_type, sample_is_file;
    int openQ, extQ, tandemQ, min_support; // for indels
    double min_frac; // for indels
    char *reg_fname, *pl_list, *fai_fname, *output_fname, *samples_list;
    int reg_is_file, record_cmd_line, n_threads;
    faidx_t *fai;
    regidx_t *bed, *reg;    // bed: skipping regions, reg: index-jump to regions
    int bed_logic;          // 1: include region, 0: exclude region
    gvcf_t *gvcf;

    // auxiliary structures for calling
    bcf_callaux_t *bca;
    bcf_callret1_t *bcr;
    bcf_call_t bc;
    bam_mplp_t iter;
    mplp_aux_t **mplp_data;
    int nfiles;
    char **files;
    mplp_pileup_t *gplp;
    int *n_plp;
    const bam_pileup1_t **plp;
    bam_sample_t *smpl;
    kstring_t buf;
    bcf1_t *bcf_rec;
    htsFile *bcf_fp;
    bcf_hdr_t *bcf_hdr;
    void *rghash_exc;   // shared by all bams

    int argc;
    char **argv;
} mplp_conf_t;

typedef struct {
    char *ref[2];
    int ref_id[2];
    int ref_len[2];
} mplp_ref_t;

#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}

// Data specific to each bam file
struct _mplp_aux_t {
    samFile *fp;
    hts_itr_t *iter;
    bam_hdr_t *h;
    mplp_ref_t *ref;
    const mplp_conf_t *conf;
    void *rghash_inc;   // whitelist for readgroups on per-bam basis, precedes blacklist
    void *rghash_exc;   // blacklist of readgroups 
    hts_idx_t *idx;     // maintained only with more than one -r regions
};

// Data passed to htslib/mpileup
struct _mplp_pileup_t {
    int n;
    int *n_plp, *m_plp;
    bam_pileup1_t **plp;
};

static int mplp_get_ref(mplp_aux_t *ma, int tid,  char **ref, int *ref_len) {
    mplp_ref_t *r = ma->ref;

    //printf("get ref %d {%d/%p, %d/%p}\n", tid, r->ref_id[0], r->ref[0], r->ref_id[1], r->ref[1]);

    if (!r || !ma->conf->fai) {
        *ref = NULL;
        return 0;
    }

    // Do we need to reference count this so multiple mplp_aux_t can
    // track which references are in use?
    // For now we just cache the last two. Sufficient?
    if (tid == r->ref_id[0]) {
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }
    if (tid == r->ref_id[1]) {
        // Last, swap over
        int tmp;
        tmp = r->ref_id[0];  r->ref_id[0]  = r->ref_id[1];  r->ref_id[1]  = tmp;
        tmp = r->ref_len[0]; r->ref_len[0] = r->ref_len[1]; r->ref_len[1] = tmp;

        char *tc;
        tc = r->ref[0]; r->ref[0] = r->ref[1]; r->ref[1] = tc;
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }

    // New, so migrate to old and load new
    free(r->ref[1]);
    r->ref[1]     = r->ref[0];
    r->ref_id[1]  = r->ref_id[0];
    r->ref_len[1] = r->ref_len[0];

    r->ref_id[0] = tid;
    r->ref[0] = faidx_fetch_seq(ma->conf->fai,
                                ma->h->target_name[r->ref_id[0]],
                                0,
                                INT_MAX,
                                &r->ref_len[0]);

    if (!r->ref[0]) {
        r->ref[0] = NULL;
        r->ref_id[0] = -1;
        r->ref_len[0] = 0;
        *ref = NULL;
        return 0;
    }

    *ref = r->ref[0];
    *ref_len = r->ref_len[0];
    return 1;
}

static int mplp_func(void *data, bam1_t *b)
{
    char *ref;
    mplp_aux_t *ma = (mplp_aux_t*)data;
    int ret, ref_len;
    while (1)
    {
        int has_ref;
        ret = ma->iter? sam_itr_next(ma->fp, ma->iter, b) : sam_read1(ma->fp, ma->h, b);
        if (ret < 0) break;
        // The 'B' cigar operation is not part of the specification, considering as obsolete.
        //  bam_remove_B(b);
        if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) continue; // exclude unmapped reads
        if (ma->conf->rflag_require && !(ma->conf->rflag_require&b->core.flag)) continue;
        if (ma->conf->rflag_filter && ma->conf->rflag_filter&b->core.flag) continue;
        if (ma->conf->bed)
        {
            // test overlap
            regitr_t itr;
            int beg = b->core.pos, end = bam_endpos(b)-1;
            int overlap = regidx_overlap(ma->conf->bed, ma->h->target_name[b->core.tid],beg,end, &itr);
            if ( !ma->conf->bed_logic && !overlap )
            {
                // exclude only reads which are fully contained in the region
                while ( REGITR_OVERLAP(itr,beg,end) )
                {
                    if ( beg < REGITR_START(itr) ) { overlap = 1; break; }
                    if ( end > REGITR_END(itr) ) { overlap = 1; break; }
                    itr.i++;
                }
            }
            if ( !overlap ) continue;
        }
        if (ma->rghash_inc)
        {
            // skip unless read group is listed
            uint8_t *rg = bam_aux_get(b, "RG");
            if ( !rg ) continue;
            if ( !khash_str2int_has_key(ma->rghash_inc, (const char*)(rg+1)) ) continue;
        }
        if ( ma->rghash_exc )
        {
            // exclude read groups
            uint8_t *rg = bam_aux_get(b, "RG");
            if ( rg && khash_str2int_has_key(ma->rghash_exc, (const char*)(rg+1)) ) continue;
        }
        if (ma->conf->flag & MPLP_ILLUMINA13) {
            int i;
            uint8_t *qual = bam_get_qual(b);
            for (i = 0; i < b->core.l_qseq; ++i)
                qual[i] = qual[i] > 31? qual[i] - 31 : 0;
        }

        if (ma->conf->fai && b->core.tid >= 0) {
            has_ref = mplp_get_ref(ma, b->core.tid, &ref, &ref_len);
            if (has_ref && ref_len <= b->core.pos) { // exclude reads outside of the reference sequence
                fprintf(stderr,"[%s] Skipping because %d is outside of %d [ref:%d]\n",
                        __func__, b->core.pos, ref_len, b->core.tid);
                continue;
            }
        } else {
            has_ref = 0;
        }

        if (has_ref && (ma->conf->flag&MPLP_REALN)) sam_prob_realn(b, ref, ref_len, (ma->conf->flag & MPLP_REDO_BAQ)? 7 : 3);
        if (has_ref && ma->conf->capQ_thres > 10) {
            int q = sam_cap_mapq(b, ref, ref_len, ma->conf->capQ_thres);
            if (q < 0) continue;    // skip
            else if (b->core.qual > q) b->core.qual = q;
        }
        if (b->core.qual < ma->conf->min_mq) continue;
        else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) continue;

        return ret;
    };
    return ret;
}

static void group_smpl(mplp_pileup_t *m, bam_sample_t *sm, kstring_t *buf,
                       int n, char *const*fn, int *n_plp, const bam_pileup1_t **plp, int ignore_rg)
{
    int i, j;
    memset(m->n_plp, 0, m->n * sizeof(int));
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n_plp[i]; ++j) {
            const bam_pileup1_t *p = plp[i] + j;
            uint8_t *q;
            int id = -1;
            q = ignore_rg? NULL : bam_aux_get(p->b, "RG");
            if (q) id = bam_smpl_rg2smid(sm, fn[i], (char*)q+1, buf);
            if (id < 0) id = bam_smpl_rg2smid(sm, fn[i], 0, buf);
            if (id < 0 || id >= m->n) {
                assert(q); // otherwise a bug
                fprintf(stderr, "[%s] Read group %s used in file %s but absent from the header or an alignment missing read group.\n", __func__, (char*)q+1, fn[i]);
                exit(EXIT_FAILURE);
            }
            if (m->n_plp[id] == m->m_plp[id]) {
                m->m_plp[id] = m->m_plp[id]? m->m_plp[id]<<1 : 8;
                m->plp[id] = (bam_pileup1_t*) realloc(m->plp[id], sizeof(bam_pileup1_t) * m->m_plp[id]);
            }
            m->plp[id][m->n_plp[id]++] = *p;
        }
    }
}

static void flush_bcf_records(mplp_conf_t *conf, htsFile *fp, bcf_hdr_t *hdr, bcf1_t *rec)
{
    if ( !conf->gvcf )
    {
        if ( rec ) bcf_write1(fp, hdr, rec);
        return;
    }

    if ( !rec )
    {
        gvcf_write(conf->gvcf, fp, hdr, NULL, 0);
        return;
    }

    int is_ref = 0;
    if ( rec->n_allele==1 ) is_ref = 1;
    else if ( rec->n_allele==2 )
    {
        // second allele is mpileup's X, not a variant
        if ( rec->d.allele[1][0]=='<' && rec->d.allele[1][1]=='*' && rec->d.allele[1][2]=='>' ) is_ref = 1;
    }
    rec = gvcf_write(conf->gvcf, fp, hdr, rec, is_ref);
    if ( rec ) bcf_write1(fp,hdr,rec);
}

static int mpileup_reg(mplp_conf_t *conf, uint32_t beg, uint32_t end)
{
    bam_hdr_t *hdr = conf->mplp_data[0]->h; // header of first file in input list

    int ret, i, tid, pos, ref_len;
    char *ref;

    while ( (ret=bam_mplp_auto(conf->iter, &tid, &pos, conf->n_plp, conf->plp)) > 0) 
    {
        if ( end && (pos<beg || pos>end) ) continue;
        if ( conf->bed && tid >= 0 )
        {
            int overlap = regidx_overlap(conf->bed, hdr->target_name[tid], pos, pos, NULL);
            if ( !conf->bed_logic ) overlap = overlap ? 0 : 1;
            if ( !overlap ) continue;
        }
        mplp_get_ref(conf->mplp_data[0], tid, &ref, &ref_len);

        int total_depth, _ref0, ref16;
        for (i = total_depth = 0; i < conf->nfiles; ++i) total_depth += conf->n_plp[i];
        group_smpl(conf->gplp, conf->smpl, &conf->buf, conf->nfiles, conf->files, conf->n_plp, conf->plp, conf->flag & MPLP_IGNORE_RG);
        _ref0 = (ref && pos < ref_len)? ref[pos] : 'N';
        ref16 = seq_nt16_table[_ref0];
        bcf_callaux_clean(conf->bca, &conf->bc);
        for (i = 0; i < conf->gplp->n; ++i)
            bcf_call_glfgen(conf->gplp->n_plp[i], conf->gplp->plp[i], ref16, conf->bca, conf->bcr + i);
        conf->bc.tid = tid; conf->bc.pos = pos;
        bcf_call_combine(conf->gplp->n, conf->bcr, conf->bca, ref16, &conf->bc);
        bcf_clear1(conf->bcf_rec);
        bcf_call2bcf(&conf->bc, conf->bcf_rec, conf->bcr, conf->fmt_flag, 0, 0);
        flush_bcf_records(conf, conf->bcf_fp, conf->bcf_hdr, conf->bcf_rec);

        // call indels; todo: subsampling with total_depth>max_indel_depth instead of ignoring?
        // check me: rghash in bcf_call_gap_prep() should have no effect, reads mplp_func already excludes them
        if (!(conf->flag&MPLP_NO_INDEL) && total_depth < conf->max_indel_depth 
            && bcf_call_gap_prep(conf->gplp->n, conf->gplp->n_plp, conf->gplp->plp, pos, conf->bca, ref) >= 0)
        {
            bcf_callaux_clean(conf->bca, &conf->bc);
            for (i = 0; i < conf->gplp->n; ++i)
                bcf_call_glfgen(conf->gplp->n_plp[i], conf->gplp->plp[i], -1, conf->bca, conf->bcr + i);
            if (bcf_call_combine(conf->gplp->n, conf->bcr, conf->bca, -1, &conf->bc) >= 0) 
            {
                bcf_clear1(conf->bcf_rec);
                bcf_call2bcf(&conf->bc, conf->bcf_rec, conf->bcr, conf->fmt_flag, conf->bca, ref);
                flush_bcf_records(conf, conf->bcf_fp, conf->bcf_hdr, conf->bcf_rec);
            }
        }
    }
    return 0;
}

static int mpileup(mplp_conf_t *conf)
{
    extern void *bcf_call_add_rg(void *rghash, const char *hdtext, const char *list);
    extern void bcf_call_del_rghash(void *rghash);

    if (conf->nfiles == 0) {
        fprintf(stderr,"[%s] no input file/data given\n", __func__);
        exit(EXIT_FAILURE);
    }

    mplp_ref_t mp_ref = MPLP_REF_INIT;
    conf->gplp = (mplp_pileup_t *) calloc(1,sizeof(mplp_pileup_t));
    conf->mplp_data = (mplp_aux_t**) calloc(conf->nfiles, sizeof(mplp_aux_t*));
    conf->plp = (const bam_pileup1_t**) calloc(conf->nfiles, sizeof(bam_pileup1_t*));
    conf->n_plp = (int*) calloc(conf->nfiles, sizeof(int));
    conf->smpl = bam_smpl_init();

    int i, samples_logic = 1;   // 1: include, 0: exclude
    void *samples_list = NULL;  // hash of sample names to include or exclude
    if ( conf->samples_list )
    {
        int nsamples = 0;
        char *list = conf->samples_list;
        if ( list[0]=='^' )
        {
            list++;
            samples_logic = 0;
        }
        char **samples = hts_readlist(list, conf->sample_is_file, &nsamples);
        if (!samples) {
            fprintf(stderr,"[%s] could not read the samples list %s\n", __func__, list);
            exit(EXIT_FAILURE);
        }
        if ( nsamples ) samples_list = khash_str2int_init();
        for (i=0; i<nsamples; i++)
        {
            khash_str2int_inc(samples_list, strdup(samples[i]));
            free(samples[i]);
        }
        free(samples);
    }

    // Allow to run mpileup on multiple regions in one go. This comes at cost: the bai index
    // must be kept in the memory for the whole time which can be a problem with many bams.
    // Therefore if none or only one region is requested, we initialize the bam iterator as
    // before and free the index. Only when multiple regions are queried, we keep the index.
    regitr_t reg_itr;
    REGITR_INIT(reg_itr);
    int nregs = 0;
    if ( conf->reg_fname )
    {
        if ( conf->reg_is_file )
        {
            conf->reg = regidx_init(conf->reg_fname,NULL,NULL,0,NULL);
            if ( !conf->reg ) {
                fprintf(stderr,"Could not parse the regions: %s\n", conf->reg_fname);
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            conf->reg = regidx_init(NULL,regidx_parse_reg,NULL,sizeof(char*),NULL);
            if ( regidx_insert_list(conf->reg,conf->reg_fname,',') !=0 ) {
                fprintf(stderr,"Could not parse the regions: %s\n", conf->reg_fname);
                exit(EXIT_FAILURE);
            }
        }
        nregs = regidx_nregs(conf->reg);
        regidx_loop(conf->reg,&reg_itr);    // position the region iterator at the first region
    }

    // read the header of each file in the list and initialize data
    // beware: mpileup has always assumed that tid's are consistent in the headers, add sanity check at least!
    bam_hdr_t *hdr = NULL;      // header of first file in input list
    for (i = 0; i < conf->nfiles; ++i) {
        bam_hdr_t *h_tmp;
        conf->mplp_data[i] = (mplp_aux_t*) calloc(1, sizeof(mplp_aux_t));
        conf->mplp_data[i]->fp = sam_open(conf->files[i], "rb");
        if ( !conf->mplp_data[i]->fp )
        {
            fprintf(stderr, "[%s] failed to open %s: %s\n", __func__, conf->files[i], strerror(errno));
            exit(EXIT_FAILURE);
        }
        if (hts_set_opt(conf->mplp_data[i]->fp, CRAM_OPT_DECODE_MD, 0)) {
            fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
            exit(EXIT_FAILURE);
        }
        if (conf->fai_fname && hts_set_fai_filename(conf->mplp_data[i]->fp, conf->fai_fname) != 0) {
            fprintf(stderr, "[%s] failed to process %s: %s\n",
                    __func__, conf->fai_fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
        conf->mplp_data[i]->conf = conf;
        conf->mplp_data[i]->ref = &mp_ref;
        h_tmp = sam_hdr_read(conf->mplp_data[i]->fp);
        if ( !h_tmp ) {
            fprintf(stderr,"[%s] fail to read the header of %s\n", __func__, conf->files[i]);
            exit(EXIT_FAILURE);
        }
        conf->mplp_data[i]->h = i? hdr : h_tmp; // for i==0, "h" has not been set yet
        conf->mplp_data[i]->rghash_exc = conf->rghash_exc;
        if ( samples_list ) conf->mplp_data[i]->rghash_inc = khash_str2int_init();
        bam_smpl_add(conf->smpl, conf->files[i], (conf->flag&MPLP_IGNORE_RG)? 0 : h_tmp->text, samples_list, samples_logic, conf->mplp_data[i]->rghash_inc);
        // Collect read group IDs with PL (platform) listed in pl_list (note: fragile, strstr search)
        conf->mplp_data[i]->rghash_exc = bcf_call_add_rg(conf->mplp_data[i]->rghash_exc, h_tmp->text, conf->pl_list);
        if (conf->reg) {
            hts_idx_t *idx = sam_index_load(conf->mplp_data[i]->fp, conf->files[i]);
            if (idx == NULL) {
                fprintf(stderr, "[%s] fail to load index for %s\n", __func__, conf->files[i]);
                exit(EXIT_FAILURE);
            }
            conf->buf.l = 0;
            ksprintf(&conf->buf,"%s:%u-%u",regitr_seqname(conf->reg,&reg_itr),REGITR_START(reg_itr)+1,REGITR_END(reg_itr)+1);
            conf->mplp_data[i]->iter = sam_itr_querys(idx, conf->mplp_data[i]->h, conf->buf.s);
            if ( !conf->mplp_data[i]->iter ) 
            {
                conf->mplp_data[i]->iter = sam_itr_querys(idx, conf->mplp_data[i]->h, regitr_seqname(conf->reg,&reg_itr));
                if ( conf->mplp_data[i]->iter ) {
                    fprintf(stderr,"[E::%s] fail to parse region '%s'\n", __func__, conf->buf.s);
                    exit(EXIT_FAILURE);
                }
                fprintf(stderr,"[E::%s] the sequence \"%s\" not found: %s\n",__func__,regitr_seqname(conf->reg,&reg_itr),conf->files[i]);
                exit(EXIT_FAILURE);
            }
            if ( nregs==1 ) // no need to keep the index in memory
               hts_idx_destroy(idx);
            else
                conf->mplp_data[i]->idx = idx;
        }

        if (i == 0) hdr = h_tmp; /* save the header of first file in list */
        else {
            // FIXME: check consistency between h and h_tmp
            bam_hdr_destroy(h_tmp);

            // we store only the first file's header; it's (alleged to be)
            // compatible with the i-th file's target_name lookup needs
            conf->mplp_data[i]->h = hdr;
        }
    }
    // allocate data storage proportionate to number of samples being studied sm->n
    conf->gplp->n = conf->smpl->n;
    conf->gplp->n_plp = (int*) calloc(conf->smpl->n, sizeof(int));
    conf->gplp->m_plp = (int*) calloc(conf->smpl->n, sizeof(int));
    conf->gplp->plp = (bam_pileup1_t**) calloc(conf->smpl->n, sizeof(bam_pileup1_t*));  

    fprintf(stderr, "[%s] %d samples in %d input files\n", __func__, conf->smpl->n, conf->nfiles);
    // write the VCF header
    conf->bcf_fp = hts_open(conf->output_fname?conf->output_fname:"-", hts_bcf_wmode(conf->output_type));
    if (conf->bcf_fp == NULL) {
        fprintf(stderr, "[%s] failed to write to %s: %s\n", __func__, conf->output_fname? conf->output_fname : "standard output", strerror(errno));
        exit(EXIT_FAILURE);
    }
    if ( conf->n_threads ) hts_set_threads(conf->bcf_fp, conf->n_threads);

    // BCF header creation
    conf->bcf_hdr = bcf_hdr_init("w");
    conf->buf.l = 0;

    if (conf->record_cmd_line)
    {
        ksprintf(&conf->buf, "##bcftoolsVersion=%s+htslib-%s\n",bcftools_version(),hts_version());
        bcf_hdr_append(conf->bcf_hdr, conf->buf.s);

        conf->buf.l = 0;
        ksprintf(&conf->buf, "##bcftoolsCommand=mpileup");
        for (i=1; i<conf->argc; i++) ksprintf(&conf->buf, " %s", conf->argv[i]);
        kputc('\n', &conf->buf);
        bcf_hdr_append(conf->bcf_hdr, conf->buf.s);
    }

    if (conf->fai_fname)
    {
        conf->buf.l = 0;
        ksprintf(&conf->buf, "##reference=file://%s\n", conf->fai_fname);
        bcf_hdr_append(conf->bcf_hdr, conf->buf.s);
    }

    // Translate BAM @SQ tags to BCF ##contig tags
    // todo: use/write new BAM header manipulation routines, fill also UR, M5
    for (i=0; i<hdr->n_targets; i++)
    {
        conf->buf.l = 0;
        ksprintf(&conf->buf, "##contig=<ID=%s,length=%d>", hdr->target_name[i], hdr->target_len[i]);
        bcf_hdr_append(conf->bcf_hdr, conf->buf.s);
    }
    conf->buf.l = 0;

    bcf_hdr_append(conf->bcf_hdr,"##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=IDV,Number=1,Type=Integer,Description=\"Maximum number of reads supporting an indel\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=IMF,Number=1,Type=Float,Description=\"Maximum fraction of reads supporting an indel\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)\",Version=\"3\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=RPB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias (bigger is better)\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=MQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias (bigger is better)\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=BQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Base Quality Bias (bigger is better)\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=MQSB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)\">");
#if CDF_MWU_TESTS
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=RPB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias [CDF] (bigger is better)\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=MQB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias [CDF] (bigger is better)\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=BQB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Base Quality Bias [CDF] (bigger is better)\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=MQSB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality vs Strand Bias [CDF] (bigger is better)\">");
#endif
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=SGB,Number=1,Type=Float,Description=\"Segregation based metric.\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Fraction of MQ0 reads (smaller is better)\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">");
    bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">");
    if ( conf->fmt_flag&B2B_FMT_DP )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">");
    if ( conf->fmt_flag&B2B_FMT_DV )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of high-quality non-reference bases\">");
    if ( conf->fmt_flag&B2B_FMT_DPR )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
    if ( conf->fmt_flag&B2B_INFO_DPR )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
    if ( conf->fmt_flag&B2B_FMT_DP4 )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases\">");
    if ( conf->fmt_flag&B2B_FMT_SP )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">");
    if ( conf->fmt_flag&B2B_FMT_AD )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">");
    if ( conf->fmt_flag&B2B_FMT_ADF )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand\">");
    if ( conf->fmt_flag&B2B_FMT_ADR )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand\">");
    if ( conf->fmt_flag&B2B_INFO_AD )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths\">");
    if ( conf->fmt_flag&B2B_INFO_ADF )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total allelic depths on the forward strand\">");
    if ( conf->fmt_flag&B2B_INFO_ADR )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total allelic depths on the reverse strand\">");
    if ( conf->gvcf )
        gvcf_update_header(conf->gvcf, conf->bcf_hdr);

    for (i=0; i<conf->smpl->n; i++)
        bcf_hdr_add_sample(conf->bcf_hdr, conf->smpl->smpl[i]);
    bcf_hdr_add_sample(conf->bcf_hdr, NULL);
    bcf_hdr_write(conf->bcf_fp, conf->bcf_hdr);

    conf->bca = bcf_call_init(-1., conf->min_baseQ);
    conf->bcr = (bcf_callret1_t*) calloc(conf->smpl->n, sizeof(bcf_callret1_t));
    conf->bca->openQ = conf->openQ, conf->bca->extQ = conf->extQ, conf->bca->tandemQ = conf->tandemQ;
    conf->bca->min_frac = conf->min_frac;
    conf->bca->min_support = conf->min_support;
    conf->bca->per_sample_flt = conf->flag & MPLP_PER_SAMPLE;

    conf->bc.bcf_hdr = conf->bcf_hdr;
    conf->bc.n = conf->smpl->n;
    conf->bc.PL = (int32_t*) malloc(15 * conf->smpl->n * sizeof(*conf->bc.PL));
    if (conf->fmt_flag)
    {
        assert( sizeof(float)==sizeof(int32_t) );
        conf->bc.DP4 = (int32_t*) malloc(conf->smpl->n * sizeof(int32_t) * 4);
        conf->bc.fmt_arr = (uint8_t*) malloc(conf->smpl->n * sizeof(float)); // all fmt_flag fields, float and int32
        if ( conf->fmt_flag&(B2B_INFO_DPR|B2B_FMT_DPR|B2B_INFO_AD|B2B_INFO_ADF|B2B_INFO_ADR|B2B_FMT_AD|B2B_FMT_ADF|B2B_FMT_ADR) )
        {
            // first B2B_MAX_ALLELES fields for total numbers, the rest per-sample
            conf->bc.ADR = (int32_t*) malloc((conf->smpl->n+1)*B2B_MAX_ALLELES*sizeof(int32_t));
            conf->bc.ADF = (int32_t*) malloc((conf->smpl->n+1)*B2B_MAX_ALLELES*sizeof(int32_t));
            for (i=0; i<conf->smpl->n; i++)
            {
                conf->bcr[i].ADR = conf->bc.ADR + (i+1)*B2B_MAX_ALLELES;
                conf->bcr[i].ADF = conf->bc.ADF + (i+1)*B2B_MAX_ALLELES;
            }
        }
    }

    // init mpileup
    conf->iter = bam_mplp_init(conf->nfiles, mplp_func, (void**)conf->mplp_data);
    if ( conf->flag & MPLP_SMART_OVERLAPS ) bam_mplp_init_overlaps(conf->iter);
    if ( conf->max_depth * conf->smpl->n > 1<<20)
        fprintf(stderr, "(%s) Max depth is above 1M. Potential memory hog!\n", __func__);
    if ( conf->max_depth * conf->smpl->n < 8000) {
        conf->max_depth = 8000 / conf->smpl->n;
        fprintf(stderr, "<%s> Set max per-file depth to %d\n", __func__, conf->max_depth);
    }
    bam_mplp_set_maxcnt(conf->iter, conf->max_depth);
    conf->max_indel_depth = conf->max_indel_depth * conf->smpl->n;
    conf->bcf_rec = bcf_init1();

    // Run mpileup for multiple regions
    if ( nregs )
    {
        int ireg = 0;
        do 
        {
            // first region is already positioned
            if ( ireg++ > 0 )
            {
                conf->buf.l = 0;
                ksprintf(&conf->buf,"%s:%u-%u",regitr_seqname(conf->reg,&reg_itr),REGITR_START(reg_itr),REGITR_END(reg_itr));

                for (i=0; i<conf->nfiles; i++) 
                {
                    hts_itr_destroy(conf->mplp_data[i]->iter);
                    conf->mplp_data[i]->iter = sam_itr_querys(conf->mplp_data[i]->idx, conf->mplp_data[i]->h, conf->buf.s);
                    if ( !conf->mplp_data[i]->iter ) 
                    {
                        conf->mplp_data[i]->iter = sam_itr_querys(conf->mplp_data[i]->idx, conf->mplp_data[i]->h, regitr_seqname(conf->reg,&reg_itr));
                        if ( conf->mplp_data[i]->iter ) {
                            fprintf(stderr,"[E::%s] fail to parse region '%s'\n", __func__, conf->buf.s);
                            exit(EXIT_FAILURE);
                        }
                        fprintf(stderr,"[E::%s] the sequence \"%s\" not found: %s\n",__func__,regitr_seqname(conf->reg,&reg_itr),conf->files[i]);
                        exit(EXIT_FAILURE);
                    }
                    bam_mplp_reset(conf->iter);
                }
            }
            mpileup_reg(conf,REGITR_START(reg_itr),REGITR_END(reg_itr));
        }
        while ( regidx_loop(conf->reg,&reg_itr) );
    }
    else
        mpileup_reg(conf,0,0);

    flush_bcf_records(conf, conf->bcf_fp, conf->bcf_hdr, NULL);

    // clean up
    free(conf->bc.tmp.s);
    bcf_destroy1(conf->bcf_rec);
    if (conf->bcf_fp)
    {
        hts_close(conf->bcf_fp);
        bcf_hdr_destroy(conf->bcf_hdr);
        bcf_call_destroy(conf->bca);
        free(conf->bc.PL);
        free(conf->bc.DP4);
        free(conf->bc.ADR);
        free(conf->bc.ADF);
        free(conf->bc.fmt_arr);
        free(conf->bcr);
    }
    if ( conf->gvcf ) gvcf_destroy(conf->gvcf);
    bam_smpl_destroy(conf->smpl); free(conf->buf.s);
    for (i = 0; i < conf->gplp->n; ++i) free(conf->gplp->plp[i]);
    free(conf->gplp->plp); free(conf->gplp->n_plp); free(conf->gplp->m_plp); free(conf->gplp);
    bam_mplp_destroy(conf->iter);
    bam_hdr_destroy(hdr);
    for (i = 0; i < conf->nfiles; ++i) {
        if ( nregs>1 ) hts_idx_destroy(conf->mplp_data[i]->idx);
        sam_close(conf->mplp_data[i]->fp);
        if ( conf->mplp_data[i]->iter) hts_itr_destroy(conf->mplp_data[i]->iter);
        if ( conf->mplp_data[i]->rghash_inc ) khash_str2int_destroy_free(conf->mplp_data[i]->rghash_inc);
        free(conf->mplp_data[i]);
    }
    free(conf->mplp_data); free(conf->plp); free(conf->n_plp);
    if ( samples_list ) khash_str2int_destroy_free(samples_list);
    free(mp_ref.ref[0]);
    free(mp_ref.ref[1]);
    return 0;
}

static int is_url(const char *s)
{
    static const char uri_scheme_chars[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+.-";
    return s[strspn(s, uri_scheme_chars)] == ':';
}

#define MAX_PATH_LEN 1024
int read_file_list(const char *file_list,int *n,char **argv[])
{
    char buf[MAX_PATH_LEN];
    int len, nfiles = 0;
    char **files = NULL;
    struct stat sb;

    *n = 0;
    *argv = NULL;

    FILE *fh = fopen(file_list,"r");
    if ( !fh )
    {
        fprintf(stderr,"%s: %s\n", file_list,strerror(errno));
        return 1;
    }

    files = (char**) calloc(nfiles,sizeof(char*));
    nfiles = 0;
    while ( fgets(buf,MAX_PATH_LEN,fh) )
    {
        // allow empty lines and trailing spaces
        len = strlen(buf);
        while ( len>0 && isspace(buf[len-1]) ) len--;
        if ( !len ) continue;

        // check sanity of the file list
        buf[len] = 0;
        if (! (is_url(buf) || stat(buf, &sb) == 0))
        {
            // no such file, check if it is safe to print its name
            int i, safe_to_print = 1;
            for (i=0; i<len; i++)
                if (!isprint(buf[i])) { safe_to_print = 0; break; }
            if ( safe_to_print )
                fprintf(stderr,"The file list \"%s\" appears broken, could not locate: %s\n", file_list,buf);
            else
                fprintf(stderr,"Does the file \"%s\" really contain a list of files and do all exist?\n", file_list);
            return 1;
        }

        nfiles++;
        files = (char**) realloc(files,nfiles*sizeof(char*));
        files[nfiles-1] = strdup(buf);
    }
    fclose(fh);
    if ( !nfiles )
    {
        fprintf(stderr,"No files read from %s\n", file_list);
        return 1;
    }
    *argv = files;
    *n    = nfiles;
    return 0;
}
#undef MAX_PATH_LEN

int parse_format_flag(const char *str)
{
    int i, flag = 0, n_tags;
    char **tags = hts_readlist(str, 0, &n_tags);
    for(i=0; i<n_tags; i++)
    {
        if ( !strcasecmp(tags[i],"DP") || !strcasecmp(tags[i],"FORMAT/DP") || !strcasecmp(tags[i],"FMT/DP") ) flag |= B2B_FMT_DP;
        else if ( !strcasecmp(tags[i],"DV") || !strcasecmp(tags[i],"FORMAT/DV") || !strcasecmp(tags[i],"FMT/DV") ) { flag |= B2B_FMT_DV; fprintf(stderr, "[warning] tag DV functional, but deprecated. Please switch to `AD` in future.\n"); }
        else if ( !strcasecmp(tags[i],"SP") || !strcasecmp(tags[i],"FORMAT/SP") || !strcasecmp(tags[i],"FMT/SP") ) flag |= B2B_FMT_SP;
        else if ( !strcasecmp(tags[i],"DP4") || !strcasecmp(tags[i],"FORMAT/DP4") || !strcasecmp(tags[i],"FMT/DP4") ) { flag |= B2B_FMT_DP4; fprintf(stderr, "[warning] tag DP4 functional, but deprecated. Please switch to `ADF` and `ADR` in future.\n"); }
        else if ( !strcasecmp(tags[i],"DPR") || !strcasecmp(tags[i],"FORMAT/DPR") || !strcasecmp(tags[i],"FMT/DPR") ) { flag |= B2B_FMT_DPR; fprintf(stderr, "[warning] tag DPR functional, but deprecated. Please switch to `AD` in future.\n"); }
        else if ( !strcasecmp(tags[i],"INFO/DPR") ) { flag |= B2B_INFO_DPR; fprintf(stderr, "[warning] tag INFO/DPR functional, but deprecated. Please switch to `INFO/AD` in future.\n"); }
        else if ( !strcasecmp(tags[i],"AD") || !strcasecmp(tags[i],"FORMAT/AD") || !strcasecmp(tags[i],"FMT/AD") ) flag |= B2B_FMT_AD;
        else if ( !strcasecmp(tags[i],"ADF") || !strcasecmp(tags[i],"FORMAT/ADF") || !strcasecmp(tags[i],"FMT/ADF") ) flag |= B2B_FMT_ADF;
        else if ( !strcasecmp(tags[i],"ADR") || !strcasecmp(tags[i],"FORMAT/ADR") || !strcasecmp(tags[i],"FMT/ADR") ) flag |= B2B_FMT_ADR;
        else if ( !strcasecmp(tags[i],"INFO/AD") ) flag |= B2B_INFO_AD;
        else if ( !strcasecmp(tags[i],"INFO/ADF") ) flag |= B2B_INFO_ADF;
        else if ( !strcasecmp(tags[i],"INFO/ADR") ) flag |= B2B_INFO_ADR;
        else
        {
            fprintf(stderr,"Could not parse tag \"%s\" in \"%s\"\n", tags[i], str);
            exit(EXIT_FAILURE);
        }
        free(tags[i]);
    }
    if (n_tags) free(tags);
    return flag;
}

static void list_annotations(FILE *fp)
{
    fprintf(fp,
"\n"
"FORMAT annotation tags available (\"FORMAT/\" prefix is optional):\n"
"\n"
"  FORMAT/AD  .. Allelic depth (Number=R,Type=Integer)\n"
"  FORMAT/ADF .. Allelic depths on the forward strand (Number=R,Type=Integer)\n"
"  FORMAT/ADR .. Allelic depths on the reverse strand (Number=R,Type=Integer)\n"
"  FORMAT/DP  .. Number of high-quality bases (Number=1,Type=Integer)\n"
"  FORMAT/SP  .. Phred-scaled strand bias P-value (Number=1,Type=Integer)\n"
"\n"
"INFO annotation tags available:\n"
"\n"
"  INFO/AD  .. Total allelic depth (Number=R,Type=Integer)\n"
"  INFO/ADF .. Total allelic depths on the forward strand (Number=R,Type=Integer)\n"
"  INFO/ADR .. Total allelic depths on the reverse strand (Number=R,Type=Integer)\n"
"\n");
}

static void print_usage(FILE *fp, const mplp_conf_t *mplp)
{
    char *tmp_require = bam_flag2str(mplp->rflag_require);
    char *tmp_filter  = bam_flag2str(mplp->rflag_filter);

    // Display usage information, formatted for the standard 80 columns.
    // (The unusual string formatting here aids the readability of this
    // source code in 80 columns, to the extent that's possible.)

    fprintf(fp,
"\n"
"Usage: bcftools mpileup [options] in1.bam [in2.bam [...]]\n"
"\n"
"Input options:\n"
"  -6, --illumina1.3+      quality is in the Illumina-1.3+ encoding\n"
"  -A, --count-orphans     do not discard anomalous read pairs\n"
"  -b, --bam-list FILE     list of input BAM filenames, one per line\n"
"  -B, --no-BAQ            disable BAQ (per-Base Alignment Quality)\n"
"  -C, --adjust-MQ INT     adjust mapping quality; recommended:50, disable:0 [0]\n"
"  -d, --max-depth INT     max per-file depth; avoids excessive memory usage [%d]\n", mplp->max_depth);
    fprintf(fp,
"  -E, --redo-BAQ          recalculate BAQ on the fly, ignore existing BQs\n"
"  -f, --fasta-ref FILE    faidx indexed reference sequence file\n"
"  -G, --exclude-RG FILE   exclude read groups listed in FILE\n"
"  -q, --min-MQ INT        skip alignments with mapQ smaller than INT [%d]\n", mplp->min_mq);
    fprintf(fp,
"  -Q, --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [%d]\n", mplp->min_baseQ);
    fprintf(fp,
"  -r, --regions REG[,...] comma separated list of regions in which pileup is generated\n"
"  -R, --regions-file FILE restrict to regions listed in a file\n"
"      --ignore-RG         ignore RG tags (one BAM = one sample)\n"
"  --rf, --incl-flags STR|INT  required flags: skip reads with mask bits unset [%s]\n", tmp_require);
    fprintf(fp,
"  --ff, --excl-flags STR|INT  filter flags: skip reads with mask bits set\n"
"                                            [%s]\n", tmp_filter);
    fprintf(fp,
"  -s, --samples LIST      comma separated list of samples to include\n"
"  -S, --samples-file FILE file of samples to include\n"
"  -t, --targets REG[,...] similar to -r but streams rather than index-jumps\n"
"  -T, --targets-file FILE similar to -R but streams rather than index-jumps\n"
"  -x, --ignore-overlaps   disable read-pair overlap detection\n"
"\n"
"Output options:\n"
"  -a, --annotate LIST     optional tags to output; '?' to list []\n"
"  -g, --gvcf INT[,...]    group non-variant sites into gVCF blocks according\n"
"                          to minimum per-sample DP\n"
"      --no-version        do not append version and command line to the header\n"
"  -o, --output FILE       write output to FILE [standard output]\n"
"  -O, --output-type TYPE  'b' compressed BCF; 'u' uncompressed BCF;\n"
"                          'z' compressed VCF; 'v' uncompressed VCF [v]\n"
"      --threads INT       number of extra output compression threads [0]\n"
"\n"
"SNP/INDEL genotype likelihoods options:\n"
"  -e, --ext-prob INT      Phred-scaled gap extension seq error probability [%d]\n", mplp->extQ);
    fprintf(fp,
"  -F, --gap-frac FLOAT    minimum fraction of gapped reads [%g]\n", mplp->min_frac);
    fprintf(fp,
"  -h, --tandem-qual INT   coefficient for homopolymer errors [%d]\n", mplp->tandemQ);
    fprintf(fp,
"  -I, --skip-indels       do not perform indel calling\n"
"  -L, --max-idepth INT    maximum per-file depth for INDEL calling [%d]\n", mplp->max_indel_depth);
    fprintf(fp,
"  -m, --min-ireads INT    minimum number gapped reads for indel candidates [%d]\n", mplp->min_support);
    fprintf(fp,
"  -o, --open-prob INT     Phred-scaled gap open seq error probability [%d]\n", mplp->openQ);
    fprintf(fp,
"  -p, --per-sample-mF     apply -m and -F per-sample for increased sensitivity\n"
"  -P, --platforms STR     comma separated list of platforms for indels [all]\n"
"\n"
"Notes: Assuming diploid individuals.\n"
"\n");

    free(tmp_require);
    free(tmp_filter);
}

int bam_mpileup(int argc, char *argv[])
{
    int c;
    const char *file_list = NULL;
    char **fn = NULL;
    int nfiles = 0, use_orphan = 0;
    mplp_conf_t mplp;
    memset(&mplp, 0, sizeof(mplp_conf_t));
    mplp.min_baseQ = 13;
    mplp.capQ_thres = 0;
    mplp.max_depth = 250; mplp.max_indel_depth = 250;
    mplp.openQ = 40; mplp.extQ = 20; mplp.tandemQ = 100;
    mplp.min_frac = 0.002; mplp.min_support = 1;
    mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_SMART_OVERLAPS;
    mplp.argc = argc; mplp.argv = argv;
    mplp.rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    mplp.output_fname = NULL;
    mplp.output_type = FT_VCF;
    mplp.record_cmd_line = 1;
    mplp.n_threads = 0;

    static const struct option lopts[] =
    {
        {"rf", required_argument, NULL, 1},   // require flag
        {"ff", required_argument, NULL, 2},   // filter flag
        {"incl-flags", required_argument, NULL, 1},
        {"excl-flags", required_argument, NULL, 2},
        {"output", required_argument, NULL, 3},
        {"open-prob", required_argument, NULL, 4},
        {"ignore-RG", no_argument, NULL, 5},
        {"ignore-rg", no_argument, NULL, 5},
        {"gvcf", required_argument, NULL, 'g'},
        {"no-version", no_argument, NULL, 8},
        {"threads",required_argument,NULL,9},
        {"illumina1.3+", no_argument, NULL, '6'},
        {"count-orphans", no_argument, NULL, 'A'},
        {"bam-list", required_argument, NULL, 'b'},
        {"no-BAQ", no_argument, NULL, 'B'},
        {"no-baq", no_argument, NULL, 'B'},
        {"adjust-MQ", required_argument, NULL, 'C'},
        {"adjust-mq", required_argument, NULL, 'C'},
        {"max-depth", required_argument, NULL, 'd'},
        {"redo-BAQ", no_argument, NULL, 'E'},
        {"redo-baq", no_argument, NULL, 'E'},
        {"fasta-ref", required_argument, NULL, 'f'},
        {"exclude-RG", required_argument, NULL, 'G'},
        {"exclude-rg", required_argument, NULL, 'G'},
        {"region", required_argument, NULL, 'r'},
        {"regions", required_argument, NULL, 'r'},
        {"regions-file", required_argument, NULL, 'R'},
        {"targets", required_argument, NULL, 't'},
        {"targets-file", required_argument, NULL, 'T'},
        {"min-MQ", required_argument, NULL, 'q'},
        {"min-mq", required_argument, NULL, 'q'},
        {"min-BQ", required_argument, NULL, 'Q'},
        {"min-bq", required_argument, NULL, 'Q'},
        {"ignore-overlaps", no_argument, NULL, 'x'},
        {"output-type", required_argument, NULL, 'O'},
        {"samples", required_argument, NULL, 's'},
        {"samples-file", required_argument, NULL, 'S'},
        {"annotate", required_argument, NULL, 'a'},
        {"ext-prob", required_argument, NULL, 'e'},
        {"gap-frac", required_argument, NULL, 'F'},
        {"tandem-qual", required_argument, NULL, 'h'},
        {"skip-indels", no_argument, NULL, 'I'},
        {"max-idepth", required_argument, NULL, 'L'},
        {"min-ireads ", required_argument, NULL, 'm'},
        {"per-sample-mF", no_argument, NULL, 'p'},
        {"per-sample-mf", no_argument, NULL, 'p'},
        {"platforms", required_argument, NULL, 'P'},
        {NULL, 0, NULL, 0}
    };
    while ((c = getopt_long(argc, argv, "Ag:f:r:R:q:Q:C:Bd:L:b:P:po:e:h:Im:F:EG:6O:xa:s:S:t:T:",lopts,NULL)) >= 0) {
        switch (c) {
        case 'x': mplp.flag &= ~MPLP_SMART_OVERLAPS; break;
        case  1 :
            mplp.rflag_require = bam_str2flag(optarg);
            if ( mplp.rflag_require<0 ) { fprintf(stderr,"Could not parse --rf %s\n", optarg); return 1; }
            break;
        case  2 :
            mplp.rflag_filter = bam_str2flag(optarg);
            if ( mplp.rflag_filter<0 ) { fprintf(stderr,"Could not parse --ff %s\n", optarg); return 1; }
            break;
        case  3 : mplp.output_fname = optarg; break;
        case  4 : mplp.openQ = atoi(optarg); break;
        case  5 : mplp.flag |= MPLP_IGNORE_RG; break;
        case 'g':
            mplp.gvcf = gvcf_init(optarg);
            if ( !mplp.gvcf ) error("Could not parse: --gvcf %s\n", optarg);
            break;
        case 'f':
            mplp.fai = fai_load(optarg);
            if (mplp.fai == NULL) return 1;
            mplp.fai_fname = optarg;
            break;
        case  8 : mplp.record_cmd_line = 0; break;
        case  9 : mplp.n_threads = strtol(optarg, 0, 0); break;
        case 'd': mplp.max_depth = atoi(optarg); break;
        case 'r': mplp.reg_fname = strdup(optarg); break;
        case 'R': mplp.reg_fname = strdup(optarg); mplp.reg_is_file = 1; break;
        case 't':
                  // In the original version the whole BAM was streamed which is inefficient
                  //  with few BED intervals and big BAMs. Todo: devise a heuristic to determine
                  //  best strategy, that is streaming or jumping.
                  if ( optarg[0]=='^' ) optarg++;
                  else mplp.bed_logic = 1;
                  mplp.bed = regidx_init(NULL,regidx_parse_reg,NULL,0,NULL);
                  if ( regidx_insert_list(mplp.bed,optarg,',') !=0 )
                  {
                      fprintf(stderr,"Could not parse the targets: %s\n", optarg);
                      exit(EXIT_FAILURE);
                  }
                  break;
        case 'T':
                  if ( optarg[0]=='^' ) optarg++;
                  else mplp.bed_logic = 1;
                  mplp.bed = regidx_init(optarg,NULL,NULL,0,NULL);
                  if (!mplp.bed) { fprintf(stderr, "bcftools mpileup: Could not read file \"%s\"", optarg); return 1; }
                  break;
        case 'P': mplp.pl_list = strdup(optarg); break;
        case 'p': mplp.flag |= MPLP_PER_SAMPLE; break;
        case 'B': mplp.flag &= ~MPLP_REALN; break;
        case 'I': mplp.flag |= MPLP_NO_INDEL; break;
        case 'E': mplp.flag |= MPLP_REDO_BAQ; break;
        case '6': mplp.flag |= MPLP_ILLUMINA13; break;
        case 's': mplp.samples_list = optarg; break;
        case 'S':
            mplp.samples_list = optarg;
            mplp.sample_is_file = 1;
            break;
        case 'O': 
            switch (optarg[0]) {
                case 'b': mplp.output_type = FT_BCF_GZ; break;
                case 'u': mplp.output_type = FT_BCF; break;
                case 'z': mplp.output_type = FT_VCF_GZ; break;
                case 'v': mplp.output_type = FT_VCF; break;
                default: error("[error] The option \"-O\" changed meaning when mpileup moved to bcftools. Did you mean: \"bcftools mpileup --output-type\" or \"samtools mpileup --output-BP\"?\n", optarg); 
            }
            break;
        case 'C': mplp.capQ_thres = atoi(optarg); break;
        case 'q': mplp.min_mq = atoi(optarg); break;
        case 'Q': mplp.min_baseQ = atoi(optarg); break;
        case 'b': file_list = optarg; break;
        case 'o': {
                char *end;
                long value = strtol(optarg, &end, 10);
                // Distinguish between -o INT and -o FILE (a bit of a hack!)
                if (*end == '\0') mplp.openQ = value;
                else mplp.output_fname = optarg;
            }
            break;
        case 'e': mplp.extQ = atoi(optarg); break;
        case 'h': mplp.tandemQ = atoi(optarg); break;
        case 'A': use_orphan = 1; break;
        case 'F': mplp.min_frac = atof(optarg); break;
        case 'm': mplp.min_support = atoi(optarg); break;
        case 'L': mplp.max_indel_depth = atoi(optarg); break;
        case 'G': {
                FILE *fp_rg;
                char buf[1024];
                mplp.rghash_exc = khash_str2int_init();
                if ((fp_rg = fopen(optarg, "r")) == NULL)
                    fprintf(stderr, "(%s) Fail to open file %s. Continue anyway.\n", __func__, optarg);
                while (!feof(fp_rg) && fscanf(fp_rg, "%s", buf) > 0) // this is not a good style, but forgive me...
                    khash_str2int_inc(mplp.rghash_exc, strdup(buf));
                fclose(fp_rg);
            }
            break;
        case 'a':
            if (optarg[0]=='?') {
                list_annotations(stderr);
                return 1;
            }
            mplp.fmt_flag |= parse_format_flag(optarg);
        break;
        default:
            fprintf(stderr,"Invalid option: '%c'\n", c);
            return 1;
        }
    }

    if ( mplp.gvcf && !(mplp.fmt_flag&B2B_FMT_DP) )
    {
        fprintf(stderr,"[warning] The -t DP option is required with --gvcf, switching on.\n");
        mplp.fmt_flag |= B2B_FMT_DP;
    }
    if ( mplp.flag&(MPLP_BCF|MPLP_VCF|MPLP_NO_COMP) )
    {
        if ( mplp.flag&MPLP_VCF )
        {
            if ( mplp.flag&MPLP_NO_COMP ) mplp.output_type = FT_VCF;
            else mplp.output_type = FT_VCF_GZ;
        }
        else if ( mplp.flag&MPLP_BCF )
        {
            if ( mplp.flag&MPLP_NO_COMP ) mplp.output_type = FT_BCF;
            else mplp.output_type = FT_BCF_GZ;
        }
    }
    if ( !(mplp.flag&MPLP_REALN) && mplp.flag&MPLP_REDO_BAQ )
    {
        fprintf(stderr,"Error: The -B option cannot be combined with -E\n");
        return 1;
    }
    if (use_orphan) mplp.flag &= ~MPLP_NO_ORPHAN;
    if (argc == 1)
    {
        print_usage(stderr, &mplp);
        return 1;
    }
    int ret;
    if (file_list) 
    {
        if ( read_file_list(file_list,&nfiles,&fn) ) return 1;
        mplp.files  = fn;
        mplp.nfiles = nfiles;
        ret = mpileup(&mplp);
        for (c=0; c<nfiles; c++) free(fn[c]);
        free(fn);
    }
    else
    {
        mplp.nfiles = argc - optind;
        mplp.files  = argv + optind;
        ret = mpileup(&mplp);
    }
    if (mplp.rghash_exc) khash_str2int_destroy_free(mplp.rghash_exc);
    free(mplp.reg_fname); free(mplp.pl_list);
    if (mplp.fai) fai_destroy(mplp.fai);
    if (mplp.bed) regidx_destroy(mplp.bed);
    if (mplp.reg) regidx_destroy(mplp.reg);
    return ret;
}
