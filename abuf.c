/* The MIT License

   Copyright (c) 2021-2026 Genome Research Ltd.

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
#include <limits.h>
#include <strings.h>
#include <htslib/vcf.h>
#include <ctype.h>
#include "bcftools.h"
#include "abuf.h"
#include "rbuf.h"

typedef enum
{
    M_FIRST, M_SUM
}
merge_rule_t;

typedef struct
{
    kstring_t ref, alt;
    int ial;        // the index of the original ALT allele, 1-based
    int beg, end;   // 0-based inclusive offsets to ref,alt
    int plen;       // the ref,alt prefix length, eg plen=1 for C>CA
}
atom_t;

typedef struct
{
    bcf1_t *rec;
    int nori, nout;     // number of ALTs in the input, and VCF rows on output
    uint8_t *tbl;       // nori columns, nout rows; indicates allele contribution to output rows, see "The atomization works as follows" below
    uint8_t *overlaps;  // is the star allele needed for this variant?
    atom_t **atoms;
    int matoms, mtbl, moverlaps;
    char *info_tag;
}
split_t;

typedef struct
{
    int beg, end, len;       // half-open REF interval and replacement length
    const char *alt;
}
join_edit_t;

typedef struct
{
    bcf1_t *rec;             // buffered input record
    hts_pos_t ref_len;       // cached literal REF length
    int32_t *gt;             // working genotypes, loaded only when joining records with GT
    int mgt, width;          // allocation and per-sample width, including padding
    int duplicate_of;       // earlier matching record, or -1
}
join_rec_t;

typedef struct
{
    join_rec_t *records; // records of the currently open block, in position order, not yet merged
    int nrec, mrec; // number of records in the open block; allocated size of records
    hts_pos_t end;  // 0-based, exclusive end of the block's REF span accumulated so far
    int snvs_only;
    uint64_t njoined; // number of joined output records, excluding singletons

    // Reusable working buffers. Active lengths are reset for each merged block.
    int *smpl_ploidy;
    uint8_t *is_phased;
    int msmpl_ploidy, mis_phased;
    int32_t *alt_indices, *joined_gt;
    int malt_indices, mjoined_gt, ngt;
    int *copy_map;   // shared allocation for the unphased allocator's map, seen and records arrays
    int mcopy_map;   // capacity in integers, not chromosome copies
    join_edit_t *edits;
    uint8_t *covered;
    int medits, mcovered;
    kstring_t ref, seq, allele_str;
    kstring_t *alt;
    int nalt, malt;
}
join_t;

struct _abuf_t
{
    abuf_opt_t mode;
    split_t split;
    join_t join;
    atom_t *atoms;
    int natoms, matoms;
    const bcf_hdr_t *hdr;
    bcf_hdr_t *out_hdr;
    bcf1_t **vcf;       // dimensions stored in rbuf
    rbuf_t rbuf;

    kstring_t tmps;
    void *tmp, *tmp2;
    int32_t *gt, *tmpi;
    int ngt, mgt, ntmpi, mtmpi, mtmp, mtmp2;
    int star_allele;
};

abuf_t *abuf_init(const bcf_hdr_t *hdr, abuf_opt_t mode)
{
    if ( mode!=SPLIT && mode!=JOIN ) error("todo\n");
    abuf_t *buf = (abuf_t*) calloc(1,sizeof(abuf_t));
    buf->hdr  = hdr;
    buf->out_hdr = (bcf_hdr_t*) hdr;
    buf->mode = mode;
    buf->star_allele = 1;
    rbuf_init(&buf->rbuf, 0);
    return buf;
}

void abuf_destroy(abuf_t *buf)
{
    int i;
    for (i=0; i<buf->matoms; i++)
    {
        free(buf->atoms[i].ref.s);
        free(buf->atoms[i].alt.s);
    }
    free(buf->atoms);
    free(buf->split.atoms);
    free(buf->split.overlaps);
    free(buf->split.tbl);
    for (i=0; i<buf->rbuf.m; i++)
        if ( buf->vcf[i] ) bcf_destroy(buf->vcf[i]);
    free(buf->vcf);
    for (i=0; i<buf->join.mrec; i++)
    {
        if ( buf->join.records[i].rec ) bcf_destroy(buf->join.records[i].rec);
        free(buf->join.records[i].gt);
    }
    free(buf->join.records);
    free(buf->join.smpl_ploidy);
    free(buf->join.is_phased);
    free(buf->join.alt_indices);
    free(buf->join.joined_gt);
    free(buf->join.copy_map);
    free(buf->join.edits);
    free(buf->join.covered);
    free(buf->join.ref.s);
    free(buf->join.seq.s);
    free(buf->join.allele_str.s);
    for (i=0; i<buf->join.malt; i++) free(buf->join.alt[i].s);
    free(buf->join.alt);
    free(buf->gt);
    free(buf->tmpi);
    free(buf->tmp);
    free(buf->tmp2);
    free(buf->tmps.s);
    free(buf);
}

void abuf_set(abuf_t *buf, abuf_opt_t key, void *value)
{
    if ( key==BCF_HDR ) { buf->out_hdr = *((bcf_hdr_t**)value); return; }
    if ( key==INFO_TAG )
    {
        buf->split.info_tag = *((char**)value);
        bcf_hdr_printf(buf->out_hdr,"##INFO=<ID=%s,Number=1,Type=String,Description=\"Original variant. Format: CHR|POS|REF|ALT|USED_ALT_IDX\">",buf->split.info_tag);
        return;
    }
    if ( key==STAR_ALLELE ) { buf->star_allele = *((int*)value); return; }
    if ( key==SNVS_ONLY ) { buf->join.snvs_only = *((int*)value); return; }
}

/*
    Split alleles into primitives, e.g.
        CC>TT  becomes  C>T,C>T
        GCGT>GTGA  becomes C>T,T>A

    There is no sequence alignment, just trimming and hungry matching
    from left side.
*/
static void _atomize_allele(abuf_t *buf, bcf1_t *rec, int ial)
{
    // Trim identical sequence from right
    char *ref = rec->d.allele[0];
    char *alt = rec->d.allele[ial];
    int rlen = strlen(ref);
    int alen = strlen(alt);
    while ( rlen>1 && alen>1 && ref[rlen-1]==alt[alen-1] ) rlen--, alen--;
    int Mlen = rlen > alen ? rlen : alen;

    atom_t *atom = NULL;
    int i;
    for (i=0; i<Mlen; i++)
    {
        char refb = i<rlen ? ref[i] : '-';
        char altb = i<alen ? alt[i] : '-';
        if ( refb!=altb )
        {
            if ( refb=='-' || altb=='-' )
            {
                assert(atom);
                if ( altb!='-' ) kputc(altb, &atom->alt);
                if ( refb!='-' ) { kputc(refb, &atom->ref); atom->end++; }
                continue;
            }
            buf->natoms++;
            hts_expand0(atom_t,buf->natoms,buf->matoms,buf->atoms);
            atom = &buf->atoms[buf->natoms-1];
            atom->ref.l = 0;
            atom->alt.l = 0;
            kputc(refb, &atom->ref);
            kputc(altb, &atom->alt);
            atom->beg = atom->end = i;
            atom->ial = ial;
            atom->plen = 0;

            if ( rlen!=alen && (i+1>=rlen || i+1>=alen) )   // the next base is an indel combined with SNV, e.g. C>GGG?
            {
                buf->natoms++;
                hts_expand0(atom_t,buf->natoms,buf->matoms,buf->atoms);
                atom = &buf->atoms[buf->natoms-1];
                atom->ref.l = 0;
                atom->alt.l = 0;
                kputc(refb, &atom->ref);
                kputc(refb, &atom->alt);
                atom->beg  = atom->end = i;
                atom->ial  = ial;
                atom->plen = 1;
            }
            continue;
        }
        if ( i+1>=rlen || i+1>=alen )   // is the next base an indel?
        {
            buf->natoms++;
            hts_expand0(atom_t,buf->natoms,buf->matoms,buf->atoms);
            atom = &buf->atoms[buf->natoms-1];
            atom->ref.l = 0;
            atom->alt.l = 0;
            kputc(refb, &atom->ref);
            kputc(altb, &atom->alt);
            atom->beg = atom->end = i;
            atom->ial = ial;
            atom->plen = 0;
        }
    }
}
static int _atoms_inconsistent(const atom_t *a, const atom_t *b)
{
    if ( a->beg < b->beg ) return -1;
    if ( a->beg > b->beg ) return 1;
    int rcmp = strcasecmp(a->ref.s,b->ref.s);
    if ( rcmp ) return rcmp;
    return strcasecmp(a->alt.s,b->alt.s);
}

// returns
//      0 .. identical beg,ref,alt
//      1 .. non-overlapping variants, but record may overlap (A>AT vs A>C)
//      2 .. overlapping (conflicting) variants
static int _atoms_overlap(const atom_t *a, const atom_t *b)
{
    if ( a->beg < b->beg ) return 2;
    if ( a->beg > b->beg ) return 2;

    // consider SNV followed by DEL as not overlapping
    //      CC > C      a.plen=1 (ref,alt prefix len=1)
    //      C  > T      b.plen=0 (ref,alt prefix len=0)
    if ( a->plen && a->plen >= b->ref.l ) return 1;
    if ( b->plen && b->plen >= a->ref.l ) return 1;

    int rcmp = strcasecmp(a->ref.s,b->ref.s);
    if ( rcmp ) return 2;

    // consider SNV followed by INS as not overlapping
    //      A > AT      a.plen=1 (ref,alt prefix len=1)
    //      A > C       b.plen=0 (ref,alt prefix len=0)
    if ( a->plen && a->plen >= b->alt.l ) return 1;
    if ( b->plen && b->plen >= a->alt.l ) return 1;

    rcmp = strcasecmp(a->alt.s,b->alt.s);
    if ( rcmp ) return 2;
    return 0;
}
/*
    For reproducibility of tests on different platforms, we need to guarantee the same order of identical
    atoms originating from different source ALTs.  Even though they are consistent, different values can be
    picked for VCF annotations as currently the values from the one that comes first are used.
*/
static int _cmp_atoms(const void *aptr, const void *bptr)
{
    const atom_t *a = (const atom_t*) aptr;
    const atom_t *b = (const atom_t*) bptr;
    int rcmp = _atoms_inconsistent(a,b);
    if ( rcmp ) return rcmp;
    if ( a->ial < b->ial ) return -1;
    if ( a->ial > b->ial ) return 1;
    return 0;
}
static void _split_table_init(abuf_t *buf, bcf1_t *rec, int natoms)
{
    buf->split.rec  = rec;
    buf->split.nori = rec->n_allele - 1;
    buf->split.nout = 0;
    hts_expand(uint8_t,buf->split.nori*natoms,buf->split.mtbl,buf->split.tbl);
    hts_expand(atom_t*,natoms,buf->split.matoms,buf->split.atoms);
    hts_expand(uint8_t,natoms,buf->split.moverlaps,buf->split.overlaps);
    memset(buf->split.overlaps,0,sizeof(*buf->split.overlaps)*natoms);
}
static void _split_table_new(abuf_t *buf, atom_t *atom)
{
    int i, iout = buf->split.nout++;
    buf->split.atoms[iout] = atom;
    uint8_t *ptr = buf->split.tbl + iout*buf->split.nori;
    for (i=0; i<buf->split.nori; i++) ptr[i] = 0;
    ptr[atom->ial-1] = 1;
}
static void _split_table_overlap(abuf_t *buf, int iout, atom_t *atom)
{
    uint8_t *ptr = buf->split.tbl + iout*buf->split.nori;
    int olap = _atoms_overlap(atom,buf->split.atoms[iout]);
    ptr[atom->ial-1] =  olap > 1 ? 2 : 1;

    // The test test/atomize.split.5.vcf shows why we sometimes can and sometimes
    // cannot remove the star allele like this
    //      buf->split.overlaps[iout] = olap > 1 ? 1 : 0;
    // I forgot the details of the code, so don't immediately see
    // if this could be made smarter
    buf->split.overlaps[iout] = 1;
}
#if 0
static void _split_table_print(abuf_t *buf)
{
    int i,j;
    for (i=0; i<buf->split.nout; i++)
    {
        atom_t *atom = buf->split.atoms[i];
        uint8_t *ptr = buf->split.tbl + i*buf->split.nori;
        fprintf(stderr,"%d\t%s\t%s",(int)buf->split.rec->pos+1+atom->beg,atom->ref.s,atom->alt.s);
        for (j=0; j<buf->split.nori; j++) fprintf(stderr,"\t%d",(int)ptr[j]);
        fprintf(stderr,"\n");
    }
}
static void _split_table_print_atoms(abuf_t *buf)
{
    int i;
    for (i=0; i<buf->natoms; i++)
    {
        atom_t *atom = &buf->atoms[i];
        fprintf(stderr,"atom%d %p: ialt=%d %s>%s %d-%d\n",i,atom,atom->ial,atom->ref.s,atom->alt.s,atom->beg,atom->end);
    }
}
#endif
static inline uint8_t _has_star_allele(abuf_t *buf, int iout)
{
    if ( !buf->star_allele ) return 0;
    return buf->split.overlaps[iout];
}
static inline int _split_table_get_ial(abuf_t *buf, int irow, int ial)
{
    if ( !ial ) return ial;
    return buf->split.tbl[irow*buf->split.nori + ial - 1];
}
static void _split_table_set_chrom_qual(abuf_t *buf)
{
    int iout,j;
    bcf1_t *rec = buf->split.rec;
    for (iout=0; iout<buf->split.nout; iout++)
    {
        rbuf_expand0(&buf->rbuf, bcf1_t*, buf->rbuf.n+1, buf->vcf);
        j = rbuf_append(&buf->rbuf);
        if ( !buf->vcf[j] ) buf->vcf[j] = bcf_init1();
        bcf1_t *out = buf->vcf[j];
        bcf_clear1(out);

        atom_t *atom = buf->split.atoms[iout];
        out->rid = rec->rid;
        out->pos = rec->pos + atom->beg;
        bcf_update_id(buf->out_hdr, out, rec->d.id);

        const char *als[3];
        als[0] = atom->ref.s;
        als[1] = atom->alt.s;
        als[2] = "*";
        int nals = _has_star_allele(buf,iout) ? 3 : 2;
        bcf_update_alleles(buf->out_hdr, out, als, nals);

        if ( bcf_float_is_missing(rec->qual) )
            bcf_float_set_missing(out->qual);
        else
            out->qual = rec->qual;

        bcf_update_filter(buf->out_hdr, out, rec->d.flt, rec->d.n_flt);
    }
}
int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst);
static void _split_table_set_info(abuf_t *buf, bcf_info_t *info, merge_rule_t mode)
{
    const char *tag = bcf_hdr_int2id(buf->hdr,BCF_DT_ID,info->key);
    int type = bcf_hdr_id2type(buf->hdr,BCF_HL_INFO,info->key);
    int len  = bcf_hdr_id2length(buf->hdr,BCF_HL_INFO,info->key);
    if ( len==BCF_VL_G ) return;                                                // todo: Number=G INFO tags
    if ( type==BCF_HT_LONG ) return;                                            // todo: 64bit integers

    bcf1_t *rec = buf->split.rec;
    int mtmp = ( type==BCF_HT_INT || type==BCF_HT_REAL ) ? buf->mtmp/4 : buf->mtmp;
    int nval = bcf_get_info_values(buf->hdr,rec,tag,&buf->tmp,&mtmp,type);
    if ( type==BCF_HT_INT || type==BCF_HT_REAL ) buf->mtmp = mtmp*4;

    // Check for incorrect number of values. Note this check does not consider all values missing
    // and will remove annotations that don't pass.
    if ( type==BCF_HT_INT || type==BCF_HT_REAL )
    {
        if ( (len==BCF_VL_A && nval != rec->n_allele - 1) || (len==BCF_VL_R && nval != rec->n_allele) ) return;
    }

    if ( buf->mtmp2 < buf->mtmp )
    {
        buf->tmp2  = realloc(buf->tmp2, buf->mtmp);
        if ( !buf->tmp2 ) error("Failed to alloc %d bytes\n", buf->mtmp);
        buf->mtmp2 = buf->mtmp;
    }

    const int num_size = 4;
    assert( num_size==sizeof(int32_t) && num_size==sizeof(float) );
    int32_t missing = bcf_int32_missing;
    void *missing_ptr = (void*)&missing;
    if ( type==BCF_HT_REAL ) bcf_float_set_missing(*((float*)missing_ptr));
    int32_t vector_end = bcf_int32_vector_end;
    void *vector_end_ptr = (void*)&vector_end;
    if ( type==BCF_HT_REAL ) bcf_float_set_vector_end(*((float*)vector_end_ptr));

    int iout,i;
    for (iout=0; iout<buf->split.nout; iout++)
    {
        bcf1_t *out = buf->vcf[rbuf_kth(&buf->rbuf,iout)];
        int star_allele = _has_star_allele(buf,iout);
        int ret = 0;
        if ( len==BCF_VL_FIXED || len==BCF_VL_VAR )
            ret = bcf_update_info(buf->out_hdr, out, tag, type==BCF_HT_FLAG ? NULL : buf->tmp, nval, type);
        else if ( len==BCF_VL_A && type!=BCF_HT_STR )
        {
            int iori = buf->split.atoms[iout]->ial - 1;
            assert( iori<nval );
            if ( !memcmp(vector_end_ptr,buf->tmp+num_size*iori,num_size) )
                memcpy(buf->tmp2,missing_ptr,num_size);
            else
                memcpy(buf->tmp2,buf->tmp+num_size*iori,num_size);
            if ( star_allele )
                memcpy(buf->tmp2+num_size,missing_ptr,num_size);
            ret = bcf_update_info(buf->out_hdr, out, tag, buf->tmp2, 1 + star_allele, type);
        }
        else if ( len==BCF_VL_A && type==BCF_HT_STR )
        {
            int iori = buf->split.atoms[iout]->ial - 1;
            kstring_t dst;
            dst.l = 0; dst.m = buf->mtmp2; dst.s = (char*)buf->tmp2;
            kputc('.',&dst);
            if ( star_allele ) kputs(",.",&dst);
            copy_string_field(buf->tmp, iori, nval, &dst, 0);
            if ( star_allele ) copy_string_field(".", 0, 1, &dst, 1);
            buf->mtmp2 = dst.m;
            buf->tmp2  = dst.s;
            ret = bcf_update_info(buf->out_hdr, out, tag, buf->tmp2, dst.l, type);
        }
        else if ( len==BCF_VL_R && type!=BCF_HT_STR )
        {
            memcpy(buf->tmp2,buf->tmp,num_size);   // REF contributes to all records
            int iori = buf->split.atoms[iout]->ial;
            assert( iori<nval && iori<=buf->split.nori );
            if ( !memcmp(vector_end_ptr,buf->tmp+num_size*iori,num_size) )
                memcpy(buf->tmp2+num_size,missing_ptr,num_size);
            else
                memcpy(buf->tmp2+num_size,buf->tmp+num_size*iori,num_size);
            if ( type==BCF_HT_INT && mode==M_SUM )
            {
                uint8_t *tbl = buf->split.tbl + iout*buf->split.nori;
                for (i=iori; i<buf->split.nori; i++)
                {
                    if ( tbl[i]==1 ) ((int32_t*)buf->tmp2)[1] += ((int32_t*)buf->tmp)[i+1];
                }
            }
            if ( star_allele )
                memcpy(buf->tmp2+2*num_size,missing_ptr,num_size);
            ret = bcf_update_info(buf->out_hdr, out, tag, buf->tmp2, 2 + star_allele, type);
        }
        else if ( len==BCF_VL_R && type==BCF_HT_STR )
        {
            int iori = buf->split.atoms[iout]->ial - 1;
            kstring_t dst;
            dst.l = 0; dst.m = buf->mtmp2; dst.s = (char*)buf->tmp2;
            kputs(".,.",&dst);
            if ( star_allele ) kputs(",.",&dst);
            copy_string_field(buf->tmp, 0, nval, &dst, 0);
            copy_string_field(buf->tmp, iori+1, nval, &dst, 1);
            if ( star_allele ) copy_string_field(".", 0, 1, &dst, 2);
            buf->mtmp2 = dst.m;
            buf->tmp2  = dst.s;
            ret = bcf_update_info(buf->out_hdr, out, tag, buf->tmp2, dst.l, type);
        }
        if ( ret!=0 ) error("An error occurred while updating INFO/%s (errcode=%d)\n",tag,ret);
    }
}
static void _split_table_set_history(abuf_t *buf)
{
    int i,j,ret;
    bcf1_t *rec = buf->split.rec;

    // Don't update if the tag already exists. This is to prevent -a from overwriting -m
    int m = 0;
    char *tmp = NULL;
    ret = bcf_get_info_string(buf->hdr,rec,buf->split.info_tag,&tmp,&m);
    free(tmp);
    if ( ret>0 ) return;

    buf->tmps.l = 0;
    ksprintf(&buf->tmps,"%s|%"PRIhts_pos"|%s|",bcf_seqname(buf->hdr,rec),rec->pos+1,rec->d.allele[0]);
    for (i=1; i<rec->n_allele; i++)
    {
        kputs(rec->d.allele[i],&buf->tmps);
        if ( i+1<rec->n_allele ) kputc(',',&buf->tmps);
        else kputc(',',&buf->tmps);
    }
    int len = buf->tmps.l;
    buf->tmps.s[buf->tmps.l-1] = '|';

    for (i=0; i<buf->split.nout; i++)
    {
        buf->tmps.l = len;
        bcf1_t *out = buf->vcf[rbuf_kth(&buf->rbuf,i)];
        uint8_t *ptr = buf->split.tbl + i*buf->split.nori;
        for (j=0; j<buf->split.nori; j++)
        {
            if ( ptr[j]!=1 ) continue;
            kputw(j+1,&buf->tmps);
            kputc(',',&buf->tmps);
        }
        buf->tmps.s[--buf->tmps.l] = 0;
        if ( (ret=bcf_update_info_string(buf->out_hdr, out, buf->split.info_tag, buf->tmps.s))!=0 )
            error("An error occurred while updating INFO/%s (errcode=%d)\n",buf->split.info_tag,ret);
    }
}
static void _split_table_set_gt(abuf_t *buf)
{
    int nsmpl = bcf_hdr_nsamples(buf->hdr);
    if ( !nsmpl ) return;

    bcf1_t *rec = buf->split.rec;
    buf->ngt = bcf_get_genotypes(buf->hdr, rec, &buf->gt, &buf->mgt);
    if ( buf->ngt<=0 ) return;
    else
        hts_expand(int32_t,buf->ngt,buf->mtmpi,buf->tmpi);

    int iout,i,j;
    for (iout=0; iout<buf->split.nout; iout++)
    {
        bcf1_t *out = buf->vcf[rbuf_kth(&buf->rbuf,iout)];
        int star_allele = _has_star_allele(buf,iout);
        int max_ploidy = buf->ngt/nsmpl;
        int32_t *src = buf->gt, *dst = buf->tmpi;
        for (i=0; i<nsmpl; i++)
        {
            for (j=0; j<max_ploidy; j++)
            {
                if ( src[j]==bcf_int32_vector_end || bcf_gt_is_missing(src[j]) )
                {
                    dst[j] = src[j];
                    continue;
                }
                int iori = bcf_gt_allele(src[j]);
                if ( iori<0 || iori>=rec->n_allele )
                    error("Out-of-bounds genotypes at %s:%"PRIhts_pos"\n",bcf_seqname(buf->hdr,rec),rec->pos+1);
                int ial = _split_table_get_ial(buf,iout,iori);
                if ( ial==2 && !star_allele )
                {
                    dst[j] = bcf_gt_missing;
                    if ( bcf_gt_is_phased(src[j]) ) dst[j] |= 1;
                }
                else
                    dst[j] = bcf_gt_is_phased(src[j]) ? bcf_gt_phased(ial) : bcf_gt_unphased(ial);
            }
            src += max_ploidy;
            dst += max_ploidy;
        }
        bcf_update_genotypes(buf->out_hdr,out,buf->tmpi,buf->ngt);
    }
}
static void _split_table_set_format(abuf_t *buf, bcf_fmt_t *fmt, merge_rule_t mode)
{
    int nsmpl = bcf_hdr_nsamples(buf->hdr);
    if ( !nsmpl ) return;

    const char *tag = bcf_hdr_int2id(buf->hdr,BCF_DT_ID,fmt->id);
    if ( tag[0]=='G' && tag[1]=='T' && !tag[2] )        // FORMAT/GT
    {
        _split_table_set_gt(buf);
        return;
    }

    int type = bcf_hdr_id2type(buf->hdr,BCF_HL_FMT,fmt->id);
    int len  = bcf_hdr_id2length(buf->hdr,BCF_HL_FMT,fmt->id);
    if ( type==BCF_HT_STR && len==BCF_VL_G ) return;                            // possible todo: Number=G for strings
    if ( type==BCF_HT_LONG ) return;                                            // todo: 64bit integers

    const int num_size = 4;
    assert( num_size==sizeof(int32_t) && num_size==sizeof(float) );
    int32_t missing = bcf_int32_missing;
    void *missing_ptr = (void*)&missing;
    if ( type==BCF_HT_REAL ) bcf_float_set_missing(*((float*)missing_ptr));
    int32_t vector_end = bcf_int32_vector_end;
    void *vector_end_ptr = (void*)&vector_end;
    if ( type==BCF_HT_REAL ) bcf_float_set_vector_end(*((float*)vector_end_ptr));

    bcf1_t *rec = buf->split.rec;
    int mtmp = ( type==BCF_HT_INT || type==BCF_HT_REAL ) ? buf->mtmp/num_size : buf->mtmp;  // number of items
    int nval = bcf_get_format_values(buf->hdr,rec,tag,&buf->tmp,&mtmp,type);
    if ( type==BCF_HT_INT || type==BCF_HT_REAL ) buf->mtmp = mtmp*num_size;                 // number of bytes

    if ( type==BCF_HT_INT || type==BCF_HT_REAL )
    {
        if ( len==BCF_VL_G && nval!=nsmpl*rec->n_allele && nval!=nsmpl*rec->n_allele*(rec->n_allele+1)/2 ) return;      // not haploid nor diploid

        // Check for incorrect number of values. Note this check does not consider all values missing
        // and will remove annotations that don't pass.
        if ( (len==BCF_VL_A && nval != nsmpl*(rec->n_allele - 1)) || (len==BCF_VL_R && nval != nsmpl*rec->n_allele) ) return;
    }

    // Increase buffer size to accommodate star allele
    int nval1 = nval / nsmpl;
    mtmp = buf->mtmp;
    if ( type==BCF_HT_INT || type==BCF_HT_REAL )
    {
        if ( (len==BCF_VL_A || len==BCF_VL_R) && mtmp < num_size*nsmpl*(nval1+1) ) mtmp = num_size*nsmpl*(nval1+1); // +1 for the possibility of the star allele
        else if ( len==BCF_VL_G && mtmp < num_size*nsmpl*(nval1+3) ) mtmp = num_size*nsmpl*(nval1+3);
    }
    else if ( type==BCF_HT_STR )
    {
        if ( (len==BCF_VL_A || len==BCF_VL_R) && mtmp < nsmpl*(nval1+2) ) mtmp = nsmpl*(nval1+2); // +2 for the possibility of the star allele, ",."
        else if ( len==BCF_VL_G && mtmp < nsmpl*(nval1+6) ) mtmp = nsmpl*(nval1+6);
    }

    if ( buf->mtmp2 < mtmp )
    {
        buf->tmp2  = realloc(buf->tmp2, mtmp);
        if ( !buf->tmp2 ) error("Failed to alloc %d bytes\n", mtmp);
        buf->mtmp2 = mtmp;
    }

    int iout, i, j;
    for (iout=0; iout<buf->split.nout; iout++)
    {
        int star_allele = _has_star_allele(buf,iout);
        bcf1_t *out = buf->vcf[rbuf_kth(&buf->rbuf,iout)];
        int ret = 0;
        if ( len==BCF_VL_FIXED || len==BCF_VL_VAR )
            ret = bcf_update_format(buf->out_hdr, out, tag, buf->tmp, nval, type);
        else if ( len==BCF_VL_A && type!=BCF_HT_STR )
        {
            int iori = buf->split.atoms[iout]->ial - 1;
            assert( iori<nval );
            for (i=0; i<nsmpl; i++)
            {
                void *src = buf->tmp  + nval1*num_size*i;
                void *dst = buf->tmp2 + num_size*i*(star_allele+1);
                if ( !memcmp(vector_end_ptr,src+iori*num_size,num_size) )
                    memcpy(dst,missing_ptr,num_size);
                else
                    memcpy(dst,src+iori*num_size,num_size);
                if ( star_allele )
                    memcpy(dst+num_size,missing_ptr,num_size);
            }
            ret = bcf_update_format(buf->out_hdr, out, tag, buf->tmp2, nsmpl*(star_allele+1), type);
        }
        else if ( (len==BCF_VL_A || len==BCF_VL_R) && type==BCF_HT_STR )
        {
            int ioff = len==BCF_VL_R ? 1 : 0;
            int iori = buf->split.atoms[iout]->ial - 1;
            int nval1_dst = star_allele ? nval1 + 2 : nval1;
            memset(buf->tmp2,0,nval1_dst*nsmpl);
            for (i=0; i<nsmpl; i++)
            {
                char *src = (char*)buf->tmp2 + nval1_dst*i;
                kstring_t dst;
                dst.l = 0; dst.m = nval1_dst + 1; dst.s = src;
                kputc_('.',&dst);
                if ( star_allele ) kputsn_(",.",2,&dst);
                if ( len==BCF_VL_R )
                {
                    kputsn_(",.",2,&dst);
                    copy_string_field(buf->tmp+nval1*i, 0, nval1, &dst, 0);
                }
                copy_string_field(buf->tmp+nval1*i, iori+ioff, nval1, &dst, 0+ioff);
                if ( star_allele ) copy_string_field(".", 0, 1, &dst, 1+ioff);
                assert( src==dst.s );   // make sure copy_string_field() did not reallocate
            }
            ret = bcf_update_format(buf->out_hdr, out, tag, buf->tmp2, nval1_dst*nsmpl, type);
        }
        else if ( len==BCF_VL_R && type!=BCF_HT_STR )
        {
            int iori = buf->split.atoms[iout]->ial;
            assert( iori<=nval );
            for (i=0; i<nsmpl; i++)
            {
                void *src = buf->tmp  + nval1*num_size*i;
                void *dst = buf->tmp2 + num_size*i*(star_allele+2);
                memcpy(dst,src,num_size);
                memcpy(dst+num_size,src+iori*num_size,num_size);
                if ( type==BCF_HT_INT && mode==M_SUM )
                {
                    uint8_t *tbl = buf->split.tbl + iout*buf->split.nori;
                    for (j=iori; j<buf->split.nori; j++)
                        if ( tbl[j]==1 ) ((int32_t*)dst)[1] += ((int32_t*)src)[j+1];
                }
                if ( star_allele )
                    memcpy(dst+num_size*2,missing_ptr,num_size);
            }
            ret = bcf_update_format(buf->out_hdr, out, tag, buf->tmp2, nsmpl*(star_allele+2), type);
        }
        else if ( len==BCF_VL_G && type!=BCF_HT_STR )
        {
            int iori = buf->split.atoms[iout]->ial;
            int i01  = bcf_alleles2gt(0,iori);
            int i11  = bcf_alleles2gt(iori,iori);
            assert( iori<nval );
            #define BRANCH(type_t, is_missing, is_vector_end, set_missing, set_vector_end) { \
                for (i=0; i<nsmpl; i++) \
                { \
                    type_t *src = (type_t*)buf->tmp + i*nval1; \
                    type_t *dst = (type_t*)buf->tmp2 + i*3*(1+star_allele); \
                    int n=0; /* determine ploidy of this genotype */ \
                    while ( n<nval1 && !(is_vector_end) ) { n++; src++; } \
                    src = (type_t*)buf->tmp + i*nval1; \
                    memcpy(dst++,src,sizeof(type)); \
                    int nmiss = 0, nend = 0; \
                    if ( n==rec->n_allele ) /* haploid */ \
                    { \
                        memcpy(dst++,src+iori,sizeof(type)); \
                        if ( star_allele ) { nmiss = 1; nend = 3; } \
                        else nend = 1; \
                    } \
                    else if ( n==nval1 ) \
                    { \
                        memcpy(dst++,src+i01,sizeof(type)); \
                        memcpy(dst++,src+i11,sizeof(type)); \
                        if ( star_allele ) nmiss = 3; \
                    } \
                    else if ( n==1 && is_missing ) \
                    { \
                        if ( star_allele ) nend = 5; \
                        else nend = 2; \
                    } \
                    else  \
                        error("Incorrect number of values at %s:%"PRIhts_pos" .. tag=FORMAT/%s Number=G nAlleles=%d nValues=%d, %d-th sample\n", \
                                bcf_seqname(buf->hdr,rec),rec->pos+1,tag,rec->n_allele,n,i+1); \
                    for (j=0; j<nmiss; j++) { set_missing; dst++; } \
                    for (j=0; j<nend; j++) { set_vector_end; dst++; } \
                } \
            }
            switch (type)
            {
                case BCF_HT_INT:  BRANCH(int32_t, *src==bcf_int32_missing, *src==bcf_int32_vector_end, *dst=bcf_int32_missing, *dst=bcf_int32_vector_end); break;
                case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(*src), bcf_float_is_vector_end(*src), bcf_float_set_missing(*dst), bcf_float_set_vector_end(*dst)); break;
                default: error("Unexpected case: %d\n", type);
            }
            #undef BRANCH
            ret = bcf_update_format(buf->out_hdr, out, tag, buf->tmp2, 3*(1+star_allele)*nsmpl, type);
        }
        if ( ret!=0 ) error("An error occurred while updating FORMAT/%s (errcode=%d)\n",tag,ret);
    }
}
static inline int _is_acgtn(char *seq)
{
    while ( *seq )
    {
        char c = toupper_c(*seq);
        if ( c!='A' && c!='C' && c!='G' && c!='T' && c!='N' ) return 0;
        seq++;
    }
    return 1;
}
/*
    The atomization works as follows:
    - Atomize each alternate allele separately by leaving out sequence identical to the reference. No
      alignment is performed, just greedy trimming of the end, then from left. This operation returns
      a list of atoms (atom_t) which carry fragments of REF,ALT and their positions as 0-based offsets
      to the original REF allele
    - Sort atoms by POS, REF and ALT. Each unique atom (POS+REF+ALT) forms a new VCF record, each
      with a single ALT.
    - For each new VCF record determine how to translate the original allele index (iori) to this new
      record:
        - 1: the original allele matches the atom
        - 0: the original allele does not overlap this atom or the overlapping part matches the REF
             allele
        - 2 (or equivalently "."): there is a mismatch between the original allele and the atom
      The mapping is encoded in a table with columns corresponding to the original ALTs and rows
      to the new POS+ALTs (atoms). The table is initialized to 0, then we set 1's for matching
      atoms and 2's for overlapping mismatching atoms.

    Note that different ALT alleles can result in the same atom (the same output line) and this code
    does not know how to reconcile possibly conflicting VCF annotations. This could be improved
    and merge logic provided, similarly to `merge -l`. For example, the allelic depths (AD) should
    be summed for the same atomized output allele. However, this level of complexity is not addressed
    in this initial draft. Higher priority for now is to provide the inverse "join" operation.

    Update 2021-04-09:
        Tags QS,AD are now automatically incremented as they should be, for both INFO and FORMAT.
        Note that the code will fail on missing values (todo) and it needs to be generalized and
        made customizable.
*/
void _abuf_split(abuf_t *buf, bcf1_t *rec)
{
    int i,j;
    if ( rec->n_allele < 2 )
    {
        rbuf_expand0(&buf->rbuf, bcf1_t*, buf->rbuf.n+1, buf->vcf);
        int j = rbuf_append(&buf->rbuf);
        if ( buf->vcf[j] ) bcf_destroy(buf->vcf[j]);
        buf->vcf[j] = bcf_dup(rec);
        return;
    }
    for (i=0; i<rec->n_allele; i++)
    {
        if ( _is_acgtn(rec->d.allele[i]) ) continue;
        rbuf_expand0(&buf->rbuf, bcf1_t*, buf->rbuf.n+1, buf->vcf);
        int j = rbuf_append(&buf->rbuf);
        if ( buf->vcf[j] ) bcf_destroy(buf->vcf[j]);
        buf->vcf[j] = bcf_dup(rec);
        return;
    }

    buf->natoms = 0;
    for (i=1; i<rec->n_allele; i++) _atomize_allele(buf,rec,i);
    qsort(buf->atoms,buf->natoms,sizeof(*buf->atoms),_cmp_atoms);
    _split_table_init(buf,rec,buf->natoms);
    for (i=0; i<buf->natoms; i++)
    {
        if ( i && _atoms_inconsistent(&buf->atoms[i-1],&buf->atoms[i])==0 ) continue;
        _split_table_new(buf, &buf->atoms[i]);  // add a new unique output atom
    }
    for (i=0; i<buf->natoms; i++)
    {
        // Looping over sorted list of all atoms with possible duplicates from different source ALT alleles
        atom_t *atom = &buf->atoms[i];
        for (j=0; j<buf->split.nout; j++)
        {
            atom_t *out = buf->split.atoms[j];
            if ( atom == out ) continue;            // table already set to 1
            if ( atom->beg > out->end ) continue;   // cannot overlap this output atom
            if ( atom->end < out->beg ) break;      // this atom is ahead of all subsequent output records
            _split_table_overlap(buf, j, atom);
        }
    }
    // _split_table_print(buf);
    // _split_table_print_atoms(buf);
    assert( !buf->rbuf.n ); // all records should be flushed first in the SPLIT mode

    // Create the output records, transferring all annotations:
    // CHROM-QUAL
    _split_table_set_chrom_qual(buf);

    // INFO
    for (i=0; i<rec->n_info; i++)
    {
        // this implementation of merging rules is temporary: generalize and made customizable through the API
        merge_rule_t mode = M_FIRST;
        const char *tag = bcf_hdr_int2id(buf->hdr,BCF_DT_ID,rec->d.info[i].key);
        if ( !strcmp(tag,"QS") || !strcmp(tag,"AD") ) mode = M_SUM;

        _split_table_set_info(buf, &rec->d.info[i], mode);
    }

    // Set INFO tag showing the original record
    if ( buf->split.info_tag )
        _split_table_set_history(buf);

    // FORMAT
    for (i=0; i<rec->n_fmt; i++)
    {
        // this implementation of merging rules is temporary: generalize and made customizable through the API
        merge_rule_t mode = M_FIRST;
        const char *tag = bcf_hdr_int2id(buf->hdr,BCF_DT_ID,rec->d.fmt[i].id);
        if ( !strcmp(tag,"QS") || !strcmp(tag,"AD") ) mode = M_SUM;

        _split_table_set_format(buf, &rec->d.fmt[i], mode);
    }

    // Check that at least one FORMAT field was added, if not, the number of samples must be set manually
    for (i=0; i<buf->split.nout; i++)
    {
        bcf1_t *out = buf->vcf[rbuf_kth(&buf->rbuf,i)];
        if ( !out->n_sample ) out->n_sample = rec->n_sample;
    }
}

/*
    JOIN groups adjacent or overlapping REF spans. It is not a general inverse of
    SPLIT: annotations cannot be reconstructed, and unphased haplotypes are not
    uniquely determined. Only GT is merged; other FORMAT fields are dropped.

    For unphased samples, non-reference alleles are assigned to the rightmost free
    copies within each position group, e.g. 1/0 -> 0/1. This deterministic heuristic
    recovers split compound heterozygotes, but can also combine variants that were
    in trans.
    Incompatible edits are rejected rather than overwriting an earlier allele.
*/
static int _abuf_join_edit_cmp(const void *ap, const void *bp)
{
    const join_edit_t *a = ap, *b = bp;
    if ( a->beg!=b->beg ) return a->beg < b->beg ? -1 : 1;
    if ( a->end!=b->end ) return a->end < b->end ? -1 : 1;
    return 0;
}

// trim sequence identical in REF and ALT, from both ends
static int _abuf_join_edit(const char *ref, int ref_end, const char *alt, int ref_offset, join_edit_t *edit)
{
    int alt_end = strlen(alt), beg = 0;
    while ( beg<ref_end && beg<alt_end && toupper_c(ref[beg])==toupper_c(alt[beg]) ) beg++;
    while ( ref_end>beg && alt_end>beg && toupper_c(ref[ref_end-1])==toupper_c(alt[alt_end-1]) ) { ref_end--; alt_end--; }
    *edit = (join_edit_t){ref_offset+beg,ref_offset+ref_end,alt_end-beg,alt+beg};
    return beg!=ref_end || beg!=alt_end;
}

// Compare two edits. Identical edits and substitutions that agree on their
// overlapping bases are compatible; insertions conflict at the same boundary.
static int _abuf_join_edits_conflict(const join_edit_t *edit, const join_edit_t *query)
{
    if ( edit->beg==query->beg && edit->end==query->end && edit->len==query->len
            && !strncasecmp(edit->alt,query->alt,edit->len) ) return 0;
    if ( edit->len==edit->end-edit->beg && query->len==query->end-query->beg )
    {
        int beg = edit->beg > query->beg ? edit->beg : query->beg;
        int end = edit->end < query->end ? edit->end : query->end;
        return beg<end && strncasecmp(edit->alt+beg-edit->beg,query->alt+beg-query->beg,end-beg)!=0;
    }
    if ( edit->beg<query->end && query->beg<edit->end ) return 1;
    return edit->beg==edit->end && query->beg==query->end && edit->beg==query->beg;
}

// Compare edit against a sorted list of edits.
static int _abuf_join_conflicts(const join_edit_t *edits, int nedit, const join_edit_t *query)
{
    int lo = 0, hi = nedit;
    while ( lo<hi )
    {
        int mid = lo+(hi-lo)/2;
        if ( edits[mid].end<query->beg ) lo = mid+1;
        else hi = mid;
    }
    for (; lo<nedit && edits[lo].beg<=query->end; lo++)
        if ( _abuf_join_edits_conflict(&edits[lo],query) ) return 1;
    return 0;
}


// find the first earlier record with the same POS, REF and ordered ALT list
static void _abuf_join_find_duplicates(join_t *blk)
{
    int irec, grp_beg = 0;
    for (irec=0; irec<blk->nrec; irec++)
    {
        bcf1_t *rec = blk->records[irec].rec;
        if ( rec->pos!=blk->records[grp_beg].rec->pos ) grp_beg = irec;
        blk->records[irec].duplicate_of = -1;
        int iprev;
        for (iprev=grp_beg; iprev<irec; iprev++)
        {
            bcf1_t *prev = blk->records[iprev].rec;
            if ( prev->n_allele!=rec->n_allele ) continue;
            int ial;
            for (ial=0; ial<rec->n_allele; ial++)
                if ( strcasecmp(prev->d.allele[ial],rec->d.allele[ial]) ) break;
            if ( ial<rec->n_allele ) continue;
            blk->records[irec].duplicate_of = iprev;
            break;
        }
    }
}

// Return sorted ALT indices and count missing copies for one sample's GT.
// Keep multiplicity: 1/1 produces [1,1], since both chromosome copies need an ALT.
// Sorting makes unphased 2/1 and 1/2 identical for assignment and duplicate checks.
static int _abuf_join_gt_alts(int32_t *gt, int ploidy, int32_t *alts, int *nmissing)
{
    int icopy, nalt = 0;
    *nmissing = 0;
    for (icopy=0; icopy<ploidy; icopy++)
    {
        int alt = bcf_gt_allele(gt[icopy]);     // alt=-1 for missing, 0 for ref
        if ( alt<0 ) (*nmissing)++;
        else if ( alt>0 )
        {
            int ial = nalt++;
            while ( ial && alts[ial-1]>alt )
                alts[ial] = alts[ial-1], ial--;
            alts[ial] = alt;
        }
    }
    return nalt;
}

// Compare edit against earlier ALTs assigned in this sample's chr copy
//      prev_beg  .. start of earlier records potentially overlapping the current position
//      grp_beg   .. start of records at the current position, i.e. earlier relevant records are [prev_beg,grp_beg)
//      query     .. the candidate edit, prepared once before trying destinations
//      dst       .. chr copy, i.e. {0,1} for diploid GTs
static int _abuf_join_copy_conflicts(join_t *blk, int ismpl, int prev_beg, int grp_beg, const join_edit_t *query, int dst)
{
    int iprev;
    for (iprev=prev_beg; iprev<grp_beg; iprev++)
    {
        int32_t call = blk->records[iprev].gt[ismpl*blk->records[iprev].width+dst];
        if ( bcf_gt_is_missing(call) || !bcf_gt_allele(call) ) continue;
        bcf1_t *prev = blk->records[iprev].rec;
        join_edit_t edit;
        if ( !_abuf_join_edit(prev->d.allele[0],blk->records[iprev].ref_len,prev->d.allele[bcf_gt_allele(call)],prev->pos-blk->records[0].rec->pos,&edit) ) continue;
        if ( _abuf_join_edits_conflict(&edit,query) ) return 1;
    }
    return 0;
}

// Assign ALT occurrence src to a compatible chromosome copy. Consider this example, with two records
// at the same POS:
//      record 0: REF=A ALT=C GT=0/1
//      record 1: REF=A ALT=G GT=0/1
//
// We have two ALTs to place, later records come first (the order is to place first records to rightmost positions):
//      src   records[src]   alt_indices[src]    Meaning
//      0     1              1                   record 1's first ALT: G
//      1     0              1                   record 0's first ALT: C
//
// The array map[0..ploidy) identifies where we placed the ALTs. Initially:
//      map = [-1,-1]
// First, place first ALT (G, src=0) on the preferred rightmost copy:
//      map = [-1, 0]   . G
// Next, try to place second ALT (C, src=1) on the rightmost copy. But as it is occupied, so first move G to 0:
//      map = [ 0, 1]   G C
//
// The array seen prevents infinite recursion:
//      seen = [0, 0]   start searching for C
//             [0, 1]   visit copy 1, occupied by G
//             [1, 1]   visit copy 0, relocate G there
//
// Mappings:
//      src = map[dst]          // which alt occupies the copy dst
//      irec = records[src]     // index of its owning record in blk->records
//      rec = blk->records[irec].rec // buffered input record
//      ial = alt_indices[src]  // which allele in that record
//
// Returns 1 after placing src, or 0 if no placement exists. Assignments at earlier positions are not modified.
//
static int _abuf_join_match_copy(join_t *blk, int ismpl, int prev_beg, int grp_beg, int *records, int src, int *map, int *seen)
{
    join_edit_t query;
    int has_edit = 0;
    if ( prev_beg<grp_beg )
    {
        bcf1_t *rec = blk->records[records[src]].rec;
        has_edit = _abuf_join_edit(rec->d.allele[0],blk->records[records[src]].ref_len,rec->d.allele[blk->alt_indices[src]],rec->pos-blk->records[0].rec->pos,&query);
    }
    int dst;
    for (dst=blk->smpl_ploidy[ismpl]-1; dst>=0; dst--)
    {
        if ( seen[dst] ) continue;
        if ( has_edit && _abuf_join_copy_conflicts(blk,ismpl,prev_beg,grp_beg,&query,dst) ) continue;

        // Visit each chr copy at most once. First try the rightmost compatible slot, even if occupied. If occupied,
        // try recursively relocate its ALT, possibly through a chain of moves ending at an empty slot.
        // Only update map after that succeeds, a failed search makes no changes
        seen[dst] = 1;
        if ( map[dst]<0 || _abuf_join_match_copy(blk,ismpl,prev_beg,grp_beg,records,map[dst],map,seen) )
        {
            map[dst] = src;
            return 1;
        }
    }
    return 0;
}

// Unphased genotypes are currently handled as follows:
//  - take ALT occurrences from all records at this POS, each needs a separate
//    chr copy. This recovers split multiallelics from 0/1 to 1/2. Duplicate
//    records are recognized and don't increase dosage
//  - match those occurrences to earlier compatible copies, favor rightmost
//    placement in cis
//  - after ALTs are placed, fill missing and ref alleles
//
static void _abuf_join_assign_unphased(abuf_t *buf, int ismpl)
{
    join_t *blk = &buf->join;
    join_rec_t *rows = blk->records;
    int smpl_ploidy = blk->smpl_ploidy[ismpl];
    int32_t *alt_indices = blk->alt_indices;

    // One reusable allocation split into three arrays, each of length P=smpl_ploidy:
    //      map[dst]     = ALT-occurrence index assigned to chr copy dst, or -1 when unassigned
    //      seen[dst]    = whether the chr copy dst was visited in the current matching search
    //      records[src] = index in blk->records of the input record owning ALT occurrence src
    //      alt_indices[src] = the src's allele index in the input record (records[src])
    //
    if ( smpl_ploidy>INT_MAX/3 ) error("Error: too many genotype copies to join\n");
    hts_expand(int,3*smpl_ploidy,blk->mcopy_map,blk->copy_map);
    int *map = blk->copy_map, *seen = map+smpl_ploidy, *records = seen+smpl_ploidy;
    int grp_end = 0, prev_beg = 0;
    while ( grp_end<blk->nrec )
    {
        // [grp_beg,grp_end) contains all records at this POS. Their ALT occurrences compete for distinct
        // copies, regardless of their REF lengths
        int grp_beg = grp_end, irec, icopy, nalt = 0;
        while ( grp_end<blk->nrec && blk->records[grp_end].rec->pos==blk->records[grp_beg].rec->pos ) grp_end++;

        // Build a list of all ALT occurrences across all records at this POS.
        for (irec=grp_end-1; irec>=grp_beg; irec--)
        {
            bcf1_t *rec = blk->records[irec].rec;
            int32_t *smpl_gt = rows[irec].gt + ismpl*rows[irec].width;

            // The arrays seen/map are reused as scratch space
            int *alts = seen, *dups = map;   // reuse the arrays temporarily, they're not needed yet
            int nmissing, nrec_alt = _abuf_join_gt_alts(smpl_gt,smpl_ploidy,alts,&nmissing);
            int idup = rows[irec].duplicate_of;
            if ( idup>=0 )
            {
                // Neither GT has been rewritten yet. Check allele multiplicities and missing counts, but
                // ignore their original unphased order.
                int prev_nmissing;
                int32_t *prev_gt = rows[idup].gt + ismpl*rows[idup].width;
                int prev_nalt = _abuf_join_gt_alts(prev_gt,smpl_ploidy,dups,&prev_nmissing);
                if ( prev_nalt!=nrec_alt || prev_nmissing!=nmissing || memcmp(dups,alts,nrec_alt*sizeof(*alts)) )
                    error("Error: inconsistent duplicate genotypes for sample %s at %s:%"PRIhts_pos"\n", buf->hdr->samples[ismpl],bcf_seqname(buf->hdr,rec),rec->pos+1);
                continue;
            }
            if ( nrec_alt>smpl_ploidy-nalt )
                error("Error: too many unphased ALT copies for sample %s at %s:%"PRIhts_pos"\n", buf->hdr->samples[ismpl],bcf_seqname(buf->hdr,rec),rec->pos+1);

            // Save the group's ALT occurrences
            for (icopy=0; icopy<nrec_alt; icopy++)
            {
                alt_indices[nalt] = alts[icopy];
                records[nalt++] = irec;
            }
        }

        // find where to start checking ealier records that could span the current pos
        while ( prev_beg<grp_beg && blk->records[prev_beg].rec->pos+blk->records[prev_beg].ref_len<blk->records[grp_beg].rec->pos ) prev_beg++;

        // match occurrences to chr copies
        for (icopy=0; icopy<smpl_ploidy; icopy++) map[icopy] = -1;
        for (icopy=0; icopy<nalt; icopy++)
        {
            // in _abuf_join_match_copy the previous placement can be also relocated, e.g. "./C" + G -> "C/G"
            memset(seen,0,smpl_ploidy*sizeof(*seen));
            if ( !_abuf_join_match_copy(blk,ismpl,prev_beg,grp_beg,records,icopy,map,seen) )
                error("Error: cannot assign unphased ALT copies to compatible haplotypes for sample %s at %s:%"PRIhts_pos"\n",
                      buf->hdr->samples[ismpl],bcf_seqname(buf->hdr,blk->records[grp_beg].rec),blk->records[grp_beg].rec->pos+1);
        }

        // matching is complete, use mapped ALT occurrences (map[]) to construct final GTs
        for (irec=grp_beg; irec<grp_end; irec++)
        {
            int32_t *smpl_gt = rows[irec].gt + ismpl*rows[irec].width;
            int iprev = rows[irec].duplicate_of;
            if ( iprev>=0 )
            {
                // duplicate records come after the first occurrence which has GT already set: simply copy them
                memcpy(smpl_gt,rows[iprev].gt+ismpl*rows[iprev].width,smpl_ploidy*sizeof(*smpl_gt));
                continue;
            }
            int nmissing = 0;
            for (icopy=0; icopy<smpl_ploidy; icopy++)
                if ( bcf_gt_is_missing(smpl_gt[icopy]) ) nmissing++;
            for (icopy=0; icopy<smpl_ploidy; icopy++)
            {
                // first fills ALT occurrences, followed by missing and REF
                if ( map[icopy]>=0 && records[map[icopy]]==irec ) smpl_gt[icopy] = bcf_gt_unphased(alt_indices[map[icopy]]);
                else if ( nmissing ) { smpl_gt[icopy] = bcf_gt_missing; nmissing--; }
                else smpl_gt[icopy] = bcf_gt_unphased(0);
            }
        }
    }
}

static bcf1_t *_abuf_new_rec(abuf_t *buf)
{
    rbuf_expand0(&buf->rbuf, bcf1_t*, buf->rbuf.n+1, buf->vcf);
    int ibuf = rbuf_append(&buf->rbuf);
    if ( !buf->vcf[ibuf] ) buf->vcf[ibuf] = bcf_init1();
    bcf1_t *out = buf->vcf[ibuf];
    bcf_clear1(out);
    return out;
}

// Copy a new ALT sequence into a reusable buffer, or return its existing index.
static int _abuf_join_add_alt(join_t *blk, const char *seq)
{
    int ial;
    for (ial=0; ial<blk->nalt; ial++)
        if ( !strcasecmp(blk->alt[ial].s,seq) ) return ial+1;
    hts_expand0(kstring_t,blk->nalt+1,blk->malt,blk->alt);
    blk->alt[blk->nalt].l = 0;
    kputs(seq,&blk->alt[blk->nalt]);
    return ++blk->nalt;
}

// no-genotype case: greedily concatenate ALTs in position order, fill gaps with
// reference bases. When two variants overlap, start a new ALT
static void _abuf_join_greedy(join_t *blk)
{
    // the joined REF covers [new_pos, new_pos + ref_len), end excluded.
    hts_pos_t new_pos = blk->records[0].rec->pos;   // 0-based
    const char *new_ref = blk->ref.s;       // new_ref[0] is at new_pos
    int ref_len = blk->ref.l;
    int irec, ial, has_alt = 0;
    kstring_t *alt_seq = &blk->seq;
    alt_seq->l = 0;
    hts_pos_t ref_end = new_pos;            // first ref pos not yet represented in alt_seq
    for (irec=0; irec<blk->nrec; irec++)
    {
        bcf1_t *rec = blk->records[irec].rec;
        for (ial=1; ial<rec->n_allele; ial++)
        {
            if ( !strcmp(rec->d.allele[ial],"*") )
                error("Error: cannot join star alleles without FORMAT/GT at position %"PRIhts_pos"\n",rec->pos+1);

            if ( has_alt && rec->pos<ref_end )  // another ALT of the same record also overlaps its preceding ALT
            {
                // append the remaining reference suffix, save the completed ALT
                kputsn(new_ref+(ref_end-new_pos), ref_len-(ref_end-new_pos), alt_seq);

                // add the new ALT unless identical to REF
                if ( strcasecmp(alt_seq->s,new_ref) ) _abuf_join_add_alt(blk,alt_seq->s);

                // clear ALT seq and reset ref_end to begin another ALT spanning the same joined block
                alt_seq->l = 0;
                ref_end = new_pos;
            }

            // copy unchanged REF up until the current variant
            kputsn(new_ref+(ref_end-new_pos), rec->pos-ref_end, alt_seq);

            // replace this record's entire REF span with the selected ALT
            kputs(rec->d.allele[ial], alt_seq);
            ref_end = rec->pos + blk->records[irec].ref_len;
            has_alt = 1;
        }
    }
    if ( has_alt ) // flush the last sequence
    {
        kputsn(new_ref+(ref_end-new_pos), ref_len-(ref_end-new_pos), alt_seq);
        if ( strcasecmp(alt_seq->s,new_ref) ) _abuf_join_add_alt(blk,alt_seq->s);
    }
}

// reconstruct called haplotypes and deduplicate their sequences into the ALT list
static void _abuf_join_genotypes(abuf_t *buf)
{
    join_t *blk = &buf->join;
    hts_pos_t new_pos = blk->records[0].rec->pos;
    const char *new_ref = blk->ref.s;
    int ref_len = blk->ref.l;
    int nsmpl = bcf_hdr_nsamples(buf->hdr);
    int irec, icopy, ismpl, iedit;
    hts_expand(int,nsmpl,blk->msmpl_ploidy,blk->smpl_ploidy);
    hts_expand(uint8_t,nsmpl,blk->mis_phased,blk->is_phased);
    join_rec_t *rows = blk->records;
    int *smpl_ploidy = blk->smpl_ploidy;
    uint8_t *is_phased = blk->is_phased;
    memset(is_phased,1,nsmpl);
    int max_ploidy = 0;

    // load genotypes, determine samples' ploidies and phase
    for (irec=0; irec<blk->nrec; irec++)
    {
        bcf1_t *rec = blk->records[irec].rec;
        int ngt = bcf_get_genotypes(buf->hdr,rec,&rows[irec].gt,&rows[irec].mgt);
        if ( ngt<=0 || !nsmpl || ngt%nsmpl )
            error("Error: could not read FORMAT/GT at %s:%"PRIhts_pos"\n",bcf_seqname(buf->hdr,rec),rec->pos+1);
        rows[irec].width = ngt/nsmpl;
        for (ismpl=0; ismpl<nsmpl; ismpl++)
        {
            int32_t *smpl_gt = rows[irec].gt + ismpl*rows[irec].width;
            int ncopy = 0;
            while ( ncopy<rows[irec].width && smpl_gt[ncopy]!=bcf_int32_vector_end )
            {
                if ( ncopy && !bcf_gt_is_phased(smpl_gt[ncopy]) ) is_phased[ismpl] = 0;
                if ( !bcf_gt_is_missing(smpl_gt[ncopy]) )
                {
                    int ial = bcf_gt_allele(smpl_gt[ncopy]);
                    if ( ial<0 || ial>=rec->n_allele )
                        error("Error: invalid GT allele index at %s:%"PRIhts_pos"\n",bcf_seqname(buf->hdr,rec),rec->pos+1);

                    // a star is an overlapping allele placeholder, never sequence; phased 0|* becomes 0|.
                    if ( !strcmp(rec->d.allele[ial],"*") ) smpl_gt[ncopy] = bcf_gt_is_phased(smpl_gt[ncopy]) ? bcf_gt_phased(-1) : bcf_gt_missing;
                }
                ncopy++;
            }
            if ( !ncopy || (irec && ncopy!=smpl_ploidy[ismpl]) )
                error("Error: inconsistent ploidy for sample %s at %s:%"PRIhts_pos"\n",buf->hdr->samples[ismpl],bcf_seqname(buf->hdr,rec),rec->pos+1);
            smpl_ploidy[ismpl] = ncopy;
            if ( ncopy>max_ploidy ) max_ploidy = ncopy;
        }
    }

    // mark duplicate VCF records (pos,ref,alt) so they don't mess up deatomization
    hts_expand(int32_t,max_ploidy,blk->malt_indices,blk->alt_indices);
    int duplicates_ready = 0;
    for (ismpl=0; ismpl<nsmpl; ismpl++)
    {
        if ( is_phased[ismpl] ) continue;
        if ( !duplicates_ready )
        {
            // run only once and only for unphased data, so that two 0/1 records with the same pos can become 1/2
            _abuf_join_find_duplicates(blk);
            duplicates_ready = 1;
        }
        // check consistency and reorder alleles in sample's GTs, favors cis placement
        _abuf_join_assign_unphased(buf,ismpl);
    }

    if ( max_ploidy>INT_MAX/nsmpl ) error("Error: too many genotype copies to join\n");
    blk->ngt = nsmpl*max_ploidy;
    hts_expand(int32_t,blk->ngt,blk->mjoined_gt,blk->joined_gt);
    hts_expand(join_edit_t,blk->nrec,blk->medits,blk->edits);
    hts_expand(uint8_t,ref_len,blk->mcovered,blk->covered);
    join_edit_t *edits = blk->edits;
    uint8_t *covered = blk->covered;
    kstring_t *seq = &blk->seq;
    for (ismpl=0; ismpl<nsmpl; ismpl++)
    {
        for (icopy=0; icopy<max_ploidy; icopy++)
        {
            int32_t *out_call = blk->joined_gt + (size_t)ismpl*max_ploidy+icopy;
            if ( icopy>=smpl_ploidy[ismpl] ) { *out_call = bcf_int32_vector_end; continue; }

            // check which ref bases are accounted for by the this chromosomal copy across records,
            // collect changes that need to be applied to ref
            int nedit = 0;
            memset(covered,0,ref_len);
            for (irec=0; irec<blk->nrec; irec++)
            {
                bcf1_t *rec = blk->records[irec].rec;
                int32_t gt_call = rows[irec].gt[ismpl*rows[irec].width+icopy];
                if ( bcf_gt_is_missing(gt_call) ) continue;
                int ref_offset = rec->pos-new_pos, record_ref_len = blk->records[irec].ref_len;
                memset(covered+ref_offset,1,record_ref_len);
                int ial = bcf_gt_allele(gt_call);
                if ( !ial ) continue;
                if ( _abuf_join_edit(rec->d.allele[0],blk->records[irec].ref_len,rec->d.allele[ial],ref_offset,&edits[nedit]) ) nedit++;
            }
            qsort(edits,nedit,sizeof(*edits),_abuf_join_edit_cmp);  // sort edits by position

            // build one chr copy from the edits
            seq->l = 0;
            int ref_end = 0, nunique = 0;   // ref_end: how much ref has been consumed; nunique: number of accepted edits
            for (iedit=0; iedit<nedit; iedit++)
            {
                join_edit_t *edit = &edits[iedit];
                join_edit_t *prev = nunique ? &edits[nunique-1] : NULL;

                // skip duplicates
                int same_span = prev && edit->beg==prev->beg && edit->end==prev->end;
                if ( same_span && edit->len==prev->len && !strncasecmp(edit->alt,prev->alt,edit->len) ) continue;

                // Overlapping substitutions can share bases if they agree.
                // Keep only the suffix not already represented in the sequence.
                if ( same_span || edit->beg<ref_end )
                {
                    if ( _abuf_join_conflicts(edits,nunique,edit) )
                        error("Error: conflicting overlapping alleles for sample %s at %s:%"PRIhts_pos"\n",
                              buf->hdr->samples[ismpl],bcf_seqname(buf->hdr,blk->records[0].rec),new_pos+edit->beg+1);
                    if ( edit->end<=ref_end ) continue;
                    int skip = ref_end-edit->beg;
                    edit->beg += skip;
                    edit->alt += skip;
                    edit->len -= skip;
                }

                // process the accepted edit
                kputsn(new_ref+ref_end,edit->beg-ref_end,seq);  // unchanged REF before the edit
                kputsn(edit->alt,edit->len,seq);                // replacement sequence
                ref_end = edit->end;                            // skip replaced REF
                edits[nunique++] = *edit;
            }
            kputsn(new_ref+ref_end,ref_len-ref_end,seq);        // append ref after the last edit

            if ( !seq->l )
                error("Error: empty joined allele for sample %s at %s:%"PRIhts_pos"\n", buf->hdr->samples[ismpl],bcf_seqname(buf->hdr,blk->records[0].rec),new_pos+1);

            // are there any gaps in ref allele coverage?
            // - if a non-star ALT candidate from a record with a missing genotype conflicts a called ALT edit on the same copy, mark it as covered
            int complete = !memchr(covered,0,ref_len);
            if ( !complete )
            {
                for (irec=0; irec<blk->nrec; irec++)
                {
                    if ( !bcf_gt_is_missing(rows[irec].gt[ismpl*rows[irec].width+icopy]) ) continue;
                    bcf1_t *rec = blk->records[irec].rec;
                    int ial, ref_offset = rec->pos-new_pos;
                    for (ial=1; ial<rec->n_allele; ial++)
                    {
                        if ( !strcmp(rec->d.allele[ial],"*") ) continue;    // '*' is a placeholder, not an actual replacement sequence
                        join_edit_t query;
                        if ( !_abuf_join_edit(rec->d.allele[0],blk->records[irec].ref_len,rec->d.allele[ial],ref_offset,&query) ) continue;
                        if ( !_abuf_join_conflicts(edits,nunique,&query) ) continue;

                        // ial of this row overlaps and therefore it can explain the gap
                        memset(covered+ref_offset,1,blk->records[irec].ref_len);
                        break;
                    }
                }
                complete = !memchr(covered,0,ref_len);
            }

            int ial = -1;   // if there is still a gap, missing allele will be used for this chr copy
            if ( complete )
                ial = !strcasecmp(seq->s,new_ref) ? 0 : _abuf_join_add_alt(blk,seq->s);
            *out_call = icopy && is_phased[ismpl] ? bcf_gt_phased(ial) : bcf_gt_unphased(ial);
        }
    }
}

// Finalize the block and queue its output. Singletons pass through unchanged.
static void _abuf_join_close(abuf_t *buf)
{
    join_t *blk = &buf->join;
    int irec, ial;

    if ( blk->nrec<=0 ) return;
    if ( blk->nrec==1 )
    {
        bcf1_t *out = _abuf_new_rec(buf);
        bcf_copy(out, blk->records[0].rec);
        bcf_unpack(out, BCF_UN_ALL);
        blk->nrec = 0;
        return;
    }

    // join receives records sorted by their normalized positions.
    bcf1_t *first_rec = blk->records[0].rec;
    hts_pos_t new_pos = first_rec->pos;
    if ( blk->end<=new_pos || blk->end-new_pos>=INT_MAX )
        error("Error: invalid joined REF span at %s:%"PRIhts_pos"\n",bcf_seqname(buf->hdr,first_rec),new_pos+1);
    int ref_len = blk->end-new_pos;

    // Only GT is copied to the merged record; warn once about discarded fields.
    static int warned = 0;
    int has_gt = bcf_get_fmt(buf->hdr,first_rec,"GT")!=NULL;
    for (irec=0; irec<blk->nrec; irec++)
    {
        bcf1_t *rec = blk->records[irec].rec;
        int rec_has_gt = bcf_get_fmt(buf->hdr,rec,"GT")!=NULL;
        if ( rec_has_gt!=has_gt )
            error("Error: records to be joined must either all have FORMAT/GT or all lack it, "
                    "see %s:%"PRIhts_pos" and %s:%"PRIhts_pos"\n",bcf_seqname(buf->hdr,first_rec),first_rec->pos+1,bcf_seqname(buf->hdr,rec),rec->pos+1);
        if ( rec->n_fmt>rec_has_gt && !warned )
        {
            fprintf(stderr,"Warning: dropping FORMAT fields other than GT from joined records\n");
            warned = 1;
        }
    }

    // Overlapping or touching REF spans cover the block without gaps.
    // Copy their bases, checking that overlapping records agree.
    if ( ks_resize(&blk->ref,(size_t)ref_len+1)<0 ) error("Error: could not allocate joined REF sequence\n");
    blk->ref.l = ref_len;
    char *ref = blk->ref.s;
    memset(ref,0,(size_t)ref_len+1);
    for (irec=0; irec<blk->nrec; irec++)
    {
        bcf1_t *rec = blk->records[irec].rec;
        int ref_offset = rec->pos - new_pos;
        int ibase;
        for (ibase=0; rec->d.allele[0][ibase]; ibase++)
        {
            char base = rec->d.allele[0][ibase];
            if ( ref[ref_offset+ibase] && toupper_c(ref[ref_offset+ibase])!=toupper_c(base) )
                error("Error: inconsistent REF at %s:%"PRIhts_pos"\n",bcf_seqname(buf->hdr,rec),rec->pos+ibase+1);
            ref[ref_offset+ibase] = base;
        }
    }

    blk->nalt = blk->ngt = 0;
    if ( has_gt )
        _abuf_join_genotypes(buf);
    else
        _abuf_join_greedy(blk);

    bcf1_t *out = _abuf_new_rec(buf);
    out->rid = first_rec->rid;
    out->pos = new_pos;
    bcf_update_id(buf->out_hdr, out, ".");

    kstring_t *allele_str = &blk->allele_str;
    allele_str->l = 0;
    kputs(ref, allele_str);
    for (ial=0; ial<blk->nalt; ial++)
    {
        kputc(',',allele_str);
        kputs(blk->alt[ial].s,allele_str);
    }
    if ( !blk->nalt ) kputs(",.", allele_str);
    bcf_update_alleles_str(buf->out_hdr, out, allele_str->s);

    bcf_float_set_missing(out->qual);
    bcf_update_filter(buf->out_hdr, out, NULL, 0);
    if ( has_gt ) bcf_update_genotypes(buf->out_hdr, out, blk->joined_gt, blk->ngt);
    else out->n_sample = first_rec->n_sample;

    blk->njoined++;
    blk->nrec = 0;
}

static void _abuf_join(abuf_t *buf, bcf1_t *rec)
{
    join_t *blk = &buf->join;
    hts_pos_t ref_len = strlen(rec->d.allele[0]);
    hts_pos_t end = rec->pos + ref_len;   // 0-based, exclusive

    if ( blk->nrec && rec->rid==blk->records[blk->nrec-1].rec->rid && rec->pos<blk->records[blk->nrec-1].rec->pos )
        error("Error: join requires position-sorted records, see %s:%"PRIhts_pos"\n",bcf_seqname(buf->hdr,rec),rec->pos+1);

    // SNVs and MNVs can extend a substitution block; other variant types end it.
    int can_join = 1;
    if ( blk->snvs_only ) can_join = bcf_has_variant_types(rec, VCF_SNP|VCF_MNP, bcf_match_subset) > 0;
    int ial;
    for (ial=0; can_join && ial<rec->n_allele; ial++)
    {
        if ( ial && !strcmp(rec->d.allele[ial],"*") ) continue;
        if ( strspn(rec->d.allele[ial],"ACGTNacgtn")!=(ial ? strlen(rec->d.allele[ial]) : ref_len) ) can_join = 0;
    }
    if ( !can_join )
    {
        _abuf_join_close(buf);
        bcf_copy(_abuf_new_rec(buf),rec);
        return;
    }

    if ( blk->nrec && (rec->rid!=blk->records[0].rec->rid || rec->pos>blk->end) ) _abuf_join_close(buf);

    hts_expand0(join_rec_t,blk->nrec+1,blk->mrec,blk->records);
    blk->records[blk->nrec].ref_len = ref_len;
    if ( !blk->records[blk->nrec].rec ) blk->records[blk->nrec].rec = bcf_init1();
    bcf_copy(blk->records[blk->nrec].rec, rec);
    bcf_unpack(blk->records[blk->nrec].rec, BCF_UN_ALL);
    if ( !blk->nrec || end>blk->end ) blk->end = end;
    blk->nrec++;
}

void abuf_push(abuf_t *buf, bcf1_t *rec)
{
    bcf_unpack(rec, BCF_UN_ALL);
    if ( buf->mode==SPLIT ) _abuf_split(buf,rec);
    else if ( buf->mode==JOIN ) _abuf_join(buf,rec);
}

void *abuf_get(abuf_t *buf, abuf_opt_t key)
{
    switch (key)
    {
        case N_JOINED: return &buf->join.njoined;
        default: hts_log_error("Unsupported abuf_get key=%d",(int)key); return NULL;
    }
}

bcf1_t *abuf_flush(abuf_t *buf, int flush_all)
{
    if ( flush_all && buf->mode==JOIN && buf->join.nrec ) _abuf_join_close(buf);

    if ( buf->rbuf.n==0 ) return NULL;
    int i = rbuf_shift(&buf->rbuf);
    return buf->vcf[i];
}

