/*  vcfnorm.c -- Left-align and normalize indels.

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
#include <htslib/synced_bcf_reader.h>
#include <htslib/faidx.h>
#include "bcftools.h"
#include "rbuf.h"

#define CHECK_REF_EXIT 0
#define CHECK_REF_WARN 1
#define CHECK_REF_SKIP 2

#define MROWS_SPLIT 1
#define MROWS_MERGE  2

typedef struct
{
    int32_t dir:4, val:28;
}
cell_t;
typedef struct
{
    int nmat, nref, nseq;
    int ipos, lref, lseq;
    cell_t *mat;
    char *ref, *seq;
    int m_arr, *ipos_arr, *lref_arr, *lseq_arr;
}
aln_aux_t;

// for -m+, mapping from allele indexes of a single input record
// to allele indexes of output record
typedef struct
{
    int nals, mals, *map;
}
map_t;

typedef struct
{
    aln_aux_t aln;
    char *tseq, *seq;
    int mseq;
    bcf1_t **lines, **tmp_lines, **alines, **blines, *mrow_out;
    int ntmp_lines, mtmp_lines, nalines, malines, nblines, mblines;
    map_t *maps;     // mrow map for each buffered record
    char **als;
    int mmaps, nals, mals;
    uint8_t *tmp_arr1, *tmp_arr2;
    int ntmp_arr1, ntmp_arr2;
    kstring_t *tmp_str;
    kstring_t *tmp_als, tmp_als_str;
    int ntmp_als;
    rbuf_t rbuf;
    int buf_win;            // maximum distance between two records to consider
    int aln_win;            // the realignment window size (maximum repeat size)
    bcf_srs_t *files;       // using the synced reader only for -r option
    bcf_hdr_t *hdr;
    faidx_t *fai;
    char **argv, *output_fname, *ref_fname, *vcf_fname, *region, *targets;
    int argc, rmdup, output_type, check_ref, strict_filter;
    int nchanged, nskipped, ntotal, mrows_op, mrows_collapse, parsimonious;
}
args_t;

#if OLD_WAY
void _vcfnorm_debug_print(aln_aux_t *aux)
{
    cell_t *mat = aux->mat;
    char *ref = aux->ref;
    char *seq = aux->seq;
    int nref  = aux->nref;
    int nseq  = aux->nseq;
    int k     = (nref+1)*(nseq+1)-1;
    int kd    = nref+2;
    int i = k/(nref+1);
    int j = k - i*(nref+1);
    assert(i>0 && j>0);
    int l = k, ialn = 0, nout_ref = 0, nout_seq = 0, ipos = 0;
    char *aln_ref = (char*) malloc(sizeof(char)*(k+1));
    char *aln_seq = (char*) malloc(sizeof(char)*(k+1));
    while ( l>0 )
    {
        if ( j<=0 || mat[l].dir==1 )    // i
        {
            aln_ref[ialn] = '-';
            aln_seq[ialn] = seq[nseq-i];
            ipos = 0;
            nout_seq++;
            l -= kd - 1;
            i--;
        }
        else if ( i<=0 || mat[l].dir==-1 )  // d
        {
            aln_ref[ialn] = ref[nref-j];
            aln_seq[ialn] = '-';
            ipos = 0;
            nout_ref++;
            l--;
            j--;
        }
        else     // match or mismatch
        {
            aln_ref[ialn] = ref[nref-j];
            aln_seq[ialn] = seq[nseq-i];
            ipos = ref[nref-j+1]==seq[nseq-i+1] ? ipos+1 : 1;
            nout_seq++;
            nout_ref++;
            l -= kd;
            i--;
            j--;
        }
        ialn++;
    }
    aln_ref[ialn] = aln_seq[ialn] = 0;
    fprintf(stderr, "ref: %s\n", ref);
    fprintf(stderr, "seq: %s\n", seq);
    fprintf(stderr, "-> %s\n", aln_ref);
    fprintf(stderr, "-> %s\n", aln_seq);
    free(aln_ref);
    free(aln_seq);

    fprintf(stderr, "      ");
    for (j=0; j<nref; j++) fprintf(stderr, "   %c ", ref[nref-j-1]); fprintf(stderr, "\n");
    for (i=0; i<=nseq; i++)
    {
        fprintf(stderr, "%c", i==0 ? ' ' : seq[nseq-i]);
        for (j=0; j<=nref; j++)
        {
            char dir = ' ';
            if ( mat[i*(nref+1)+j].dir==1 ) dir = 'i';
            else if ( mat[i*(nref+1)+j].dir==-1 ) dir = 'd';
            fprintf(stderr, " %3d%c", (int)mat[i*(nref+1)+j].val, dir);
        }
        fprintf(stderr,"\n");
    }
}

static int align(args_t *args, aln_aux_t *aux)
{
    // Needleman-Wunsch global alignment. Note that the sequences are aligned from
    //  the end where matches are preferred, gaps are pushed to the front (left-aligned)
    char *ref = aux->ref;
    char *seq = aux->seq;
    int nref  = aux->nref;
    int nseq  = aux->nseq;
    if ( (nref+1)*(nseq+1) > aux->nmat )
    {
        aux->nmat = (nref+1)*(nseq+1);
        aux->mat  = (cell_t *) realloc(aux->mat, sizeof(cell_t)*aux->nmat);
        if ( !aux->mat )
            error("Could not allocate %ld bytes of memory at %d\n", sizeof(cell_t)*aux->nmat, args->files->readers[0].buffer[0]->pos+1);
    }
    const int GAP_OPEN = -1, GAP_CLOSE = -1, GAP_EXT = 0, MATCH = 1, MISM = -1, DI = 1, DD = -1, DM = 0;
    cell_t *mat = aux->mat;
    int i, j, k = nref+2, kd = nref+2;
    mat[0].val = 20; mat[0].dir = DM;   // the last ref and alt bases match
    for (j=1; j<=nref; j++) { mat[j].val = 0; mat[j].dir = DM; }
    for (i=1; i<=nseq; i++)
    {
        mat[k-1].val = 0;
        mat[k-1].dir = DM;
        int jmax = i-1 < nref ? i-1 : nref;
        for (j=1; j<=jmax; j++)
        {
            // prefer insertions to deletions and mismatches
            int max, dir, score;
            if ( ref[nref-j]==seq[nseq-i] )
            {
                // match
                max = mat[k-kd].val + MATCH;
                if ( mat[k-kd].dir!=DM ) max += GAP_CLOSE;
                dir = DM;

                // insertion
                score = mat[k-kd+1].dir == DI ? mat[k-kd+1].val + GAP_EXT: mat[k-kd+1].val + GAP_OPEN;
                if ( max < score )  { max = score; dir = DI; }

                // deletion
                score = mat[k-1].dir == DD ? mat[k-1].val + GAP_EXT : mat[k-1].val + GAP_OPEN;
                if ( max < score ) { max = score; dir = DD; }
            }
            else
            {
                // insertion
                max = mat[k-kd+1].dir == DI ? mat[k-kd+1].val + GAP_EXT : mat[k-kd+1].val + GAP_OPEN;
                dir = DI;

                // deletion
                score = mat[k-1].dir == DD ? mat[k-1].val + GAP_EXT : mat[k-1].val + GAP_OPEN;
                if ( max < score ) { max = score; dir = DD; }

                // mismatch
                score = mat[k-kd].val + MISM;
                if ( mat[k-kd].dir!=DM ) score += GAP_CLOSE;
                if ( max < score ) { max = score; dir = DM; }
            }
            mat[k].val = max;
            mat[k].dir = dir;
            k++;
        }
        for (j=jmax+1; j<=nref; j++)
        {
            // prefer deletions to insertions and mismatches
            int max, dir, score;
            if ( ref[nref-j]==seq[nseq-i] )
            {
                // match
                max = mat[k-kd].val + MATCH;
                if ( mat[k-kd].dir!=DM ) max += GAP_CLOSE;
                dir = DM;

                // deletion
                score = mat[k-1].dir == DD ? mat[k-1].val + GAP_EXT: mat[k-1].val + GAP_OPEN;
                if ( max < score ) { max = score; dir = DD; }

                // insertion
                score = mat[k-kd+1].dir == DI ? mat[k-kd+1].val + GAP_EXT : mat[k-kd+1].val + GAP_OPEN;
                if ( max < score )  { max = score; dir = DI; }
            }
            else
            {
                // deletion
                max = mat[k-1].dir == DD ? mat[k-1].val + GAP_EXT : mat[k-1].val + GAP_OPEN;
                dir = DD;

                // insertion
                score = mat[k-kd+1].dir == DI ? mat[k-kd+1].val + GAP_EXT : mat[k-kd+1].val + GAP_OPEN;
                if ( max < score ) { max = score; dir = DI; }

                // mismatch
                score = mat[k-kd].val + MISM;
                if ( mat[k-kd].dir!=DM ) score += GAP_CLOSE;
                if ( max < score ) { max = score; dir = DM; }
            }
            mat[k].val = max;
            mat[k].dir = dir;
            k++;
        }
        k++;
    }

    // _vcfnorm_debug_print(aux);

    // Skip as much of the matching sequence at the beggining as possible. (Note, the sequence
    // is reversed, thus skipping from the end.)
    k = (nref+1)*(nseq+1)-1;
    int kmin = nref>nseq ? 2*(nref+1) - nseq : (nseq-nref)*(nref+1);    // skip the first row and column of the matrix, which are 0s
    int ipos = 0;
    while (k>kmin && mat[k].dir==DM && ref[ipos]==seq[ipos]) { k -= kd; ipos++; }

    i = k/(nref+1);             // seq[nseq-i]
    j = k - i*(nref+1);         // ref[nref-j]

    if ( !i && !j )
    {
        // no gaps: this is a legitimate case, consider MNPs
        int l = 0;
        while ( l<nseq && l<nref && ref[l]==seq[l] ) l++;
        if ( l==nseq || l==nref ) return -2;    // ALT is identical to REF
        assert( l>0 );
        while ( ipos>l && ref[ipos]==seq[ipos] ) ipos--;
        aux->ipos = l;          // position of the last matching base (buffer coordinates)
        aux->lref = ipos + 1;   // position of the first base in the suffix
        aux->lseq = ipos + 1;
        return 0;
    }
    assert(i>0 && j>0);

    int l = k, nout_ref = ipos, nout_seq = ipos, nsuffix = 0;
    while ( l>0 )
    {
        if ( j<=0 || mat[l].dir==DI )    // insertion
        {
            nsuffix = 0;
            nout_seq++;
            l -= kd - 1;
            i--;
        }
        else if ( i<=0 || mat[l].dir==DD )  // deletion
        {
            nsuffix = 0;
            nout_ref++;
            l--;
            j--;
        }
        else     // match/mismatch
        {
            nsuffix = ref[nref-j]==seq[nseq-i] ? nsuffix + 1 : 0;
            nout_seq++;
            nout_ref++;
            l -= kd;
            i--;
            j--;
        }
    }
    if ( !ipos ) return -1; // the window is too small

    aux->ipos = ipos - 1;
    aux->lref = nout_ref - nsuffix;
    aux->lseq = nout_seq - nsuffix;

    // The indels and complex events do not have to be padded
    if ( aux->lref - aux->ipos > 1 && aux->lseq - aux->ipos > 1 && ref[aux->ipos]==seq[aux->ipos] ) aux->ipos++;
    return 0;
}

int realign(args_t *args, bcf1_t *line)
{
    bcf_unpack(line, BCF_UN_STR);

    char *tmp;
    int i, ref_len = strlen(line->d.allele[0]), len = ref_len;
    for (i=1; i<line->n_allele; i++)
    {
        int l = strlen(line->d.allele[i]);
        if ( len < l ) len = l;
    }
    for (i=0; i<line->n_allele; i++)
    {
        tmp = line->d.allele[i];
        while (*tmp) { *tmp = toupper(*tmp); tmp++; }
    }

    int ref_winlen;
    char *ref = NULL;
    if ( len==1 || line->n_allele==1 )
    {
        // SNP
        int reflen = strlen(line->d.allele[0]);
        ref = faidx_fetch_seq(args->fai, (char*)args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos, line->pos+reflen-1, &ref_winlen);
        if ( !ref ) error("faidx_fetch_seq failed at %s:%d\n", args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1);
        if ( !strncmp(ref,line->d.allele[0],reflen) )
        {
            free(ref);
            return 0;
        }
        if ( args->check_ref==CHECK_REF_EXIT )
            error("Reference allele mismatch at %s:%d .. '%c' vs '%c'\n", bcf_seqname(args->hdr,line),line->pos+1,ref[0],line->d.allele[0][0]);
        if ( args->check_ref & CHECK_REF_WARN )
            fprintf(stderr,"REF_MISMATCH\t%s\t%d\t%s\n", bcf_seqname(args->hdr,line),line->pos+1,line->d.allele[0]);
        free(ref);
        return -1;
    }

    // Sanity check: exclude broken VCFs with long stretches of N's
    if ( len>1000 )
    {
        for (i=0; i<ref_len; i++)
            if ( line->d.allele[0][i]=='N' ) return -1;
    }

    int win = line->pos < args->aln_win ? line->pos : args->aln_win;
    len += win + 2;
    if ( args->mseq < len*(line->n_allele-1) )
    {
        args->mseq = len*(line->n_allele-1);
        args->seq  = (char*) realloc(args->seq, sizeof(char)*args->mseq);
    }
    ref = faidx_fetch_seq(args->fai, (char*)args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos-win, line->pos+ref_len, &ref_winlen);
    if ( !ref ) error("faidx_fetch_seq failed at %s:%d\n", args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos-win);
    assert( ref_winlen==ref_len+win+1 );

    for (i=0; i<ref_len; i++) ref[i] = toupper(ref[i]);

    // Sanity check: the reference sequence must match the REF allele
    if ( strncmp(&ref[win],line->d.allele[0],ref_len) )
    {
        for (i=0; i<ref_len; i++)
            if ( ref[win+i]!=line->d.allele[0][i] ) break;
        if ( args->check_ref==CHECK_REF_EXIT )
            error("\nSanity check failed, the reference sequence differs at %s:%d[%d] .. '%c' vs '%c'\n",
                    args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1, i+1,ref[win+i],line->d.allele[0][i]);
        if ( args->check_ref & CHECK_REF_WARN )
            fprintf(stderr,"REF_MISMATCH\t%s\t%d\t%s\n", bcf_seqname(args->hdr,line),line->pos+1,line->d.allele[0]);
        free(ref);
        return -1;
    }

    if ( args->aln.m_arr < line->n_allele )
    {
        args->aln.m_arr = line->n_allele;
        args->aln.ipos_arr = (int*) realloc(args->aln.ipos_arr, sizeof(int)*args->aln.m_arr);
        args->aln.lref_arr = (int*) realloc(args->aln.lref_arr, sizeof(int)*args->aln.m_arr);
        args->aln.lseq_arr = (int*) realloc(args->aln.lseq_arr, sizeof(int)*args->aln.m_arr);
    }
    int *ipos = args->aln.ipos_arr;
    int *iref = args->aln.lref_arr;
    int *iseq = args->aln.lseq_arr;
    int min_pos = INT_MAX, max_ref = 0;

    static int aln_win_warned = 0;
    int j, k;
    for (j=1; j<line->n_allele; j++)
    {
        // get the ALT ready for alignment
        args->tseq = args->seq + (j-1)*len;
        for (i=0; i<win; i++) args->tseq[i] = ref[i];
        char *t = line->d.allele[j];
        while (*t) { args->tseq[i++] = *t; t++; }
        args->tseq[i++] = ref[ref_winlen-1];
        args->tseq[i]   = 0;

        args->aln.ref  = ref;
        args->aln.seq  = args->tseq;
        args->aln.nref = ref_winlen;    // reflen which goes into realignment: VCF ref + #win bases which precede
        args->aln.nseq = i;             // same as reflen but for alt

        int ret = align(args, &args->aln);
        if ( ret<0 )
        {
            free(ref);
            if ( ret==-1 )
            {
                // todo: better error analysis, see 2:1 in test/norm.vcf
                // fprintf(stderr,"Warning: The -w alignment window too small for %s:%d\n", bcf_seqname(args->hdr,line),line->pos+1);
                return 0;   // leave as is
            }
            if ( ret==-2 ) fprintf(stderr,"Warning: REF allele is identical to ALT at %s:%d\n", bcf_seqname(args->hdr,line),line->pos+1);
            return -1;
        }

        #if 0
            fprintf(stderr, "%s  \t nref=%d\n", ref, ref_winlen);
            fprintf(stderr, "%s  \t nseq=%d\n", args->tseq, i);
            fprintf(stderr, "pos=%d win=%d  ipos=%d lref=%d lseq=%d\n", line->pos+1, win, args->aln.ipos, args->aln.lref, args->aln.lseq);
            fprintf(stderr, "-> "); for (k=args->aln.ipos; k<args->aln.lref; k++) fprintf(stderr, "%c", ref[k]); fprintf(stderr, "\n");
            fprintf(stderr, "-> "); for (k=args->aln.ipos; k<args->aln.lseq; k++) fprintf(stderr, "%c", args->tseq[k]); fprintf(stderr, "\n");
            fprintf(stderr, "\n");
        #endif

        if ( !args->aln.ipos && !aln_win_warned )
        {
            fprintf(stderr,"Warning: bigger -w is needed [todo: improve me, %s:%d %s,%s]\n",
                    bcf_seqname(args->hdr,line),line->pos+1, line->d.allele[0],line->d.allele[j]);
            aln_win_warned = 1;
        }
        ipos[j] = args->aln.ipos;   // 0-based position before the first difference (w.r.t. the window)
        iref[j] = args->aln.lref;   // length of the REF alignment (or the index after the last aligned position)
        iseq[j] = args->aln.lseq;
        if ( max_ref < iref[j] ) max_ref = iref[j];
        if ( min_pos > ipos[j] ) min_pos = ipos[j];
        assert( iseq[j]<=len );
    }

    // Check if the record's position must be changed
    int nmv = win - min_pos;
    // This assertion that we will never want align more to the right is too
    // strong, consider cases like GATG -> GACT
    //  assert( nmv>=0 );
    line->pos -= nmv;

    hts_expand0(kstring_t,line->n_allele,args->ntmp_als,args->tmp_als);
    // REF
    args->tmp_als[0].l = 0;
    kputsn(&ref[min_pos], max_ref-min_pos, &args->tmp_als[0]);
    // ALTs
    int min_len = args->tmp_als[0].l;
    for (k=1; k<line->n_allele; k++)
    {
        args->tmp_als[k].l = 0;

        // prefix the sequence with REF bases if the other alleles were aligned more to the left
        int nprefix = ipos[k] - min_pos;
        if ( nprefix ) kputsn(&ref[min_pos], nprefix, &args->tmp_als[k]);

        // the ALT sequence
        int nseq = iseq[k] - ipos[k];
        if ( nseq )
        {
            char *alt = args->seq + (k-1)*len + ipos[k];
            kputsn(alt, nseq, &args->tmp_als[k]);
        }

        // suffix invoked by other deletions which must be added to match the REF
        int nsuffix = max_ref - iref[k];
        if ( nsuffix ) kputsn(&ref[iref[k]], nsuffix, &args->tmp_als[k]);

        if ( min_len > args->tmp_als[k].l ) min_len = args->tmp_als[k].l;
    }
    free(ref);

    // create new block of alleles
    args->tmp_als_str.l = 0;
    if ( args->parsimonious )
    {
        // check if we can trim from the left
        for (i=0; i<args->tmp_als[0].l; i++)
        {
            for (k=1; k<line->n_allele; k++)
                if ( args->tmp_als[0].s[i]!=args->tmp_als[k].s[i] ) break;
            if ( k!=line->n_allele ) break;
        }
        if ( i>=min_len ) i = min_len - 1;
        if ( i>0 )  // can be left trimmed
        {
            for (k=0; k<line->n_allele; k++)
            {
                if (k>0) kputc(',',&args->tmp_als_str);
                kputsn(args->tmp_als[k].s+i,args->tmp_als[k].l-i,&args->tmp_als_str);
            }
            kputc(0,&args->tmp_als_str);
            line->pos += i;
        }
        // if REF allele has not changed, ALTs must be unchanged as well
        else if ( !strcmp(line->d.allele[0],args->tmp_als[0].s) ) return 1;
    }
    else if ( !strcmp(line->d.allele[0],args->tmp_als[0].s) ) return 1;
    if ( !args->tmp_als_str.l )
    {
        for (k=0; k<line->n_allele; k++)
        {
            if (k>0) kputc(',',&args->tmp_als_str);
            kputsn(args->tmp_als[k].s,args->tmp_als[k].l,&args->tmp_als_str);
        }
    }
    args->nchanged++;
    bcf_update_alleles_str(args->hdr,line,args->tmp_als_str.s);

    return 1;
}
#else

static int realign(args_t *args, bcf1_t *line)
{
    bcf_unpack(line, BCF_UN_STR);

    // Sanity check REF
    int nref, reflen = strlen(line->d.allele[0]);
    char *ref = faidx_fetch_seq(args->fai, (char*)args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos, line->pos+reflen-1, &nref);
    if ( !ref ) error("faidx_fetch_seq failed at %s:%d\n", args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1);
    if ( strcasecmp(ref,line->d.allele[0]) )
    {
        if ( args->check_ref==CHECK_REF_EXIT )
            error("Reference allele mismatch at %s:%d .. '%s' vs '%s'\n", bcf_seqname(args->hdr,line),line->pos+1,ref,line->d.allele[0]);
        if ( args->check_ref & CHECK_REF_WARN )
            fprintf(stderr,"REF_MISMATCH\t%s\t%d\t%s\n", bcf_seqname(args->hdr,line),line->pos+1,line->d.allele[0]);
        free(ref);
        return -1;
    }
    free(ref);
    ref = NULL;

    if ( line->n_allele == 1 ) return 0;    // a REF

    // make a copy of each allele for trimming
    int i;
    hts_expand0(kstring_t,line->n_allele,args->ntmp_als,args->tmp_als);
    kstring_t *als = args->tmp_als;
    for (i=0; i<line->n_allele; i++)
    {
        if ( line->d.allele[i][0]=='<' ) return 0;  // symbolic allele

        als[i].l = 0;
        kputs(line->d.allele[i], &als[i]);
    }


    // trim from right
    int ori_pos = line->pos;
    while (1)
    {
        // is the rightmost base identical in all alleles?
        for (i=1; i<line->n_allele; i++)
        {
            if ( als[0].s[ als[0].l-1 ]!=als[i].s[ als[i].l-1 ] ) break;
        }
        if ( i!=line->n_allele ) break; // there are differences, cannot be trimmed

        int pad_from_left = 0;
        for (i=0; i<line->n_allele; i++) // trim all alleles
        {
            als[i].l--;
            if ( !als[i].l ) pad_from_left = 1;
        }
        if ( pad_from_left )
        {
            int npad = line->pos >= args->aln_win ? args->aln_win : line->pos;
            free(ref);
            ref = faidx_fetch_seq(args->fai, (char*)args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos-npad, line->pos-1, &nref);
            if ( !ref ) error("faidx_fetch_seq failed at %s:%d\n", args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos-npad+1);
            for (i=0; i<line->n_allele; i++)
            {
                ks_resize(&als[i], als[i].l + npad);
                if ( als[i].l ) memmove(als[i].s+npad,als[i].s,als[i].l);
                memcpy(als[i].s,ref,npad);
                als[i].l += npad;
            }
            line->pos -= npad;
        }
    }
    free(ref);

    // trim from left
    int ntrim_left = 0;
    while (1)
    {
        // is the first base identical in all alleles?
        int min_len = als[0].l - ntrim_left;
        for (i=1; i<line->n_allele; i++)
        {
            if ( als[0].s[ntrim_left]!=als[i].s[ntrim_left] ) break;
            if ( min_len > als[i].l - ntrim_left ) min_len = als[i].l - ntrim_left;
        }
        if ( i!=line->n_allele || min_len==1 ) break; // there are differences, cannot be trimmed
        ntrim_left++;
    }
    if ( ntrim_left )
    {
        for (i=0; i<line->n_allele; i++)
        {
            memmove(als[i].s,als[i].s+ntrim_left,als[i].l-ntrim_left);
            als[i].l -= ntrim_left;
        }
        line->pos += ntrim_left;
    }

    // Have the alleles changed?
    als[0].s[ als[0].l ] = 0;  // in order for strcmp to work
    if ( ori_pos==line->pos && !strcasecmp(line->d.allele[0],als[0].s) ) return 1;

    // Create new block of alleles and update
    args->tmp_als_str.l = 0;
    for (i=0; i<line->n_allele; i++)
    {
        if (i>0) kputc(',',&args->tmp_als_str);
        kputsn(als[i].s,als[i].l,&args->tmp_als_str);
    }
    args->tmp_als_str.s[ args->tmp_als_str.l ] = 0;
    bcf_update_alleles_str(args->hdr,line,args->tmp_als_str.s);
    args->nchanged++;

    return 1;
}

#endif

static void split_info_numeric(args_t *args, bcf1_t *src, bcf_info_t *info, int ialt, bcf1_t *dst)
{
    #define BRANCH_NUMERIC(type,type_t) \
    { \
        const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key); \
        int ntmp = args->ntmp_arr1 / sizeof(type_t); \
        int ret = bcf_get_info_##type(args->hdr,src,tag,&args->tmp_arr1,&ntmp); \
        args->ntmp_arr1 = ntmp * sizeof(type_t); \
        assert( ret>0 ); \
        type_t *vals = (type_t*) args->tmp_arr1; \
        int len = bcf_hdr_id2length(args->hdr,BCF_HL_INFO,info->key); \
        if ( len==BCF_VL_A ) \
        { \
            assert( ret==src->n_allele-1); \
            bcf_update_info_##type(args->hdr,dst,tag,vals+ialt,1); \
        } \
        else if ( len==BCF_VL_R ) \
        { \
            assert( ret==src->n_allele); \
            if ( ialt!=0 ) vals[1] = vals[ialt+1]; \
            bcf_update_info_##type(args->hdr,dst,tag,vals,2); \
        } \
        else if ( len==BCF_VL_G ) \
        { \
            assert( ret==src->n_allele*(src->n_allele+1)/2 ); \
            if ( ialt!=0 ) \
            { \
                vals[1] = vals[bcf_alleles2gt(0,ialt+1)]; \
                vals[2] = vals[bcf_alleles2gt(ialt+1,ialt+1)]; \
            } \
            bcf_update_info_##type(args->hdr,dst,tag,vals,3); \
        } \
        else \
            bcf_update_info_##type(args->hdr,dst,tag,vals,ret); \
    }
    switch (bcf_hdr_id2type(args->hdr,BCF_HL_INFO,info->key))
    {
        case BCF_HT_INT:  BRANCH_NUMERIC(int32, int32_t); break;
        case BCF_HT_REAL: BRANCH_NUMERIC(float, float); break;
    }
    #undef BRANCH_NUMERIC
}
// Find n-th field in a comma-separated list and move it to dst.
// The memory areas may overlap.
#define STR_MOVE_NTH(dst,src,end,nth,len) \
{ \
    char *ss = src, *se = src; \
    int j = 0; \
    while ( *se && se<(end) ) \
    { \
        if ( *se==',' ) \
        { \
            if ( j==nth ) break; \
            j++; \
            ss = se+1; \
        } \
        se++; \
    } \
    if ( j==nth ) \
    { \
        int n = se - ss; \
        memmove((dst),ss,n); \
        src = se; \
        len += n; \
    } \
    else len = -1; \
}
static void split_info_string(args_t *args, bcf1_t *src, bcf_info_t *info, int ialt, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key);
    int ret = bcf_get_info_string(args->hdr,src,tag,&args->tmp_arr1,&args->ntmp_arr1);
    assert( ret>0 );

    kstring_t str;
    str.m = args->ntmp_arr1;
    str.l = ret;
    str.s = (char*) args->tmp_arr1;

    int len = bcf_hdr_id2length(args->hdr,BCF_HL_INFO,info->key);
    if ( len==BCF_VL_A )
    {
        char *tmp = str.s;
        int len = 0;
        STR_MOVE_NTH(str.s,tmp,str.s+str.l,ialt,len);
        if ( len<0 ) return;   // wrong number of fields: skip
        str.s[len] = 0;
        bcf_update_info_string(args->hdr,dst,tag,str.s);
    }
    else if ( len==BCF_VL_R )
    {
        char *tmp = str.s;
        int len = 0;
        STR_MOVE_NTH(str.s,tmp,str.s+str.l,0,len);
        str.s[len]=','; tmp++; len++;
        STR_MOVE_NTH(&str.s[len],tmp,str.s+str.l,ialt,len);
        if ( len<0 ) return;   // wrong number of fields: skip
        str.s[len] = 0;
        bcf_update_info_string(args->hdr,dst,tag,str.s);
    }
    else if ( len==BCF_VL_G )
    {
        int i0a = bcf_alleles2gt(0,ialt+1), iaa = bcf_alleles2gt(ialt+1,ialt+1);
        char *tmp = str.s;
        int len = 0;
        STR_MOVE_NTH(str.s,tmp,str.s+str.l,0,len);
        str.s[len]=','; tmp++; len++;
        STR_MOVE_NTH(&str.s[len],tmp,str.s+str.l,i0a-1,len);
        if ( len<0 ) return;   // wrong number of fields: skip
        str.s[len]=','; tmp++; len++;
        STR_MOVE_NTH(&str.s[len],tmp,str.s+str.l,iaa-i0a-1,len);
        if ( len<0 ) return;   // wrong number of fields: skip
        str.s[len] = 0;
        bcf_update_info_string(args->hdr,dst,tag,str.s);
    }
    else
        bcf_update_info_string(args->hdr,dst,tag,str.s);
}
static void split_info_flag(args_t *args, bcf1_t *src, bcf_info_t *info, int ialt, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key);
    int ret = bcf_get_info_flag(args->hdr,src,tag,&args->tmp_arr1,&args->ntmp_arr1);
    bcf_update_info_flag(args->hdr,dst,tag,NULL,ret);
}

static void split_format_genotype(args_t *args, bcf1_t *src, bcf_fmt_t *fmt, int ialt, bcf1_t *dst)
{
    int ntmp = args->ntmp_arr1 / 4;
    int ngts = bcf_get_genotypes(args->hdr,src,&args->tmp_arr1,&ntmp);
    args->ntmp_arr1 = ntmp * 4;
    assert( ngts >0 );

    int32_t *gt = (int32_t*) args->tmp_arr1;
    int i, j, nsmpl = bcf_hdr_nsamples(args->hdr);
    ngts /= nsmpl;
    for (i=0; i<nsmpl; i++)
    {
        for (j=0; j<ngts; j++)
        {
            if ( gt[j]==bcf_int32_vector_end ) break;
            if ( bcf_gt_is_missing(gt[j]) || bcf_gt_allele(gt[j])==0 ) continue; // missing allele or ref: leave as is
            if ( bcf_gt_allele(gt[j])==ialt+1 )
                gt[j] = bcf_gt_unphased(1) | bcf_gt_is_phased(gt[j]); // set to first ALT
            else
                gt[j] = bcf_gt_unphased(0) | bcf_gt_is_phased(gt[j]); // set to REF
        }
        gt += ngts;
    }
    bcf_update_genotypes(args->hdr,dst,args->tmp_arr1,ngts*nsmpl);
}
static void split_format_numeric(args_t *args, bcf1_t *src, bcf_fmt_t *fmt, int ialt, bcf1_t *dst)
{
    #define BRANCH_NUMERIC(type,type_t,is_vector_end,set_vector_end) \
    { \
        const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id); \
        int ntmp = args->ntmp_arr1 / sizeof(type_t); \
        int nvals = bcf_get_format_##type(args->hdr,src,tag,&args->tmp_arr1,&ntmp); \
        args->ntmp_arr1 = ntmp * sizeof(type_t); \
        assert( nvals>0 ); \
        type_t *vals = (type_t *) args->tmp_arr1; \
        int len = bcf_hdr_id2length(args->hdr,BCF_HL_FMT,fmt->id); \
        int i, nsmpl = bcf_hdr_nsamples(args->hdr); \
        if ( nvals==nsmpl ) /* all values are missing */ \
        { \
            bcf_update_format_##type(args->hdr,dst,tag,vals,nsmpl); \
            return; \
        } \
        if ( len==BCF_VL_A ) \
        { \
            assert( nvals==(src->n_allele-1)*nsmpl); \
            nvals /= nsmpl; \
            type_t *src_vals = vals, *dst_vals = vals; \
            for (i=0; i<nsmpl; i++) \
            { \
                dst_vals[0] = src_vals[ialt]; \
                dst_vals += 1; \
                src_vals += nvals; \
            } \
            bcf_update_format_##type(args->hdr,dst,tag,vals,nsmpl); \
        } \
        else if ( len==BCF_VL_R ) \
        { \
            assert( nvals==src->n_allele*nsmpl); \
            nvals /= nsmpl; \
            type_t *src_vals = vals, *dst_vals = vals; \
            for (i=0; i<nsmpl; i++) \
            { \
                dst_vals[0] = src_vals[0]; \
                dst_vals[1] = src_vals[ialt+1]; \
                dst_vals += 2; \
                src_vals += nvals; \
            } \
            bcf_update_format_##type(args->hdr,dst,tag,vals,nsmpl*2); \
        } \
        else if ( len==BCF_VL_G ) \
        { \
            if ( nvals!=src->n_allele*(src->n_allele+1)/2*nsmpl && nvals!=src->n_allele*nsmpl ) \
                error("Error at %s:%d, the tag %s has wrong number of fields\n", bcf_seqname(args->hdr,src),src->pos+1,bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id)); \
            nvals /= nsmpl; \
            int all_haploid = nvals==src->n_allele ? 1 : 0; \
            type_t *src_vals = vals, *dst_vals = vals; \
            for (i=0; i<nsmpl; i++) \
            { \
                int haploid = all_haploid; \
                if ( !haploid ) \
                { \
                    int j; \
                    for (j=0; j<nvals; j++) if ( is_vector_end ) break; \
                    if ( j!=nvals ) haploid = 1; \
                } \
                dst_vals[0] = src_vals[0]; \
                if ( haploid ) \
                { \
                    dst_vals[1] = src_vals[ialt+1]; \
                    if ( !all_haploid ) set_vector_end; \
                } \
                else \
                { \
                    dst_vals[1] = src_vals[bcf_alleles2gt(0,ialt+1)]; \
                    dst_vals[2] = src_vals[bcf_alleles2gt(ialt+1,ialt+1)]; \
                } \
                dst_vals += all_haploid ? 2 : 3; \
                src_vals += nvals; \
            } \
            bcf_update_format_##type(args->hdr,dst,tag,vals,all_haploid ? nsmpl*2 : nsmpl*3); \
        } \
        else \
            bcf_update_format_##type(args->hdr,dst,tag,vals,nvals); \
    }
    switch (bcf_hdr_id2type(args->hdr,BCF_HL_FMT,fmt->id))
    {
        case BCF_HT_INT:  BRANCH_NUMERIC(int32, int32_t, src_vals[j]==bcf_int32_vector_end, dst_vals[2]=bcf_int32_vector_end); break;
        case BCF_HT_REAL: BRANCH_NUMERIC(float, float, bcf_float_is_vector_end(src_vals[j]), bcf_float_set_vector_end(dst_vals[2])); break;
    }
    #undef BRANCH_NUMERIC
}
static void squeeze_format_char(char *str, int src_blen, int dst_blen, int n)
{
    int i, isrc = 0, idst = 0;
    for (i=0; i<n; i++)
    {
        memmove(str+idst,str+isrc,dst_blen);
        idst += dst_blen;
        isrc += src_blen;
    }
}
static void split_format_string(args_t *args, bcf1_t *src, bcf_fmt_t *fmt, int ialt, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id);
    int ret = bcf_get_format_char(args->hdr,src,tag,&args->tmp_arr1,&args->ntmp_arr1);
    assert( ret>0 );

    kstring_t str;
    str.m = args->ntmp_arr1;
    str.l = ret;
    str.s = (char*) args->tmp_arr1;

    int nsmpl = bcf_hdr_nsamples(args->hdr);
    int len = bcf_hdr_id2length(args->hdr,BCF_HL_FMT,fmt->id);
    if ( len==BCF_VL_A )
    {
        int i, blen = ret/nsmpl, maxlen = 0;
        char *ptr = str.s;
        for (i=0; i<nsmpl; i++)
        {
            char *tmp = ptr;
            int len = 0;
            STR_MOVE_NTH(tmp,tmp,ptr+blen,ialt,len);
            if ( len<0 ) return;   // wrong number of fields: skip
            if ( maxlen < len ) maxlen = len;
            ptr += blen;
        }
        if ( maxlen<blen ) squeeze_format_char(str.s,blen,maxlen,nsmpl);
        bcf_update_format_char(args->hdr,dst,tag,str.s,nsmpl*maxlen);
    }
    else if ( len==BCF_VL_R )
    {
        int i, blen = ret/nsmpl, maxlen = 0;
        char *ptr = str.s;
        for (i=0; i<nsmpl; i++)
        {
            char *tmp = ptr;
            int len = 0;
            STR_MOVE_NTH(ptr,tmp,ptr+blen,0,len);
            ptr[len]=','; tmp++; len++;
            STR_MOVE_NTH(&ptr[len],tmp,ptr+blen,ialt,len);
            if ( len<0 ) return;   // wrong number of fields: skip
            if ( maxlen < len ) maxlen = len;
            ptr += blen;
        }
        if ( maxlen<blen ) squeeze_format_char(str.s,blen,maxlen,nsmpl);
        bcf_update_format_char(args->hdr,dst,tag,str.s,nsmpl*maxlen);
    }
    else if ( len==BCF_VL_G )
    {
        int i, blen = ret/nsmpl, maxlen = 0, i0a = bcf_alleles2gt(0,ialt+1), iaa = bcf_alleles2gt(ialt+1,ialt+1);
        char *ptr = str.s;
        for (i=0; i<nsmpl; i++)
        {
            char *se = ptr, *sx = ptr+blen;
            int nfields = 1;
            while ( *se && se<sx )
            {
                if ( *se==',' ) nfields++;
                se++;
            }
            assert( nfields==src->n_allele*(src->n_allele+1)/2 || nfields==src->n_allele );
            int len = 0;
            if ( nfields==src->n_allele )   // haploid
            {
                char *tmp = ptr;
                STR_MOVE_NTH(&ptr[len],tmp,ptr+blen,0,len);
                ptr[len]=','; tmp++; len++;
                STR_MOVE_NTH(&ptr[len],tmp,ptr+blen,ialt,len);
                if ( len<0 ) return;   // wrong number of fields: skip
            }
            else    // diploid
            {
                char *tmp = ptr;
                STR_MOVE_NTH(&ptr[len],tmp,ptr+blen,0,len);
                ptr[len]=','; tmp++; len++;
                STR_MOVE_NTH(&ptr[len],tmp,ptr+blen,i0a-1,len);
                if ( len<0 ) return;   // wrong number of fields: skip
                ptr[len]=','; tmp++; len++;
                STR_MOVE_NTH(&ptr[len],tmp,ptr+blen,iaa-i0a-1,len);
                if ( len<0 ) return;   // wrong number of fields: skip
            }
            if ( maxlen < len ) maxlen = len;
            ptr += blen;
        }
        if ( maxlen<blen ) squeeze_format_char(str.s,blen,maxlen,nsmpl);
        bcf_update_format_char(args->hdr,dst,tag,str.s,nsmpl*maxlen);
    }
    else
        bcf_update_format_char(args->hdr,dst,tag,str.s,str.l);
}


static void split_multiallelic_to_biallelics(args_t *args, bcf1_t *line)
{
    int i;

    bcf_unpack(line, BCF_UN_ALL);

    // Init the target biallelic lines
    args->ntmp_lines = line->n_allele-1;
    if ( args->mtmp_lines < args->ntmp_lines )
    {
        args->tmp_lines = (bcf1_t **)realloc(args->tmp_lines,sizeof(bcf1_t*)*args->ntmp_lines);
        for (i=args->mtmp_lines; i<args->ntmp_lines; i++)
            args->tmp_lines[i] = bcf_init1();
        args->mtmp_lines = args->ntmp_lines;
    }
    kstring_t tmp = {0,0,0};
    kputs(line->d.allele[0], &tmp);
    kputc(',', &tmp);
    int rlen  = tmp.l;
    int gt_id = bcf_hdr_id2int(args->hdr,BCF_DT_ID,"GT");
    for (i=0; i<args->ntmp_lines; i++)  // for each ALT allele
    {
        bcf1_t *dst = args->tmp_lines[i];
        bcf_clear(dst);

        dst->rid  = line->rid;
        dst->pos  = line->pos;
        dst->qual = line->qual;

        tmp.l = rlen;
        kputs(line->d.allele[i+1],&tmp);
        bcf_update_alleles_str(args->hdr,dst,tmp.s);

        if ( line->d.n_flt ) bcf_update_filter(args->hdr, dst, line->d.flt, line->d.n_flt);

        int j;
        for (j=0; j<line->n_info; j++)
        {
            bcf_info_t *info = &line->d.info[j];
            int type = bcf_hdr_id2type(args->hdr,BCF_HL_INFO,info->key);
            if ( type==BCF_HT_INT || type==BCF_HT_REAL ) split_info_numeric(args, line, info, i, dst);
            else if ( type==BCF_HT_FLAG ) split_info_flag(args, line, info, i, dst);
            else split_info_string(args, line, info, i, dst);
        }

        dst->n_sample = line->n_sample;
        for (j=0; j<line->n_fmt; j++)
        {
            bcf_fmt_t *fmt = &line->d.fmt[j];
            int type = bcf_hdr_id2type(args->hdr,BCF_HL_FMT,fmt->id);
            if ( fmt->id==gt_id ) split_format_genotype(args, line, fmt, i, dst);
            else if ( type==BCF_HT_INT || type==BCF_HT_REAL ) split_format_numeric(args, line, fmt, i, dst);
            else split_format_string(args, line, fmt, i, dst);
        }
    }
    free(tmp.s);
}

// Enlarge FORMAT array containing nsmpl samples each with nals_ori values
// to accommodate nvals values for each sample, filling the gaps with missing
// values. Works also for INFO arrays, with nsmpl set to 1.
#define ENLARGE_ARRAY(type_t,set_missing,arr,narr_bytes,nsmpl,nvals_ori,nvals) \
{ \
    int nbytes_new = (nsmpl)*(nvals)*sizeof(type_t); \
    hts_expand(uint8_t,nbytes_new,narr_bytes,arr); \
    int ismpl, k; \
    for (ismpl=nsmpl-1; ismpl>=0; ismpl--) \
    { \
        type_t *dst_ptr = ((type_t*)arr) + ismpl*(nvals); \
        type_t *src_ptr = ((type_t*)arr) + ismpl*nvals_ori; \
        memmove(dst_ptr,src_ptr,sizeof(type_t)*nvals_ori); \
        for (k=nvals_ori; k<nvals; k++) set_missing; \
    } \
}
static void merge_info_numeric(args_t *args, bcf1_t **lines, int nlines, bcf_info_t *info, bcf1_t *dst)
{
    #define BRANCH_NUMERIC(type,type_t,set_missing,is_vector_end) \
    { \
        const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key); \
        int ntmp = args->ntmp_arr1 / sizeof(type_t); \
        int nvals_ori = bcf_get_info_##type(args->hdr,lines[0],tag,&args->tmp_arr1,&ntmp); \
        args->ntmp_arr1 = ntmp * sizeof(type_t); \
        assert( nvals_ori>0 ); \
        type_t *vals = (type_t*) args->tmp_arr1, *vals2; \
        int i,k,len = bcf_hdr_id2length(args->hdr,BCF_HL_INFO,info->key);  \
        if ( len==BCF_VL_A ) \
        { \
            assert( nvals_ori==lines[0]->n_allele - 1); \
            int nvals = dst->n_allele - 1; \
            ENLARGE_ARRAY(type_t,set_missing,args->tmp_arr1,args->ntmp_arr1,1,nvals_ori,nvals); \
            vals = (type_t*) args->tmp_arr1; \
            for (i=1; i<nlines; i++) \
            { \
                int ntmp2 = args->ntmp_arr2 / sizeof(type_t); \
                int nvals2 = bcf_get_info_##type(args->hdr,lines[i],tag,&args->tmp_arr2,&ntmp2); \
                args->ntmp_arr2 = ntmp2 * sizeof(type_t); \
                assert( nvals2==lines[i]->n_allele-1 ); \
                vals2 = (type_t*) args->tmp_arr2; \
                for (k=0; k<nvals2; k++) \
                { \
                    if ( is_vector_end ) break; \
                    vals[ args->maps[i].map[k+1] - 1 ] = vals2[k]; \
                } \
            } \
            bcf_update_info_##type(args->hdr,dst,tag,args->tmp_arr1,nvals); \
        } \
        else if ( len==BCF_VL_R ) \
        { \
            assert( nvals_ori==lines[0]->n_allele ); \
            int nvals = dst->n_allele; \
            ENLARGE_ARRAY(type_t,set_missing,args->tmp_arr1,args->ntmp_arr1,1,nvals_ori,nvals); \
            vals = (type_t*) args->tmp_arr1; \
            for (i=1; i<nlines; i++) \
            { \
                int ntmp2 = args->ntmp_arr2 / sizeof(type_t); \
                int nvals2 = bcf_get_info_##type(args->hdr,lines[i],tag,&args->tmp_arr2,&ntmp2); \
                args->ntmp_arr2 = ntmp2 * sizeof(type_t); \
                assert( nvals2==lines[i]->n_allele ); \
                vals2 = (type_t*) args->tmp_arr2; \
                for (k=0; k<nvals2; k++) \
                { \
                    if ( is_vector_end ) break; \
                    vals[ args->maps[i].map[k] ] = vals2[k]; \
                } \
            } \
            bcf_update_info_##type(args->hdr,dst,tag,args->tmp_arr1,nvals); \
        } \
        else if ( len==BCF_VL_G ) \
        { \
            assert( nvals_ori==lines[0]->n_allele*(lines[0]->n_allele+1)/2 );   /* expecting diploid gt in INFO */ \
            int nvals = dst->n_allele*(dst->n_allele+1)/2; \
            ENLARGE_ARRAY(type_t,set_missing,args->tmp_arr1,args->ntmp_arr1,1,nvals_ori,nvals); \
            vals = (type_t*) args->tmp_arr1; \
            for (i=1; i<nlines; i++) \
            { \
                int ntmp2 = args->ntmp_arr2 / sizeof(type_t); \
                int nvals2 = bcf_get_info_##type(args->hdr,lines[i],tag,&args->tmp_arr2,&ntmp2); \
                args->ntmp_arr2 = ntmp2 * sizeof(type_t); \
                assert( nvals2==lines[i]->n_allele*(lines[i]->n_allele+1)/2 ); \
                vals2 = (type_t*) args->tmp_arr2; \
                int ia,ib; \
                k = 0; \
                for (ia=0; ia<lines[i]->n_allele; ia++) \
                { \
                    for (ib=0; ib<=ia; ib++) \
                    { \
                        if ( is_vector_end ) break; \
                        int l = bcf_alleles2gt(args->maps[i].map[ia],args->maps[i].map[ib]); \
                        vals[l] = vals2[k];  \
                        k++; \
                    } \
                } \
            } \
            bcf_update_info_##type(args->hdr,dst,tag,args->tmp_arr1,nvals); \
        } \
        else \
            bcf_update_info_##type(args->hdr,dst,tag,vals,nvals_ori); \
    }
    switch (bcf_hdr_id2type(args->hdr,BCF_HL_INFO,info->key))
    {
        case BCF_HT_INT:  BRANCH_NUMERIC(int32, int32_t, dst_ptr[k]=bcf_int32_missing, vals2[k]==bcf_int32_vector_end); break;
        case BCF_HT_REAL: BRANCH_NUMERIC(float, float, bcf_float_set_missing(dst_ptr[k]), bcf_float_is_vector_end(vals2[k])); break;
    }
    #undef BRANCH_NUMERIC
}
static void merge_info_flag(args_t *args, bcf1_t **lines, int nlines, bcf_info_t *info, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key);
    int ret = bcf_get_info_flag(args->hdr,lines[0],tag,&args->tmp_arr1,&args->ntmp_arr1);
    bcf_update_info_flag(args->hdr,dst,tag,NULL,ret);
}
int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst); // see vcfmerge.c
static void merge_info_string(args_t *args, bcf1_t **lines, int nlines, bcf_info_t *info, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key);

    kstring_t str;
    str.m = args->ntmp_arr1;
    str.l = 0;
    str.s = (char*) args->tmp_arr1;

    int i, j, len = bcf_hdr_id2length(args->hdr,BCF_HL_INFO,info->key);
    if ( len==BCF_VL_A || len==BCF_VL_R )
    {
        int jfrom = len==BCF_VL_A ? 1 : 0;
        kputc('.',&str);
        for (i=jfrom+1; i<dst->n_allele; i++) kputs(",.",&str);
        for (i=0; i<nlines; i++)
        {
            bcf_info_t *src = bcf_get_info(args->hdr,lines[i],tag);
            for (j=jfrom; j<lines[i]->n_allele; j++)
                copy_string_field((char*)src->vptr, j-jfrom, src->len, &str, args->maps[i].map[j]-jfrom);
        }
        str.s[str.l] = 0;
        args->tmp_arr1  = (uint8_t*) str.s;
        args->ntmp_arr1 = str.m;
        bcf_update_info_string(args->hdr,dst,tag,str.s);
    }
    else if ( len==BCF_VL_G )
    {
        int ngts = dst->n_allele*(dst->n_allele+1)/2;
        kputc('.',&str);
        for (i=1; i<ngts; i++) kputs(",.",&str);
        for (i=0; i<nlines; i++)
        {
            bcf_info_t *src = bcf_get_info(args->hdr,lines[i],tag);
            int iori, jori, kori = 0;
            for (iori=0; iori<lines[i]->n_allele; iori++)
            {
                int inew = args->maps[i].map[iori];
                for (jori=0; jori<=iori; jori++)
                {
                    int jnew = args->maps[i].map[jori];
                    int knew = bcf_alleles2gt(inew,jnew);
                    copy_string_field((char*)src->vptr,kori,src->len,&str,knew);
                    kori++;
                }
            }
        }
        str.s[str.l] = 0;
        args->tmp_arr1  = (uint8_t*) str.s;
        args->ntmp_arr1 = str.m;
        bcf_update_info_string(args->hdr,dst,tag,str.s);
    }
    else
    {
        bcf_get_info_string(args->hdr,lines[0],tag,&args->tmp_arr1,&args->ntmp_arr1);
        bcf_update_info_string(args->hdr,dst,tag,args->tmp_arr1);
    }

}
static void merge_format_genotype(args_t *args, bcf1_t **lines, int nlines, bcf_fmt_t *fmt, bcf1_t *dst)
{
    int ntmp = args->ntmp_arr1 / 4;
    int ngts = bcf_get_genotypes(args->hdr,lines[0],&args->tmp_arr1,&ntmp);
    args->ntmp_arr1 = ntmp * 4;
    assert( ngts >0 );

    int nsmpl = bcf_hdr_nsamples(args->hdr);
    ngts /= nsmpl;

    int i, j, k;
    for (i=1; i<nlines; i++)
    {
        int ntmp2 = args->ntmp_arr2 / 4;
        int ngts2 = bcf_get_genotypes(args->hdr,lines[i],&args->tmp_arr2,&ntmp2);
        args->ntmp_arr2 = ntmp2 * 4;
        ngts2 /= nsmpl;
        assert( ngts==ngts2 );

        int32_t *gt  = (int32_t*) args->tmp_arr1;
        int32_t *gt2 = (int32_t*) args->tmp_arr2;
        for (j=0; j<nsmpl; j++)
        {
            for (k=0; k<ngts; k++)
            {
                if ( gt2[k]==bcf_int32_vector_end ) break;
                if ( bcf_gt_is_missing(gt2[k]) || bcf_gt_allele(gt2[k])==0 ) continue;
                if ( gt2[k]==0 ) gt[k] = 0; // missing genotype
                else
                {
                    int ial = bcf_gt_allele(gt2[k]);
                    assert( ial<args->maps[i].nals );
                    gt[k] = bcf_gt_unphased( args->maps[i].map[ial] ) | bcf_gt_is_phased(gt[k]);
                }
            }
            gt  += ngts;
            gt2 += ngts;
        }
    }
    bcf_update_genotypes(args->hdr,dst,args->tmp_arr1,ngts*nsmpl);
}
static void merge_format_numeric(args_t *args, bcf1_t **lines, int nlines, bcf_fmt_t *fmt, bcf1_t *dst)
{
    #define BRANCH_NUMERIC(type,type_t,set_missing,is_vector_end) \
    { \
        const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id); \
        int ntmp = args->ntmp_arr1 / sizeof(type_t); \
        int nvals_ori = bcf_get_format_##type(args->hdr,lines[0],tag,&args->tmp_arr1,&ntmp); \
        args->ntmp_arr1 = ntmp * sizeof(type_t); \
        assert( nvals_ori>0 ); \
        type_t *vals2, *vals = (type_t *) args->tmp_arr1; \
        int len = bcf_hdr_id2length(args->hdr,BCF_HL_FMT,fmt->id); \
        int i, j, k, nsmpl = bcf_hdr_nsamples(args->hdr); \
        nvals_ori /= nsmpl; \
        if ( len==BCF_VL_A ) \
        { \
            int nvals = dst->n_allele - 1; \
            ENLARGE_ARRAY(type_t,set_missing,args->tmp_arr1,args->ntmp_arr1,nsmpl,nvals_ori,nvals); \
            for (i=1; i<nlines; i++) \
            { \
                int ntmp2 = args->ntmp_arr2 / sizeof(type_t); \
                int nvals2 = bcf_get_format_##type(args->hdr,lines[i],tag,&args->tmp_arr2,&ntmp2); \
                args->ntmp_arr2 = ntmp2 * sizeof(type_t); \
                nvals2 /= nsmpl; \
                assert( nvals2==lines[i]->n_allele-1 ); \
                vals  = (type_t*) args->tmp_arr1; \
                vals2 = (type_t*) args->tmp_arr2; \
                for (j=0; j<nsmpl; j++) \
                { \
                    for (k=0; k<nvals2; k++) \
                    { \
                        if ( is_vector_end ) break; \
                        vals[ args->maps[i].map[k+1] - 1 ] = vals2[k]; \
                    } \
                    vals  += nvals; \
                    vals2 += nvals2; \
                } \
            } \
            bcf_update_format_##type(args->hdr,dst,tag,args->tmp_arr1,nvals*nsmpl); \
        } \
        else if ( len==BCF_VL_R ) \
        { \
            int nvals = dst->n_allele; \
            ENLARGE_ARRAY(type_t,set_missing,args->tmp_arr1,args->ntmp_arr1,nsmpl,nvals_ori,nvals); \
            for (i=1; i<nlines; i++) \
            { \
                int ntmp2 = args->ntmp_arr2 / sizeof(type_t); \
                int nvals2 = bcf_get_format_##type(args->hdr,lines[i],tag,&args->tmp_arr2,&ntmp2); \
                args->ntmp_arr2 = ntmp2 * sizeof(type_t); \
                nvals2 /= nsmpl; \
                assert( nvals2==lines[i]->n_allele ); \
                vals  = (type_t*) args->tmp_arr1; \
                vals2 = (type_t*) args->tmp_arr2; \
                for (j=0; j<nsmpl; j++) \
                { \
                    for (k=0; k<nvals2; k++) \
                    { \
                        if ( is_vector_end ) break; \
                        vals[ args->maps[i].map[k] ] = vals2[k]; \
                    } \
                    vals  += nvals; \
                    vals2 += nvals2; \
                } \
            } \
            bcf_update_format_##type(args->hdr,dst,tag,args->tmp_arr1,nvals*nsmpl); \
        } \
        else if ( len==BCF_VL_G ) \
        { \
            int all_haploid = nvals_ori==lines[0]->n_allele ? 1 : 0; \
            int nvals = all_haploid ? dst->n_allele : dst->n_allele*(dst->n_allele+1)/2; \
            ENLARGE_ARRAY(type_t,set_missing,args->tmp_arr1,args->ntmp_arr1,nsmpl,nvals_ori,nvals); \
            for (i=1; i<nlines; i++) \
            { \
                int ntmp2 = args->ntmp_arr2 / sizeof(type_t); \
                int nvals2 = bcf_get_format_##type(args->hdr,lines[i],tag,&args->tmp_arr2,&ntmp2); \
                args->ntmp_arr2 = ntmp2 * sizeof(type_t); \
                nvals2 /= nsmpl; \
                assert( nvals2==lines[i]->n_allele || nvals2==lines[i]->n_allele*(lines[i]->n_allele+1)/2); \
                vals  = (type_t*) args->tmp_arr1; \
                vals2 = (type_t*) args->tmp_arr2; \
                for (j=0; j<nsmpl; j++) \
                { \
                    int haploid = all_haploid; \
                    if ( !haploid ) \
                    { \
                        for (k=0; k<nvals2; k++) if ( is_vector_end ) break; \
                        if ( k!=nvals2 ) haploid = 1; \
                    } \
                    if ( haploid ) \
                    { \
                        for (k=0; k<nvals2; k++) \
                        { \
                            if ( is_vector_end ) break; \
                            vals[ args->maps[i].map[k] ] = vals2[k]; \
                        } \
                    } \
                    else \
                    { \
                        k = 0; \
                        int ia,ib; \
                        for (ia=0; ia<lines[i]->n_allele; ia++) \
                        { \
                            for (ib=0; ib<=ia; ib++) \
                            { \
                                int l = bcf_alleles2gt(args->maps[i].map[ia],args->maps[i].map[ib]); \
                                vals[l] = vals2[k]; \
                                k++; \
                            } \
                        } \
                    } \
                    vals  += nvals; \
                    vals2 += nvals2; \
                } \
            } \
            bcf_update_format_##type(args->hdr,dst,tag,args->tmp_arr1,nvals*nsmpl); \
        } \
        else \
            bcf_update_format_##type(args->hdr,dst,tag,args->tmp_arr1,nvals_ori*nsmpl); \
    }
    switch (bcf_hdr_id2type(args->hdr,BCF_HL_FMT,fmt->id))
    {
        case BCF_HT_INT:  BRANCH_NUMERIC(int32, int32_t, dst_ptr[k]=bcf_int32_missing, vals2[k]==bcf_int32_vector_end); break;
        case BCF_HT_REAL: BRANCH_NUMERIC(float, float, bcf_float_set_missing(dst_ptr[k]), bcf_float_is_vector_end(vals2[k])); break;
    }
    #undef BRANCH_NUMERIC
}
static void merge_format_string(args_t *args, bcf1_t **lines, int nlines, bcf_fmt_t *fmt, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id);

    int i, j, k, len = bcf_hdr_id2length(args->hdr,BCF_HL_FMT,fmt->id);
    if ( len!=BCF_VL_A && len!=BCF_VL_R && len!=BCF_VL_G )
    {
        int nret = bcf_get_format_char(args->hdr,lines[0],tag,&args->tmp_arr1,&args->ntmp_arr1);
        bcf_update_format_char(args->hdr,dst,tag,args->tmp_arr1,nret);
        return;
    }

    int nsmpl = bcf_hdr_nsamples(args->hdr);
    for (i=0; i<nsmpl; i++) args->tmp_str[i].l = 0;

    if ( len==BCF_VL_A || len==BCF_VL_R )
    {
        int jfrom = len==BCF_VL_A ? 1 : 0;
        for (i=0; i<nsmpl; i++)
        {
            kstring_t *tmp = &args->tmp_str[i];
            kputc('.',tmp);
            for (k=jfrom+1; k<dst->n_allele; k++) kputs(",.",tmp);
        }
        for (i=0; i<nlines; i++)
        {
            int nret = bcf_get_format_char(args->hdr,lines[i],tag,&args->tmp_arr1,&args->ntmp_arr1);
            nret /= nsmpl;
            for (k=0; k<nsmpl; k++)
            {
                kstring_t *tmp = &args->tmp_str[k];
                char *src = (char*)args->tmp_arr1 + k*nret;
                for (j=jfrom; j<lines[i]->n_allele; j++)
                    copy_string_field(src, j-jfrom, nret, tmp, args->maps[i].map[j]-jfrom);
            }
        }
    }
    else if ( len==BCF_VL_G )
    {
        hts_expand(uint8_t,nsmpl,args->ntmp_arr2,args->tmp_arr2);
        uint8_t *haploid = args->tmp_arr2;
        int nret = bcf_get_format_char(args->hdr,lines[0],tag,&args->tmp_arr1,&args->ntmp_arr1);
        nret /= nsmpl;
        for (i=0; i<nsmpl; i++)
        {
            char *ss = (char*)args->tmp_arr1 + i*nret, *se = ss+nret;
            int nfields = 1;
            while ( *ss && ss<se )
            {
                if ( *ss==',' ) nfields++;
                ss++;
            }
            if ( nfields==lines[0]->n_allele )
            {
                haploid[i] = 1;
                nfields = dst->n_allele;
            }
            else if ( nfields==lines[0]->n_allele*(lines[0]->n_allele+1)/2 )
            {
                haploid[i] = 0;
                nfields = dst->n_allele*(dst->n_allele+1)/2;
            }
            else error("The field %s at %s:%d neither diploid nor haploid?\n", tag,bcf_seqname(args->hdr,dst),dst->pos+1);

            kstring_t *tmp = &args->tmp_str[i];
            kputc('.',tmp);
            for (j=1; j<nfields; j++) kputs(",.",tmp);
        }
        for (i=0; i<nlines; i++)
        {
            if ( i ) // we already have a copy
            {
                nret = bcf_get_format_char(args->hdr,lines[i],tag,&args->tmp_arr1,&args->ntmp_arr1);
                nret /= nsmpl;
            }
            for (k=0; k<nsmpl; k++)
            {
                kstring_t *tmp = &args->tmp_str[k];
                char *src = (char*)args->tmp_arr1 + k*nret;
                if ( haploid[k] )
                {
                    for (j=0; j<lines[i]->n_allele; j++)
                        copy_string_field(src,j,nret, tmp, args->maps[i].map[j]);
                }
                else
                {
                    int iori, jori, kori = 0;
                    for (iori=0; iori<lines[i]->n_allele; iori++)
                    {
                        int inew = args->maps[i].map[iori];
                        for (jori=0; jori<=iori; jori++)
                        {
                            int jnew = args->maps[i].map[jori];
                            int knew = bcf_alleles2gt(inew,jnew);
                            copy_string_field(src,kori,nret,tmp,knew);
                            kori++;
                        }
                    }
                }
            }
        }
    }
    kstring_t str;
    str.m = args->ntmp_arr2;
    str.l = 0;
    str.s = (char*) args->tmp_arr2;

    int max_len = 0;
    for (i=0; i<nsmpl; i++)
        if ( max_len < args->tmp_str[i].l ) max_len = args->tmp_str[i].l;
    for (i=0; i<nsmpl; i++)
    {
        kstring_t *tmp = &args->tmp_str[i];
        kputsn(tmp->s,tmp->l,&str);
        for (j=tmp->l; j<max_len; j++) kputc(0,tmp);
    }
    args->ntmp_arr2 = str.m;
    args->tmp_arr2  = (uint8_t*)str.s;
    bcf_update_format_char(args->hdr,dst,tag,str.s,str.l);
}

char **merge_alleles(char **a, int na, int *map, char **b, int *nb, int *mb);   // see vcfmerge.c
static void merge_biallelics_to_multiallelic(args_t *args, bcf1_t *dst, bcf1_t **lines, int nlines)
{
    int i;
    for (i=0; i<nlines; i++)
        bcf_unpack(lines[i], BCF_UN_ALL);

    dst->rid  = lines[0]->rid;
    dst->pos  = lines[0]->pos;

    // take max for QUAL
    bcf_float_set_missing(dst->qual);
    for (i=0; i<nlines; i++) {
        if (bcf_float_is_missing(lines[i]->qual)) continue;
        if (bcf_float_is_missing(dst->qual) || dst->qual<lines[i]->qual)
            dst->qual = lines[i]->qual;
    }

    bcf_update_id(args->hdr, dst, lines[0]->d.id);

    // Merge and set the alleles, create a mapping from source allele indexes to dst idxs
    hts_expand0(map_t,nlines,args->mmaps,args->maps);   // a mapping for each line
    args->nals = args->maps[0].nals = lines[0]->n_allele;
    hts_expand(int,args->maps[0].nals,args->maps[0].mals,args->maps[0].map);
    hts_expand(char*,args->nals,args->mals,args->als);
    for (i=0; i<args->maps[0].nals; i++)
    {
        args->maps[0].map[i] = i;
        args->als[i] = strdup(lines[0]->d.allele[i]);
    }
    for (i=1; i<nlines; i++)
    {
        if (lines[i]->d.id[0]!='.' || lines[i]->d.id[1]) {
            kstring_t tmp = {0,0,0};
            if (dst->d.id[0]=='.' && !dst->d.id[1])
                kputs(lines[i]->d.id, &tmp);
            else
                ksprintf(&tmp, "%s;%s", dst->d.id, lines[i]->d.id);
            bcf_update_id(args->hdr, dst, tmp.s);
            free(tmp.s);
        }
        args->maps[i].nals = lines[i]->n_allele;
        hts_expand(int,args->maps[i].nals,args->maps[i].mals,args->maps[i].map);
        args->als = merge_alleles(lines[i]->d.allele, lines[i]->n_allele, args->maps[i].map, args->als, &args->nals, &args->mals);
        if ( !args->als ) error("Failed to merge alleles at %s:%d\n", bcf_seqname(args->hdr,dst),dst->pos+1);
    }
    bcf_update_alleles(args->hdr, dst, (const char**)args->als, args->nals);
    for (i=0; i<args->nals; i++) free(args->als[i]);

    if ( lines[0]->d.n_flt ) bcf_update_filter(args->hdr, dst, lines[0]->d.flt, lines[0]->d.n_flt);
    for (i=1; i<nlines; i++) {
        int j;
        for (j=0; j<lines[i]->d.n_flt; j++) {
            // if strict_filter, set FILTER to PASS if any site PASS
            // otherwise accumulate FILTERs
            if (lines[i]->d.flt[j] == bcf_hdr_id2int(args->hdr, BCF_DT_ID, "PASS")) {
                if (args->strict_filter) {
                    bcf_update_filter(args->hdr, dst, lines[i]->d.flt, lines[i]->d.n_flt);
                    break;
                }
                else
                    continue;
            }
            bcf_add_filter(args->hdr, dst, lines[i]->d.flt[j]);
        }
    }

    // merge info
    for (i=0; i<lines[0]->n_info; i++)
    {
        bcf_info_t *info = &lines[0]->d.info[i];
        int type = bcf_hdr_id2type(args->hdr,BCF_HL_INFO,info->key);
        if ( type==BCF_HT_INT || type==BCF_HT_REAL ) merge_info_numeric(args, lines, nlines, info, dst);
        else if ( type==BCF_HT_FLAG ) merge_info_flag(args, lines, nlines, info, dst);
        else merge_info_string(args, lines, nlines, info, dst);
    }

    // merge format
    int gt_id = bcf_hdr_id2int(args->hdr,BCF_DT_ID,"GT");
    dst->n_sample = lines[0]->n_sample;
    for (i=0; i<lines[0]->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &lines[0]->d.fmt[i];
        int type = bcf_hdr_id2type(args->hdr,BCF_HL_FMT,fmt->id);
        if ( fmt->id==gt_id ) merge_format_genotype(args, lines, nlines, fmt, dst);
        else if ( type==BCF_HT_INT || type==BCF_HT_REAL ) merge_format_numeric(args, lines, nlines, fmt, dst);
        else merge_format_string(args, lines, nlines, fmt, dst);
    }
}

#define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }
static void mrows_schedule(args_t *args, bcf1_t **line)
{
    int i,m;
    if ( args->mrows_collapse==COLLAPSE_ANY ||  bcf_get_variant_types(*line)&COLLAPSE_SNPS )
    {
        args->nalines++;
        m = args->malines;
        hts_expand(bcf1_t*,args->nalines,args->malines,args->alines);
        for (i=m; i<args->malines; i++) args->alines[i] = bcf_init1();
        SWAP(bcf1_t*, args->alines[args->nalines-1], *line);
    }
    else
    {
        args->nblines++;
        m = args->mblines;
        hts_expand(bcf1_t*,args->nblines,args->mblines,args->blines);
        for (i=m; i<args->mblines; i++) args->blines[i] = bcf_init1();
        SWAP(bcf1_t*, args->blines[args->nblines-1], *line);
    }
}
static int mrows_ready_to_flush(args_t *args, bcf1_t *line)
{
    if ( args->nalines && (args->alines[0]->rid!=line->rid || args->alines[0]->pos!=line->pos) ) return 1;
    if ( args->nblines && (args->blines[0]->rid!=line->rid || args->blines[0]->pos!=line->pos) ) return 1;
    return 0;
}
static bcf1_t *mrows_flush(args_t *args)
{
    if ( args->nalines )
    {
        if ( args->nalines==1 )
        {
            args->nalines = 0;
            return args->alines[0];
        }
        bcf_clear(args->mrow_out);
        merge_biallelics_to_multiallelic(args, args->mrow_out, args->alines, args->nalines);
        args->nalines = 0;
        return args->mrow_out;
    }
    else if ( args->nblines )
    {
        if ( args->nblines==1 )
        {
            args->nblines = 0;
            return args->blines[0];
        }
        bcf_clear(args->mrow_out);
        merge_biallelics_to_multiallelic(args, args->mrow_out, args->blines, args->nblines);
        args->nblines = 0;
        return args->mrow_out;
    }
    return NULL;
}
static void flush_buffer(args_t *args, htsFile *file, int n)
{
    bcf1_t *line;
    int i, k, prev_rid = -1, prev_pos = 0, prev_type = 0;
    for (i=0; i<n; i++)
    {
        k = rbuf_shift(&args->rbuf);
        if ( args->mrows_op==MROWS_SPLIT )
        {
            int split = 1;
            if ( args->mrows_collapse!=COLLAPSE_BOTH && args->mrows_collapse!=COLLAPSE_ANY )
            {
                if ( !(bcf_get_variant_types(args->lines[k]) & args->mrows_collapse) ) split = 0;
            }
            if ( split && args->lines[k]->n_allele>2 )
            {
                split_multiallelic_to_biallelics(args, args->lines[k]);
                int j;
                for (j=0; j<args->ntmp_lines; j++)
                    bcf_write1(file, args->hdr, args->tmp_lines[j]);
                continue;
            }
        }
        if ( args->mrows_op==MROWS_MERGE )
        {
            if ( mrows_ready_to_flush(args, args->lines[k]) )
            {
                while ( (line=mrows_flush(args)) ) bcf_write1(file, args->hdr, line);
            }
            int merge = 1;
            if ( args->mrows_collapse!=COLLAPSE_BOTH && args->mrows_collapse!=COLLAPSE_ANY )
            {
                if ( !(bcf_get_variant_types(args->lines[k]) & args->mrows_collapse) ) merge = 0;
            }
            if ( merge )
            {
                mrows_schedule(args, &args->lines[k]);
                continue;
            }
        }
        // todo: merge with next record if POS and the type are same. For now, just discard if asked to do so.
        if ( args->rmdup )
        {
            int line_type = bcf_get_variant_types(args->lines[k]);
            if ( prev_rid>=0 && prev_rid==args->lines[k]->rid && prev_pos==args->lines[k]->pos && prev_type==line_type )
                continue;
            prev_rid  = args->lines[k]->rid;
            prev_pos  = args->lines[k]->pos;
            prev_type = line_type;
        }
        bcf_write1(file, args->hdr, args->lines[k]);
    }
    if ( args->mrows_op==MROWS_MERGE && !args->rbuf.n )
    {
        while ( (line=mrows_flush(args)) ) bcf_write1(file, args->hdr, line);
    }
}

static void init_data(args_t *args)
{
    args->hdr = args->files->readers[0].header;
    rbuf_init(&args->rbuf, 100);
    args->lines = (bcf1_t**) calloc(args->rbuf.m, sizeof(bcf1_t*));
    if ( args->ref_fname )
    {
        args->fai = fai_load(args->ref_fname);
        if ( !args->fai ) error("Failed to load the fai index: %s\n", args->ref_fname);
    }
    if ( args->mrows_op==MROWS_MERGE )
    {
        args->mrow_out = bcf_init1();
        args->tmp_str = (kstring_t*) calloc(bcf_hdr_nsamples(args->hdr),sizeof(kstring_t));
    }
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->rbuf.m; i++)
        if ( args->lines[i] ) bcf_destroy1(args->lines[i]);
    free(args->lines);
    for (i=0; i<args->mtmp_lines; i++)
        bcf_destroy1(args->tmp_lines[i]);
    free(args->tmp_lines);
    for (i=0; i<args->nalines; i++)
        bcf_destroy1(args->alines[i]);
    free(args->alines);
    for (i=0; i<args->nblines; i++)
        bcf_destroy1(args->blines[i]);
    free(args->blines);
    for (i=0; i<args->mmaps; i++)
        free(args->maps[i].map);
    for (i=0; i<args->ntmp_als; i++)
        free(args->tmp_als[i].s);
    free(args->tmp_als);
    free(args->tmp_als_str.s);
    if ( args->tmp_str )
    {
        for (i=0; i<bcf_hdr_nsamples(args->hdr); i++) free(args->tmp_str[i].s);
        free(args->tmp_str);
    }
    free(args->maps);
    free(args->als);
    free(args->tmp_arr1);
    free(args->tmp_arr2);
    if ( args->mrow_out ) bcf_destroy1(args->mrow_out);
    if ( args->fai ) fai_destroy(args->fai);
    if ( args->mseq ) free(args->seq);
    if ( args->aln.nmat ) free(args->aln.mat);
    if ( args->aln.m_arr ) { free(args->aln.ipos_arr); free(args->aln.lref_arr); free(args->aln.lseq_arr); }
}


static void normalize_vcf(args_t *args)
{
    htsFile *out = hts_open(args->output_fname, hts_bcf_wmode(args->output_type));
    if ( out == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    bcf_hdr_append_version(args->hdr, args->argc, args->argv, "bcftools_norm");
    bcf_hdr_write(out, args->hdr);

    while ( bcf_sr_next_line(args->files) )
    {
        args->ntotal++;

        bcf1_t *line = args->files->readers[0].buffer[0];
        if ( args->fai )
        {
            if ( realign(args, line)<0 && args->check_ref & CHECK_REF_SKIP )
            {
                args->nskipped++;
                continue;   // exclude broken VCF lines
            }
        }

        // still on the same chromosome?
        int i, j, ilast = rbuf_last(&args->rbuf);
        if ( ilast>=0 && line->rid != args->lines[ilast]->rid ) flush_buffer(args, out, args->rbuf.n); // new chromosome

        // insert into sorted buffer
        i = j = ilast = rbuf_append(&args->rbuf);
        if ( !args->lines[i] ) args->lines[i] = bcf_init1();
        SWAP(bcf1_t*, args->files->readers[0].buffer[0], args->lines[i]);
        while ( rbuf_prev(&args->rbuf,&i) )
        {
            if ( args->lines[i]->pos > args->lines[j]->pos ) SWAP(bcf1_t*, args->lines[i], args->lines[j]);
            j = i;
        }

        // find out how many sites to flush
        j = 0;
        for (i=-1; rbuf_next(&args->rbuf,&i); )
        {
            if ( args->lines[ilast]->pos - args->lines[i]->pos < args->buf_win ) break;
            j++;
        }
        if ( args->rbuf.n==args->rbuf.m ) j = 1;
        if ( j>0 ) flush_buffer(args, out, j);
    }
    flush_buffer(args, out, args->rbuf.n);
    hts_close(out);

    fprintf(stderr,"Lines total/modified/skipped:\t%d/%d/%d\n", args->ntotal,args->nchanged,args->nskipped);
}

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Left-align and normalize indels; check if REF alleles match the reference;\n");
    fprintf(stderr, "         split multiallelic sites into multiple rows; recover multiallelics from\n");
    fprintf(stderr, "         multiple rows.\n");
    fprintf(stderr, "Usage:   bcftools norm [options] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -c, --check-ref <e|w|x>           check REF alleles and exit (e), warn (w), exclude (x) bad sites [e]\n");
    fprintf(stderr, "    -D, --remove-duplicates           remove duplicate lines of the same type.\n");
    fprintf(stderr, "    -f, --fasta-ref <file>            reference sequence\n");
    fprintf(stderr, "    -m, --multiallelics <-|+>[type]   split multiallelics (-) or join biallelics (+), type: snps|indels|both|any [both]\n");
    fprintf(stderr, "    -o, --output <file>               write output to a file [standard output]\n");
    fprintf(stderr, "    -O, --output-type <type>          'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n");
    fprintf(stderr, "    -r, --regions <region>            restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>         restrict to regions listed in a file\n");
    fprintf(stderr, "    -s, --strict-filter               when merging (-m+), merged site is PASS only if all sites being merged PASS\n");
    fprintf(stderr, "    -t, --targets <region>            similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file <file>         similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "    -w, --site-win <int>              buffer for sorting lines which changed position during realignment [1000]\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfnorm(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->output_fname = "-";
    args->output_type = FT_VCF;
    args->aln_win = 100;
    args->buf_win = 1000;
    args->mrows_collapse = COLLAPSE_BOTH;
    int region_is_file  = 0;
    int targets_is_file = 0;

    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"fasta-ref",1,0,'f'},
        {"multiallelics",1,0,'m'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"site-win",1,0,'W'},
        {"remove-duplicates",0,0,'D'},
        {"output",1,0,'o'},
        {"output-type",1,0,'O'},
        {"check-ref",1,0,'c'},
        {"strict-filter",0,0,'s'},
        {0,0,0,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "hr:R:f:w:Do:O:c:m:t:T:s",loptions,NULL)) >= 0) {
        switch (c) {
            case 'm':
                if ( optarg[0]=='-' ) args->mrows_op = MROWS_SPLIT;
                else if ( optarg[0]=='+' ) args->mrows_op = MROWS_MERGE;
                else error("Expected '+' or '-' with -m\n");
                if ( optarg[1]!=0 )
                {
                    if ( !strcmp("snps",optarg+1) ) args->mrows_collapse = COLLAPSE_SNPS;
                    else if ( !strcmp("indels",optarg+1) ) args->mrows_collapse = COLLAPSE_INDELS;
                    else if ( !strcmp("both",optarg+1) ) args->mrows_collapse = COLLAPSE_BOTH;
                    else if ( !strcmp("any",optarg+1) ) args->mrows_collapse = COLLAPSE_ANY;
                    else error("The argument to -m not recognised: %s\n", optarg);
                }
                break;
            case 'c':
                if ( strchr(optarg,'w') ) args->check_ref |= CHECK_REF_WARN;
                if ( strchr(optarg,'x') ) args->check_ref |= CHECK_REF_SKIP;
                if ( strchr(optarg,'e') ) args->check_ref = CHECK_REF_EXIT; // overrides the above
                break;
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                }
                break;
            case 'o': args->output_fname = optarg; break;
            case 'D': args->rmdup = 1; break;
            case 's': args->strict_filter = 1; break;
            case 'f': args->ref_fname = optarg; break;
            case 'r': args->region = optarg; break;
            case 'R': args->region = optarg; region_is_file = 1; break;
            case 't': args->targets = optarg; break;
            case 'T': args->targets = optarg; targets_is_file = 1; break;
            case 'w':
                args->buf_win = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --site-win %s\n", optarg);
                break;
            case 'h':
            case '?': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( argc>optind+1 ) usage();
    if ( !args->ref_fname && !args->mrows_op && !args->rmdup ) usage();
    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else usage();
    }
    else fname = argv[optind];

    if ( args->region )
    {
        if ( bcf_sr_set_regions(args->files, args->region,region_is_file)<0 )
            error("Failed to read the regions: %s\n", args->region);
    }
    if ( args->targets )
    {
        if ( bcf_sr_set_targets(args->files, args->targets,targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets);
    }

    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->files->errnum));
    if ( args->mrows_op&MROWS_SPLIT && args->rmdup ) error("Cannot combine -D and -m-\n");
    init_data(args);
    normalize_vcf(args);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}

