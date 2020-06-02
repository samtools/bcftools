/*  ext-sort.h -- sort on disk

   Copyright (C) 2020 Genome Research Ltd.

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
#include <unistd.h>     // for unlink()
#include <sys/stat.h>   // for chmod()
#include <assert.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include "bcftools.h"
#include "extsort.h"
#include "kheap.h"

typedef struct
{
    extsort_t *es;  // this is to get access to extsort_cmp_f from kheap
    FILE *fp;
    char *fname;
    void *dat;
}
blk_t;

static inline int blk_is_smaller(blk_t **aptr, blk_t **bptr);
KHEAP_INIT(blk, blk_t*, blk_is_smaller)     /* defines khp_blk_t */

struct _extsort_t
{
    size_t dat_size, mem, max_mem;
    char *tmp_dir;
    extsort_cmp_f cmp;

    size_t nbuf, mbuf, nblk;
    blk_t **blk;
    void **buf, *tmp_dat;
    khp_blk_t *bhp;
};

static inline int blk_is_smaller(blk_t **aptr, blk_t **bptr)
{
    blk_t *a = *aptr;
    blk_t *b = *bptr;
    int ret = a->es->cmp(&a->dat,&b->dat);
    if ( ret < 0 ) return 1;
    return 0;
}

static void _clean_files(extsort_t *es)
{
    int i;
    for (i=0; i<es->nblk; i++)
    {
        blk_t *blk = es->blk[i];
        if ( blk->fname )
        {
            unlink(blk->fname);
            free(blk->fname);
        }
        free(blk->dat);
        free(blk);
    }
    rmdir(es->tmp_dir);
}

static void clean_files_and_throw(extsort_t *es, const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    _clean_files(es);
    exit(-1);
}

size_t parse_mem_string(const char *str);
void mkdir_p(const char *fmt, ...);

static void _init_tmpdir(extsort_t *es, const char *tmp_dir)
{
#ifdef _WIN32
    char tmp_path[MAX_PATH];
    int ret = GetTempPath(MAX_PATH, tmp_path);
    if (!ret || ret > MAX_PATH)
        error("Could not get the path to the temporary folder\n");
    if (strlen(tmp_path) + strlen("/bcftools-sort.XXXXXX") >= MAX_PATH)
        error("Full path to the temporary folder is too long\n");
    strcat(tmp_path, "/bcftools-sort.XXXXXX");
    es->tmp_dir = strdup(tmp_path);
#else
    es->tmp_dir = tmp_dir ? strdup(tmp_dir) : strdup("/tmp/bcftools-sort.XXXXXX");
#endif
    size_t len = strlen(es->tmp_dir);
    if ( !strcmp("XXXXXX",es->tmp_dir+len-6) )
    {
#ifdef _WIN32
        int ret = mkdir(mktemp(es->tmp_dir), 0700);
        if ( ret ) error("mkdir(%s) failed: %s\n", es->tmp_dir,strerror(errno));
#else
        char *tmp = mkdtemp(es->tmp_dir);
        if ( !tmp ) error("mkdtemp(%s) failed: %s\n",  es->tmp_dir,strerror(errno));
        int ret = chmod(tmp, S_IRUSR|S_IWUSR|S_IXUSR);
        if ( ret ) error("chmod(%s,S_IRUSR|S_IWUSR|S_IXUSR) failed: %s\n", es->tmp_dir,strerror(errno));
#endif
    }
    else
        mkdir_p("%s/",es->tmp_dir);

    fprintf(stderr,"Writing to %s\n",es->tmp_dir);
}

void extsort_set(extsort_t *es, extsort_opt_t key, void *value)
{
    if ( key==DAT_SIZE ) { es->dat_size = *((size_t*)value); return; }
    if ( key==MAX_MEM )
    {
        es->max_mem = parse_mem_string(*((const char**)value));
        if ( es->max_mem <=0 ) clean_files_and_throw(es,"Could not parse the memory string, expected positive number: %s\n",*((const char**)value));
        return;
    }
    if ( key==TMP_DIR ) { _init_tmpdir(es, *((const char**)value)); return; }
    if ( key==FUNC_CMP ) { es->cmp = *((extsort_cmp_f*)value); return; }
}

extsort_t *extsort_alloc(void)
{
    extsort_t *es = (extsort_t*) calloc(1,sizeof(*es));
    es->max_mem = 100e6;
    return es;
}
void extsort_init(extsort_t *es)
{
    assert( es->cmp );
    assert( es->dat_size );
    if ( !es->tmp_dir ) _init_tmpdir(es, NULL);
    es->tmp_dat = malloc(es->dat_size);
}

void extsort_destroy(extsort_t *es)
{
    _clean_files(es);
    free(es->tmp_dat);
    free(es->tmp_dir);
    free(es->blk);
    khp_destroy(blk, es->bhp);
    free(es);
}

static void _buf_flush(extsort_t *es)
{
    if ( !es->nbuf ) return;

    qsort(es->buf, es->nbuf, sizeof(void*), es->cmp);

    es->nblk++;
    es->blk = (blk_t**) realloc(es->blk, sizeof(blk_t*)*es->nblk);
    es->blk[es->nblk-1] = (blk_t*) calloc(1,sizeof(blk_t));
    blk_t *blk = es->blk[es->nblk-1];
    blk->dat = malloc(es->dat_size);
    kstring_t str = {0,0,0};
    ksprintf(&str, "%s/%05d.tmp", es->tmp_dir, (int)es->nblk);
    blk->fname = str.s;
    blk->es    = es;

    FILE *fp = fopen(blk->fname, "w");
    if ( fp == NULL ) clean_files_and_throw(es, "Cannot write %s: %s\n", blk->fname, strerror(errno));
    int i;
    for (i=0; i<es->nbuf; i++)
    {
        if ( fwrite(es->buf[i], es->dat_size, 1, fp)!=1 ) clean_files_and_throw(es, "[%s] Error: cannot write to %s\n", __func__,blk->fname);
        free(es->buf[i]);
    }
    if ( fclose(fp)!=0 ) clean_files_and_throw(es, "[%s] Error: close failed .. %s\n", __func__,blk->fname);

    es->nbuf = 0;
    es->mem  = 0;
}

void extsort_push(extsort_t *es, void *dat)
{
    int delta = sizeof(void*) + es->dat_size;
    if ( es->nbuf && es->mem + delta > es->max_mem ) _buf_flush(es);
    es->nbuf++;
    es->mem += delta;
    hts_expand(void*, es->nbuf, es->mbuf, es->buf);
    es->buf[es->nbuf-1] = dat;
}

// return number of elements read
static int _blk_read(extsort_t *es, blk_t *blk)
{
    int ret = 0;
    if ( !blk->fp ) return ret;
    ret = fread(blk->dat, es->dat_size, 1, blk->fp);
    if ( ret < -1 ) clean_files_and_throw(es, "Error reading %s\n", blk->fname);
    if ( ret == -1 )
    {
        if ( fclose(blk->fp)!=0 ) clean_files_and_throw(es, "Close failed: %s\n", blk->fname);
        blk->fp = 0;
        return ret;
    }
    if ( !ret && fclose(blk->fp)!=0 ) clean_files_and_throw(es, "[%s] Error: close failed .. %s\n", __func__,blk->fname);
    return ret;
}

void extsort_sort(extsort_t *es)
{
    _buf_flush(es);
    free(es->buf);
    es->buf = NULL;

    // fprintf(stderr,"[%s] Merging %d temporary files\n", __func__,(int)es->nblk);

    es->bhp = khp_init(blk);

    // open all blocks, read one record from each, create a heap
    int i;
    for (i=0; i<es->nblk; i++)
    {
        blk_t *blk = es->blk[i];
        blk->fp = fopen(blk->fname, "r");
        if ( !blk->fp ) clean_files_and_throw(es, "Could not read %s: %s\n", blk->fname, strerror(errno));
        int ret = _blk_read(es, blk);
        if ( ret==1 ) khp_insert(blk, es->bhp, &blk);
    }
}

void *extsort_shift(extsort_t *es)
{
    if ( !es->bhp->ndat ) return NULL;
    blk_t *blk = es->bhp->dat[0];

    // swap the pointer which keeps the location of user data so that it is not overwritten by the next read
    void *tmp = es->tmp_dat; es->tmp_dat = blk->dat; blk->dat = tmp;
    khp_delete(blk, es->bhp);

    int ret = _blk_read(es, blk);
    if ( ret==1 ) khp_insert(blk, es->bhp, &blk);

    return es->tmp_dat;
}

