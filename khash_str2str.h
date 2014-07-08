#ifndef __KHASH_STR2STR_H__
#define __KHASH_STR2STR_H__

#include <htslib/khash.h>

KHASH_MAP_INIT_STR(str2str, const char*)

/*
 *  Wrappers for khash dictionaries used by mpileup. 
 */

static inline void *khash_str2str_init(void)
{
    return kh_init(str2str);
}

/*
 *  Destroy the hash structure, but not the keys
 */ 
static inline void khash_str2str_destroy(void *_hash)
{
    khash_t(str2str) *hash = (khash_t(str2str)*)_hash;
    if (hash) kh_destroy(str2str, hash); // Note that strings are not freed.
}

/*
 *  Destroys both the hash structure and the keys
 */ 
static inline void khash_str2str_destroy_free(void *_hash)
{
    khash_t(str2str) *hash = (khash_t(str2str)*)_hash;
    khint_t k;
    if (hash == 0) return;
    for (k = 0; k < kh_end(hash); ++k)
        if (kh_exist(hash, k)) free((char*)kh_key(hash, k));
    kh_destroy(str2str, hash);
}

/*
 *  Returns value if key exists or NULL if not
 */
static inline char *khash_str2str_get(void *_hash, const char *str)
{
    khash_t(str2str) *hash = (khash_t(str2str)*)_hash;
    khint_t k = kh_get(str2str, hash, str);
    if ( k == kh_end(hash) ) return NULL;
    return (char*)kh_val(hash, k);
}

/*
 *  Set a new key,value pair. On success returns the bin index, on
 *  error -1 is returned.
 */
static inline int khash_str2str_set(void *_hash, const char *str, const char *value)
{
    khint_t k;
    int ret;
    khash_t(str2str) *hash = (khash_t(str2str)*)_hash;
    if ( !hash ) return -1;
    k = kh_put(str2str, hash, str, &ret);
    kh_val(hash,k) = value;
    return k;
}

#endif
