// VariantKey
//
// variantkey.c
//
// @category   Tools
// @author     Nicola Asuni <nicola.asuni@genomicsplc.com>
// @copyright  2017-2018 GENOMICS plc
// @license    MIT (see LICENSE)
// @link       https://github.com/genomicsplc/variantkey
//
// LICENSE
//
// Copyright (c) 2017-2018 GENOMICS plc
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
// VariantKey by Nicola Asuni

#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include "variantkey.h"

static inline int aztoupper(int c)
{
    if (c >= 'a')
    {
        return (c ^ ('a' - 'A'));
    }
    return c;
}

inline uint8_t encode_chrom(const char *chrom, size_t size)
{
    // X > 23 ; Y > 24 ; M > 25
    static const uint8_t onecharmap[] =
    {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        /*                                    M                                X  Y                  */
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,23,24, 0, 0, 0, 0, 0, 0,
        /*                                    m                                x  y                  */
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,23,24, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };
    // remove "chr" prefix
    if ((size > 3)
            && ((chrom[0] == 'c') || (chrom[0] == 'C'))
            && ((chrom[1] == 'h') || (chrom[1] == 'H'))
            && ((chrom[2] == 'r') || (chrom[2] == 'R')))
    {
        chrom += 3;
        size -= 3;
    }
    if (size == 0)
    {
        return 0;
    }
    if ((chrom[0] <= '9') && (chrom[0] >= '0')) // Number
    {
        size_t i;
        uint8_t v = (chrom[0] - '0');
        for (i = 1; i < size; i++)
        {
            if ((chrom[i] > '9') || (chrom[i] < '0'))
            {
                return 0; // NA
            }
            v = ((v * 10) + (chrom[i] - '0'));
        }
        return v;
    }
    if ((size == 1) || ((size == 2) && ((chrom[1] == 'T') || (chrom[1] == 't'))))
    {
        return onecharmap[((unsigned char)chrom[0])];
    }
    return 0; // NA
}

inline size_t decode_chrom(uint8_t code, char *chrom)
{
    if ((code < 1) || (code > 25))
    {
        return sprintf(chrom, "NA");
    }
    if (code < 23)
    {
        return sprintf(chrom, "%" PRIu8, code);
    }
    static const char *map[] = {"X", "Y", "MT"};
    return sprintf(chrom, "%s", map[(code - 23)]);
}

static inline uint32_t encode_base(const unsigned char c)
{
    /*
      Encode base:

      A > 0
      C > 1
      G > 2
      T > 3
    */
    static const uint32_t map[] = {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                                   4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                                   /*A   C       G                         T*/
                                   4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
                                   /*a   c       g                         t*/
                                   4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
                                   4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                                   4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                                   4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                                   4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                                  };
    return map[c];
}

static inline int encode_allele(uint32_t *h, uint8_t *bitpos, const char *str, size_t size)
{
    int c;
    uint32_t v;
    while ((c = *str++) && (size--))
    {
        v = encode_base(c);
        if (v > 3)
        {
            return -1;
        }
        *bitpos -= 2;
        *h |= (v << *bitpos); // A will be coded as 1
    }
    return 0;
}

static inline uint32_t encode_refalt_rev(const char *ref, size_t sizeref, const char *alt, size_t sizealt)
{
    //[******** ******** ******** ******** *RRRRAAA A1122334 45566778 8990011*]
    uint32_t h = 0;
    h |= ((uint32_t)(sizeref) << 27); // RRRR: length of (REF - 1)
    h |= ((uint32_t)(sizealt) << 23); // AAAA: length of (ALT - 1)
    uint8_t bitpos = 23;
    if ((encode_allele(&h, &bitpos, ref, sizeref) < 0) || (encode_allele(&h, &bitpos, alt, sizealt) < 0))
    {
        return 0; // error code
    }
    return h;
}

static inline uint32_t pack_chars(const char *str, size_t size)
{
    int c;
    uint32_t h = 0;
    uint8_t bitpos = VKSHIFT_POS;
    while ((c = aztoupper(*str++)) && (size--))
    {
        if (c == '*')
        {
            c = ('Z' + 1);
        }
        bitpos -= 5;
        h |= ((c - 'A' + 1) << bitpos); // 'A' will be coded as 1
    }
    return h;
}

// Mix two 32 bit hash numbers using the MurmurHash3 algorithm
static inline uint32_t muxhash(uint32_t k, uint32_t h)
{
    k *= 0xcc9e2d51;
    k = (k >> 17) | (k << (32 - 17));
    k *= 0x1b873593;
    h ^= k;
    h = (h >> 19) | (h << (32 - 19));
    return ((h * 5) + 0xe6546b64);
}

// Return a 32 bit hash of a nucleotide string
static inline uint32_t hash32(const char *str, size_t size)
{
    uint32_t h = 0;
    size_t len = 6;
    while (size > 0)
    {
        if (size < len)
        {
            len = size;
        }
        // [ 01111122 22233333 44444555 55666660 ]
        // pack blocks of 6 characters in 32 bit (6 x 5 bit + 2 spare bit)
        h = muxhash(pack_chars(str, len), h);
        size -= len;
        str += len;
    }
    return h;
}

static inline uint32_t encode_refalt_hash(const char *ref, size_t sizeref, const char *alt, size_t sizealt)
{
    // 0x3 is the separator character between REF and ALT [00000000 00000000 00000000 00000011]
    uint32_t h = muxhash(hash32(alt, sizealt), muxhash(0x3, hash32(ref, sizeref)));
    // finalization mix - MurmurHash3 algorithm
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return ((h >> 1) | 0x1); // 0x1 is the set bit to indicate HASH mode [00000000 00000000 00000000 00000001]
}

inline uint32_t encode_refalt(const char *ref, size_t sizeref, const char *alt, size_t sizealt)
{
    if ((sizeref + sizealt) <= 11)
    {
        uint32_t h = encode_refalt_rev(ref, sizeref, alt, sizealt);
        if (h != 0)
        {
            return h;
        }
    }
    return encode_refalt_hash(ref, sizeref, alt, sizealt);
}

static inline char decode_base(uint32_t code, int bitpos)
{
    static const char base[4] = {'A', 'C', 'G', 'T'};
    return base[((code >> bitpos) & 0x3)]; // 0x3 is the 2 bit mask [00000011]
}

static inline size_t decode_refalt_rev(uint32_t code, char *ref, size_t *sizeref, char *alt, size_t *sizealt)
{
    *sizeref = (size_t)((code & 0x78000000) >> 27); // [01111000 00000000 00000000 00000000]
    *sizealt = (size_t)((code & 0x07800000) >> 23); // [00000111 10000000 00000000 00000000]
    uint8_t bitpos = 23;
    size_t i = 0;
    for(i = 0; i < *sizeref; i++)
    {
        bitpos -= 2;
        ref[i] = decode_base(code, bitpos);
    }
    ref[i] = 0;
    for(i = 0; i < *sizealt; i++)
    {
        bitpos -= 2;
        alt[i] = decode_base(code, bitpos);
    }
    alt[i] = 0;
    return (*sizeref + *sizealt);
}

inline size_t decode_refalt(uint32_t code, char *ref, size_t *sizeref, char *alt, size_t *sizealt)
{
    if (code & 0x1) // check last bit
    {
        return 0; // non-reversible encoding
    }
    return decode_refalt_rev(code, ref, sizeref, alt, sizealt);
}

inline uint64_t encode_variantkey(uint8_t chrom, uint32_t pos, uint32_t refalt)
{
    return (((uint64_t)chrom << VKSHIFT_CHROM) | ((uint64_t)pos << VKSHIFT_POS) | (uint64_t)refalt);
}

inline uint8_t extract_variantkey_chrom(uint64_t vk)
{
    return (uint8_t)((vk & VKMASK_CHROM) >> VKSHIFT_CHROM);
}

inline uint32_t extract_variantkey_pos(uint64_t vk)
{
    return (uint32_t)((vk & VKMASK_POS) >> VKSHIFT_POS);
}

inline uint32_t extract_variantkey_refalt(uint64_t vk)
{
    return (uint32_t)(vk & VKMASK_REFALT);
}

inline void decode_variantkey(uint64_t code, variantkey_t *vk)
{
    vk->chrom = extract_variantkey_chrom(code);
    vk->pos = extract_variantkey_pos(code);
    vk->refalt = extract_variantkey_refalt(code);
}

inline uint64_t variantkey(const char *chrom, size_t sizechrom, uint32_t pos, const char *ref, size_t sizeref, const char *alt, size_t sizealt)
{
    return encode_variantkey(encode_chrom(chrom, sizechrom), pos, encode_refalt(ref, sizeref, alt, sizealt));
}

inline void variantkey_range(uint8_t chrom, uint32_t pos_min, uint32_t pos_max, vkrange_t *range)
{
    uint64_t c = ((uint64_t)chrom << VKSHIFT_CHROM);
    range->min = (c | ((uint64_t)pos_min << VKSHIFT_POS));
    range->max = (c | ((uint64_t)pos_max << VKSHIFT_POS) | VKMASK_REFALT);
}

static inline int compare_uint64_t(uint64_t a, uint64_t b)
{
    return (a < b) ? -1 : (a > b);
}

inline int compare_variantkey_chrom(uint64_t vka, uint64_t vkb)
{
    return compare_uint64_t((vka >> VKSHIFT_CHROM), (vkb >> VKSHIFT_CHROM));
}

inline int compare_variantkey_chrom_pos(uint64_t vka, uint64_t vkb)
{
    return compare_uint64_t((vka >> VKSHIFT_POS), (vkb >> VKSHIFT_POS));
}

inline size_t variantkey_hex(uint64_t vk, char *str)
{
    return sprintf(str, "%016" PRIx64, vk);
}

inline uint64_t parse_variantkey_hex(const char *vs)
{
    uint64_t v = 0;
    uint8_t b;
    size_t i;
    for (i = 0; i < 16; i++)
    {
        b = vs[i];
        if (b >= 'a')
        {
            b -= ('a' - 10); // a-f
        }
        else
        {
            if (b >= 'A')
            {
                b -= ('A' - 10); // A-F
            }
            else
            {
                b -= '0'; // 0-9
            }
        }
        v = ((v << 4) | b);
    }
    return v;
}
