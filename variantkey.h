// VariantKey
//
// variantkey.h
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

/**
 * @file variantkey.h
 * @brief VariantKey main functions.
 *
 * The functions provided here allows to generate and process a 64 bit Unsigned Integer Keys for Human Genetic Variants.
 * The VariantKey is sortable for chromosome and position,
 * and it is also fully reversible for variants with up to 11 bases between Reference and Alternate alleles.
 * It can be used to sort, search and match variant-based data easily and very quickly.
 */

#ifndef VARIANTKEY_H
#define VARIANTKEY_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

#define VKMASK_CHROM    0xF800000000000000  //!< VariantKey binary mask for CHROM     [ 11111000 00000000 00000000 00000000 00000000 00000000 00000000 00000000 ]
#define VKMASK_POS      0x07FFFFFF80000000  //!< VariantKey binary mask for POS       [ 00000111 11111111 11111111 11111111 10000000 00000000 00000000 00000000 ]
#define VKMASK_CHROMPOS 0xFFFFFFFF80000000  //!< VariantKey binary mask for CHROM+POS [ 11111111 11111111 11111111 11111111 10000000 00000000 00000000 00000000 ]
#define VKMASK_REFALT   0x000000007FFFFFFF  //!< VariantKey binary mask for REF+ALT   [ 00000000 00000000 00000000 00000000 01111111 11111111 11111111 11111111 ]
#define VKSHIFT_CHROM   59 //!< CHROM LSB position from the VariantKey LSB
#define VKSHIFT_POS     31 //!< POS LSB position from the VariantKey LSB

/**
 * VariantKey struct.
 * Contains the numeically encoded VariantKey components (CHROM, POS, REF+ALT).
 */
typedef struct variantkey_t
{
    uint8_t chrom;   //!< Chromosome encoded number (only the LSB 5 bit are used)
    uint32_t pos;    //!< Reference position, with the first base having position 0 (only the LSB 28 bit are used)
    uint32_t refalt; //!< Code for Reference and Alternate allele (only the LSB 31 bits are used)
} variantkey_t;

/**
 * Struct containing the minimum and maximum VariantKey values for range searches.
 */
typedef struct vkrange_t
{
    uint64_t min; //!< Minimum VariantKey value for any given REF+ALT encoding
    uint64_t max; //!< Maximum VariantKey value for any given REF+ALT encoding
} vkrange_t;

/** @brief Returns chromosome numerical encoding.
 *
 * @param chrom  Chromosome. An identifier from the reference genome, no white-space permitted.
 * @param size   Length of the chrom string, excluding the terminating null byte.
 *
 * @return CHROM code
 */
uint8_t encode_chrom(const char *chrom, size_t size);

/** @brief Decode the chromosome numerical code.
 *
 * @param code   CHROM code.
 * @param chrom  CHROM string buffer to be returned. Its size should be enough to contain the results (max 4 bytes).
 *
 * @return If successful, the total number of characters written is returned,
 *         excluding the null-character appended at the end of the string,
 *         otherwise a negative number is returned in case of failure.
 */
size_t decode_chrom(uint8_t code, char *chrom);

/** @brief Returns reference+alternate numerical encoding.
 *
 * @param ref      Reference allele. String containing a sequence of nucleotide letters.
 *                 The value in the pos field refers to the position of the first nucleotide in the String.
 *                 Characters must be A-Z, a-z or *
 * @param sizeref  Length of the ref string, excluding the terminating null byte.
 * @param alt      Alternate non-reference allele string.
 *                 Characters must be A-Z, a-z or *
 * @param sizealt  Length of the alt string, excluding the terminating null byte.
 *
 * @return REF+ALT code
 */
uint32_t encode_refalt(const char *ref, size_t sizeref, const char *alt, size_t sizealt);

/** @brief Decode the 32 bit REF+ALT code if reversible (if it has 11 or less bases in total and only contains ACGT letters).
 *
 * @param code     REF+ALT code
 * @param ref      REF string buffer to be returned.
 * @param sizeref  Pointer to the size of the ref buffer, excluding the terminating null byte.
 *                 This will contain the final ref size.
 * @param alt      ALT string buffer to be returned.
 * @param sizealt  Pointer to the size of the alt buffer, excluding the terminating null byte.
 *                 This will contain the final alt size.
 *
 * @return      If the code is reversible, then the total number of characters of REF+ALT is returned.
 *              Otherwise 0 is returned.
 */
size_t decode_refalt(uint32_t code, char *ref, size_t *sizeref, char *alt, size_t *sizealt);

/** @brief Returns a 64 bit variant key based on the pre-encoded CHROM, POS (0-based) and REF+ALT.
 *
 * @param chrom      Encoded Chromosome (see encode_chrom).
 * @param pos        Position. The reference position, with the first base having position 0.
 * @param refalt     Encoded Reference + Alternate (see encode_refalt).
 *
 * @return      VariantKey 64 bit code.
 */
uint64_t encode_variantkey(uint8_t chrom, uint32_t pos, uint32_t refalt);

/** @brief Extract the CHROM code from VariantKey.
 *
 * @param vk VariantKey code.
 *
 * @return CHROM code.
 */
uint8_t extract_variantkey_chrom(uint64_t vk);

/** @brief Extract the POS code from VariantKey.
 *
 * @param vk VariantKey code.
 *
 * @return POS code.
 */
uint32_t extract_variantkey_pos(uint64_t vk);

/** @brief Extract the REF+ALT code from VariantKey.
 *
 * @param vk VariantKey code.
 *
 * @return REF+ALT code.
 */
uint32_t extract_variantkey_refalt(uint64_t vk);

/** @brief Decode a VariantKey code and returns the components as variantkey_t structure.
 *
 * @param code VariantKey code.
 * @param vk   Decoded variantkey structure.
 *
 * @return A variantkey_t structure.
 */
void decode_variantkey(uint64_t code, variantkey_t *vk);

/** @brief Returns a 64 bit variant key based on CHROM, POS (0-based), REF, ALT.
 *
 * @param chrom      Chromosome. An identifier from the reference genome, no white-space or leading zeros permitted.
 * @param sizechrom  Length of the chrom string, excluding the terminating null byte.
 * @param pos        Position. The reference position, with the first base having position 0.
 * @param ref        Reference allele. String containing a sequence of nucleotide letters.
 *                   The value in the pos field refers to the position of the first nucleotide in the String.
 *                   Characters must be A-Z, a-z or *
 * @param sizeref    Length of the ref string, excluding the terminating null byte.
 * @param alt        Alternate non-reference allele string.
 *                   Characters must be A-Z, a-z or *
 * @param sizealt    Length of the alt string, excluding the terminating null byte.
 *
 * @return      VariantKey 64 bit code.
 */
uint64_t variantkey(const char *chrom, size_t sizechrom, uint32_t pos, const char *ref, size_t sizeref, const char *alt, size_t sizealt);

/** @brief Returns minimum and maximum VariantKeys for range searches.
 *
 * @param chrom     Chromosome encoded number.
 * @param pos_min   Start reference position, with the first base having position 0.
 * @param pos_max   End reference position, with the first base having position 0.
 * @param range     VariantKey range values.
 *
 * @return      Min and Max variant keys for any given REF+ALT encoding.
 */
void variantkey_range(uint8_t chrom, uint32_t pos_min, uint32_t pos_max, vkrange_t *range);

/** @brief Compares two VariantKeys by chromosome only.
 *
 * @param vka    The first VariantKey to be compared.
 * @param vkb    The second VariantKey to be compared.
 *
 * @return -1 if the first chromosome is smaller than the second, 0 if they are equal and 1 if the first is greater than the second.
 */
int compare_variantkey_chrom(uint64_t vka, uint64_t vkb);

/** @brief Compares two VariantKeys by chromosome and position.
 *
 * @param vka    The first VariantKey to be compared.
 * @param vkb    The second VariantKey to be compared.
 *
 * @return -1 if the first CHROM+POS is smaller than the second, 0 if they are equal and 1 if the first is greater than the second.
 */
int compare_variantkey_chrom_pos(uint64_t vka, uint64_t vkb);

/** @brief Returns VariantKey hexadecimal string (16 characters).
 *
 * The string represent a 64 bit number or:
 *   -  5 bit for CHROM
 *   - 28 bit for POS
 *   - 31 bit for REF+ALT
 *
 * @param vk    VariantKey code.
 * @param str   String buffer to be returned (it must be sized 17 bytes at least).
 *
 * @return      Upon successful return, these function returns the number of characters processed
 *              (excluding the null byte used to end output to strings).
 *              If the buffer size is not sufficient, then the return value is the number of characters required for
 *              buffer string, including the terminating null byte.
 */
size_t variantkey_hex(uint64_t vk, char *str);

/** @brief Parses a VariantKey hexadecimal string and returns the code.
 *
 * @param vs    VariantKey hexadecimal string (it must contain 16 hexadecimal characters).
 *
 * @return A VariantKey code.
 */
uint64_t parse_variantkey_hex(const char *vs);

#ifdef __cplusplus
}
#endif

#endif  // VARIANTKEY_H
