/* The MIT License

   Copyright (c) 2016 Genome Research Ltd.

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

/*
    Read about transcript types here
        http://vega.sanger.ac.uk/info/about/gene_and_transcript_types.html
        http://www.ensembl.org/info/genome/variation/predicted_data.html
        http://www.gencodegenes.org/gencode_biotypes.html

    List of supported biotypes
        antisense
        IG_C_gene
        IG_D_gene
        IG_J_gene
        IG_LV_gene
        IG_V_gene
        lincRNA
        macro_lncRNA
        miRNA
        misc_RNA
        Mt_rRNA
        Mt_tRNA
        polymorphic_pseudogene
        processed_transcript
        protein_coding
        ribozyme
        rRNA
        sRNA
        scRNA
        scaRNA
        sense_intronic
        sense_overlapping
        snRNA
        snoRNA
        TR_C_gene
        TR_D_gene
        TR_J_gene
        TR_V_gene

    The gff parsing logic
        We collect features such by combining gff lines A,B,C as follows:
            A .. gene line with a supported biotype
                    A.ID=~/^gene:/

            B .. transcript line referencing A
                    B.ID=~/^transcript:/ && B.Parent=~/^gene:A.ID/

            C .. corresponding CDS, exon, and UTR lines:
                    C[3] in {"CDS","exon","three_prime_UTR","five_prime_UTR"} && C.Parent=~/^transcript:B.ID/ 

        For coding biotypes ("protein_coding" or "polymorphic_pseudogene") the
        complete chain link C -> B -> A is required. For the rest, link B -> A suffices.
        
                
    The supported consequence types, sorted by impact:
        splice_acceptor_variant .. end region of an intron changed (2bp at the 3' end of an intron)
        splice_donor_variant    .. start region of an intron changed (2bp at the 5' end of an intron)
        stop_gained             .. DNA sequence variant resulting in a stop codon
        frameshift_variant      .. number of inserted/deleted bases not a multiple of three, disrupted translational frame
        stop_lost               .. elongated transcript, stop codon changed
        start_lost              .. the first codon changed
        inframe_insertion       .. inserted coding sequence, unchanged reading frame
        inframe_deletion        .. deleted coding sequence, unchanged reading frame
        missense_variant        .. amino acid (aa) change, unchanged length
        splice_region_variant   .. change within 1-3 bases of the exon or 3-8 bases of the intron
        synonymous_variant      .. DNA sequence variant resulting in no amino acid change
        stop_retained_variant   .. different stop codon
        non_coding_variant      .. variant in non-coding sequence, such as RNA gene
        5_prime_UTR_variant
        3_prime_UTR_variant
        intron_variant          .. reported only if none of the above
        intergenic_variant      .. reported only if none of the above


    The annotation algorithm.
        The algorithm checks if the variant falls in a region of a supported type. The
        search is performed in the following order, until a match is found:
            1. idx_cds(gf_cds_t) - lookup CDS by position, create haplotypes, call consequences
            2. idx_utr(gf_utr_t) - check UTR hits
            3. idx_exon(gf_exon_t) - check for splice variants
            4. idx_tscript(tscript_t) - check for intronic variants, RNAs, etc.

        These regidx indexes are created by parsing a gff3 file as follows:
            1.  create the array "ftr" of all UTR, CDS, exons. This will be
            processed later and pruned based on transcript types we want to keep.
            In the same go, create the hash "id2tr" of transcripts to keep
            (based on biotype) which maps from transcript_id to a transcript. At
            the same time also build the hash "gid2gene" which maps from gene_id to
            gf_gene_t pointer.
            
            2.  build "idx_cds", "idx_tscript", "idx_utr" and "idx_exon" indexes.
            Use only features from "ftr" which are present in "id2tr".

            3.  clean data that won't be needed anymore: ftr, id2tr, gid2gene.
        
    Data structures.
        idx_cds, idx_utr, idx_exon, idx_tscript:
            as described above, regidx structures for fast lookup of exons/transcripts
            overlapping a region, the payload is a pointer to tscript.cds
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/khash.h>
#include <htslib/khash_str2int.h>
#include <htslib/kseq.h>
#include <htslib/faidx.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include "bcftools.h"
#include "filter.h"
#include "regidx.h"
#include "kheap.h"
#include "smpl_ilist.h"
#include "rbuf.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

// Definition of splice_region, splice_acceptor and splice_donor
#define N_SPLICE_DONOR         2      
#define N_SPLICE_REGION_EXON   3 
#define N_SPLICE_REGION_INTRON 8 

// Ensembl ID format, e.g. 
//     ENST00000423372 for human .. ENST%011d
//  ENSMUST00000120394 for mouse .. ENSMUST%011d
char  ENSID_BUF[32], *ENSID_FMT = NULL;
static inline char *ENSID(uint32_t id)
{
    sprintf(ENSID_BUF,ENSID_FMT,id);
    return ENSID_BUF;
}


#define STRAND_REV 0
#define STRAND_FWD 1

// How to treat phased/unphased genotypes
#define PHASE_REQUIRE 0     // --phase r
#define PHASE_MERGE   1     // --phase m
#define PHASE_AS_IS   2     // --phase a
#define PHASE_SKIP    3     // --phase s
#define PHASE_DROP_GT 4     // --samples -

// Node types in the haplotype tree
#define HAP_CDS   0
#define HAP_ROOT  1 
#define HAP_SSS   2     // start/stop/splice

#define CSQ_PRINTED_UPSTREAM    (1<<0)
#define CSQ_SYNONYMOUS_VARIANT  (1<<1)
#define CSQ_MISSENSE_VARIANT    (1<<2)
#define CSQ_STOP_LOST           (1<<3)      // todo: example
#define CSQ_STOP_GAINED         (1<<4)      // todo: example
#define CSQ_INFRAME_DELETION    (1<<5)
#define CSQ_INFRAME_INSERTION   (1<<6)
#define CSQ_FRAMESHIFT_VARIANT  (1<<7)
#define CSQ_SPLICE_ACCEPTOR     (1<<8)
#define CSQ_SPLICE_DONOR        (1<<9)
#define CSQ_START_LOST          (1<<10)
#define CSQ_SPLICE_REGION       (1<<11)     // todo: introns, not the first/last exon
//#define CSQ_STOP_RETAINED       (1<<12)
#define CSQ_UTR5                (1<<13)
#define CSQ_UTR3                (1<<14)
#define CSQ_NON_CODING          (1<<15)
#define CSQ_INTRON              (1<<16)
//#define CSQ_INTERGENIC          (1<<17)

// see kput_csq()
const char *csq_strings[] = 
{
    NULL, 
    "synonymous", 
    "missense", 
    "stop_lost", 
    "stop_gained", 
    "inframe_deletion", 
    "inframe_insertion", 
    "frameshift", 
    "splice_acceptor", 
    "splice_donor", 
    "start_lost", 
    "splice_region", 
    "stop_retained", 
    "5_prime_utr", 
    "3_prime_utr", 
    "non_coding", 
    "intron", 
    "intergenic" 
};


// GFF line types
#define GFF_TSCRIPT_LINE 1
#define GFF_GENE_LINE    2


/* 
    Genomic features, for fast lookup by position to overlapping features
*/
#define GF_coding_bit 6
#define GF_is_coding(x) ((x) & (1<<GF_coding_bit))
#define GF_MT_rRNA                       1
#define GF_MT_tRNA                       2
#define GF_lincRNA                       3
#define GF_miRNA                         4
#define GF_MISC_RNA                      5
#define GF_rRNA                          6
#define GF_snRNA                         7
#define GF_snoRNA                        8
#define GF_PROCESSED_TRANSCRIPT          9
#define GF_ANTISENSE                    10
#define GF_macro_lncRNA                 11
#define GF_ribozyme                     12
#define GF_sRNA                         13
#define GF_scRNA                        14
#define GF_scaRNA                       15
#define GF_SENSE_INTRONIC               16
#define GF_SENSE_OVERLAPPING            17
#define GF_PROTEIN_CODING               (1|(1<<GF_coding_bit))  // 65, 66, ...
#define GF_POLYMORPHIC_PSEUDOGENE       (2|(1<<GF_coding_bit))
#define GF_IG_C                         (3|(1<<GF_coding_bit))
#define GF_IG_D                         (4|(1<<GF_coding_bit))
#define GF_IG_J                         (5|(1<<GF_coding_bit))
#define GF_IG_LV                        (6|(1<<GF_coding_bit))
#define GF_IG_V                         (7|(1<<GF_coding_bit))
#define GF_TR_C                         (8|(1<<GF_coding_bit))
#define GF_TR_D                         (9|(1<<GF_coding_bit))
#define GF_TR_J                        (10|(1<<GF_coding_bit))
#define GF_TR_V                        (11|(1<<GF_coding_bit))
#define GF_CDS      ((1<<(GF_coding_bit+1))+1)      // 129, 130, ...
#define GF_EXON     ((1<<(GF_coding_bit+1))+2)
#define GF_UTR3     ((1<<(GF_coding_bit+1))+3)
#define GF_UTR5     ((1<<(GF_coding_bit+1))+4)
// GF_MAX = (1<<30)-1, see hap_node_t

typedef struct _tscript_t tscript_t;
typedef struct
{
    tscript_t *tr;      // transcript
    uint32_t beg;       // the start coordinate of the CDS (on the reference strand, 0-based)
    uint32_t pos;       // 0-based index of the first exon base within the transcript (only to
                        //  update hap_node_t.sbeg in hap_init, could be calculated on the fly)
    uint32_t len;       // exon length
    uint32_t icds:30,   // exon index within the transcript
             phase:2;   // offset of the CDS
}
gf_cds_t;
typedef struct
{
    char *name;           // human readable name, e.g. ORF45
    uint8_t iseq;
}
gf_gene_t;
typedef struct
{
    uint32_t beg,end;
    tscript_t *tr;
}
gf_exon_t;
typedef enum { prime3, prime5 } utr_t;
typedef struct
{
    utr_t which;
    uint32_t beg,end;
    tscript_t *tr;
}
gf_utr_t;


/*
    Structures related to VCF output:

    vcrec_t 
        single VCF record

    csq_t
        a consequence tied to a VCF record

    vbuf_t
    pos2vbuf
        VCF records with the same position clustered together for a fast lookup via pos2vbuf
*/
typedef struct _vbuf_t vbuf_t;
typedef struct
{
    bcf1_t *line;
    uint32_t *smpl;     // bitmask of sample consequences with first/second haplotype interleaved
    int ncsq, mcsq;
    char **csq;         // pointers to hap_node.csq.csq or args.csq_buf
}
vrec_t;
typedef struct
{
    uint32_t pos;
    int type;       // one of CSQ_* types
    char *str;      // e.g. 9T>9KT|126C>CTT,127G>GG
    vrec_t *vrec;   // the csq is tied to this vcf line
    int idx;        // 0-based index of the csq at the VCF line
}
csq_t;
struct _vbuf_t
{
    vrec_t **vrec;   // buffer of VCF lines with the same position
    int n, m;
};
KHASH_MAP_INIT_INT(pos2vbuf, vbuf_t*)


/*
    Structures related to haplotype-aware consequences in coding regions

    hap_node_t
        node of a haplotype tree. Each transcript has one tree

    tscript_t
        despite its general name, it is intended for coding transcripts only

    hap_t
    hstack_t
        for traversal of the haplotype tree and braking combined
        consequences into independent parts
*/
typedef struct _hap_node_t hap_node_t;
struct _hap_node_t
{
    char *seq;          // cds segment [parent_node,this_node)
    char *var;          // variant "ref>alt"
    uint32_t type:2,    // HAP_ROOT or HAP_CDS
             csq:30;    // this node's consequence
    int dlen;           // alt minus ref length: <0 del, >0 ins, 0 substitution
    uint32_t rbeg;      // variant's VCF position (0-based, inclusive)
    int32_t rlen;       // variant's rlen; alen=rlen+dlen; fake for non CDS types
    uint32_t sbeg;      // variant's position on the spliced reference transcript (0-based, inclusive)
    uint32_t icds;      // which exon does this node's variant overlaps
    hap_node_t **child, *prev;  // children haplotypes and previous coding node
    int nchild, mchild;
    bcf1_t *cur_rec, *rec;      // current VCF record and node's VCF record
    uint32_t nend;              // number of haplotypes ending in this node
    int *cur_child, mcur_child; // mapping from the allele to the currently active child
    csq_t *csq_list;            // list of haplotype's consequences, broken by position
    int ncsq_list, mcsq_list;
};
struct _tscript_t
{
    uint32_t id;        // transcript id
    uint32_t beg,end;   // transcript's beg and end coordinate (ref strand, 0-based, inclusive)
    uint32_t strand:1,  // STRAND_REV or STRAND_FWD
             ncds:31,   // number of exons
             mcds;
    gf_cds_t **cds;     // ordered list of exons
    char *ref;          // reference sequence
    char *sref;         // spliced reference sequence
    hap_node_t *root;   // root of the haplotype tree
    hap_node_t **hap;   // pointer to haplotype leaves, two for each sample
    int nhap, nsref;    // number of haplotypes and length of sref
    uint32_t type;      // one of GF_* types
    gf_gene_t *gene;
};
static inline int cmp_tscript(tscript_t **a, tscript_t **b)
{
    return ( (*a)->end  < (*b)->end ) ? 1 : 0;
}
KHEAP_INIT(trhp, tscript_t*, cmp_tscript)
typedef khp_trhp_t tr_heap_t;
typedef struct
{
    hap_node_t *node;   // current node
    int ichild;         // current child in the active node
    int dlen;           // total dlen, from the root to the active node
    size_t slen;        // total sequence length, from the root to the active node
}
hstack_t;
typedef struct
{
    int mstack;
    hstack_t *stack;
    tscript_t *tr;      // tr->ref: spliced transcript on ref strand
    kstring_t sseq;     // spliced haplotype sequence on ref strand
    kstring_t tseq;     // the variable part of translated haplotype transcript, coding strand
    kstring_t tref;     // the variable part of translated reference transcript, coding strand
    int upstream_stop;
    uint32_t sbeg;      // stack's sbeg, for cases first node's type is HAP_SSS
}
hap_t;


/*
    Helper structures, only for initialization
    
    ftr_t
        temporary list of all exons, CDS, UTRs 
*/
KHASH_MAP_INIT_INT(int2tscript, tscript_t*)
KHASH_MAP_INIT_INT(int2int, int)
KHASH_MAP_INIT_INT(int2gene, gf_gene_t*)
typedef struct
{
    int type;       // GF_CDS, GF_EXON, GF_5UTR, GF_3UTR
    uint32_t beg;
    uint32_t end;
    uint32_t trid;
    uint32_t strand:1;   // STRAND_REV,STRAND_FWD
    uint32_t phase:2;    // 0, 1 or 2
    uint32_t iseq:29;
}
ftr_t;
typedef struct
{
    // all exons, CDS, UTRs
    ftr_t *ftr;
    int nftr, mftr;

    // mapping from transcript ensembl id to gene id
    kh_int2gene_t *gid2gene;

    // mapping from transcript id to tscript, for quick CDS anchoring
    kh_int2tscript_t *id2tr;

    // sequences
    void *seq2int;
    char **seq;
    int nseq, mseq;
}
aux_t;

typedef struct _args_t
{
    // the main regidx lookups, from chr:beg-end to overlapping features and
    // index iterator
    regidx_t *idx_cds, *idx_utr, *idx_exon, *idx_tscript;
    regitr_t *itr;

    // temporary structures, deleted after initializtion
    aux_t init;

    // text tab-delimited output (out) or vcf/bcf output (out_fh)
    FILE *out;
    htsFile *out_fh;

    // vcf
    bcf_srs_t *sr;
    bcf_hdr_t *hdr;
    int hdr_nsmpl;          // actual number of samples in the vcf, for bcf_update_format_values()

    // include or exclude sites which match the filters
    filter_t *filter;
    char *filter_str;
    int filter_logic;       // FLT_INCLUDE or FLT_EXCLUDE

    // samples to process
    int sample_is_file;
    char *sample_list;
    smpl_ilist_t *smpl;

    char *outdir, **argv, *fa_fname, *gff_fname, *output_fname;
    int argc, output_type;
    int phase, quiet;
    
    int rid;                    // current chromosome
    tr_heap_t *active_tr;       // heap of active transcripts for quick flushing
    hap_t *hap;                 // transcript haplotype recursion
    vbuf_t **vcf_buf;           // buffered VCF lines to annotate with CSQ and flush
    rbuf_t vcf_rbuf;            // round buffer indexes to vcf_buf
    kh_pos2vbuf_t *pos2vbuf;    // fast lookup of buffered lines by position
    tscript_t **rm_tr;          // buffer of transcripts to clean
    int nrm_tr, mrm_tr;
    csq_t *csq_buf;             // pool of csq not managed by hap_node_t, i.e. non-CDS csqs
    int ncsq_buf, mcsq_buf;

    faidx_t *fai;
    kstring_t str, str2;
    int32_t *gt_arr, mgt_arr;
}
args_t;

// AAA, AAC, ...
const char *gencode = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
const uint8_t nt4[] =
{
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 3,4,4,4, 4,4,4,4, 4,4,4,4,
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 3
};
const uint8_t cnt4[] =
{
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,3,4,2, 4,4,4,1, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 0,4,4,4, 4,4,4,4, 4,4,4,4,
    4,3,4,2, 4,4,4,1, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 0
};
#define dna2aa(x)  gencode[  nt4[(uint8_t)(x)[0]]<<4 |  nt4[(uint8_t)(x)[1]]<<2 |  nt4[(uint8_t)(x)[2]] ]
#define cdna2aa(x) gencode[ cnt4[(uint8_t)(x)[2]]<<4 | cnt4[(uint8_t)(x)[1]]<<2 | cnt4[(uint8_t)(x)[0]] ]


const char *gf_type2gff_string(int type)
{
    if ( !GF_is_coding(type) )
    {
        if ( type < (1<<(GF_coding_bit+1)) )
        {
            char *strings[] = 
            { 
                "", "MT_rRNA", "MT_tRNA", "lincRNA", "miRNA", "misc_RNA", "rRNA", "snRNA", "snoRNA", "processed_transcript",
                "antisense", "macro_lncRNA", "ribozyme", "sRNA", "scRNA", "scaRNA", "sense_intronic", "sense_overlapping" 
            };
            return strings[type];
        }
        char *strings[] = { "", "CDS", "exon", "3_prime_UTR", "5_prime_UTR" };
        type &= (1<<(GF_coding_bit+1)) - 1;
        return strings[type];
    }
    char *strings[] = { "", "protein_coding", "polymorphic_pseudogene", "IG_C", "IG_D", "IG_J", "IG_V", "TR_C", "TR_D", "TR_J", "TR_V", };
    type &= (1<<GF_coding_bit) - 1;
    return strings[type];
}

/*
    gff parsing functions
*/
static inline int feature_set_seq(args_t *args, char *chr_beg, char *chr_end)
{
    aux_t *aux = &args->init;
    char c = chr_end[1];
    chr_end[1] = 0;
    int iseq;
    if ( khash_str2int_get(aux->seq2int, chr_beg, &iseq)!=0 )
    {
        hts_expand(char*, aux->nseq+1, aux->mseq, aux->seq);
        aux->seq[aux->nseq] = strdup(chr_beg);
        iseq = khash_str2int_inc(aux->seq2int, aux->seq[aux->nseq]);
        aux->nseq++;
        assert( aux->nseq < 256 );  // see gf_gene_t.iseq
    }
    chr_end[1] = c;
    return iseq;
}
static inline char *gff_skip(const char *line, char *ss)
{
    while ( *ss && *ss!='\t' ) ss++;
    if ( !*ss ) error("[%s:%d %s] Could not parse the line: %s\n",__FILE__,__LINE__,__FUNCTION__,line);
    return ss+1;
}
static inline void gff_parse_chr(const char *line, char **chr_beg, char **chr_end)
{
    char *se = (char*) line;
    while ( *se && *se!='\t' ) se++;
    if ( !*se ) error("[%s:%d %s] Could not parse the line: %s\n",__FILE__,__LINE__,__FUNCTION__,line);
    *chr_beg = (char*) line;
    *chr_end = se-1;
}
static inline char *gff_parse_beg_end(const char *line, char *ss, uint32_t *beg, uint32_t *end)
{
    char *se = ss;
    *beg = strtol(ss, &se, 10) - 1;
    if ( ss==se ) error("[%s:%d %s] Could not parse the line:\n\t%s\n\t%s\n",__FILE__,__LINE__,__FUNCTION__,line,ss);
    ss = se+1;
    *end = strtol(ss, &se, 10) - 1;
    if ( ss==se ) error("[%s:%d %s] Could not parse the line: %s\n",__FILE__,__LINE__,__FUNCTION__,line);
    return se+1;
}
static inline uint32_t gff_parse_id(const char *line, const char *needle, char *ss)
{
    ss = strstr(ss,needle);
    if ( !ss ) error("[%s:%d %s] Could not parse the line, \"%s\" not present: %s\n",__FILE__,__LINE__,__FUNCTION__,needle,line);
    ss += strlen(needle);
    while ( *ss && !isdigit(*ss) ) ss++;
    if ( !ss ) error("[%s:%d %s] Could not parse the line: %s\n",__FILE__,__LINE__,__FUNCTION__, line);
    char *se;
    uint32_t id = strtol(ss, &se, 10);
    if ( ss==se ) error("[%s:%d %s] Could not parse the line: %s\n",__FILE__,__LINE__,__FUNCTION__, line);
    if ( *se && *se!=';' && *se!='\t' ) error("[%s:%d %s] Could not parse the line: %s\n",__FILE__,__LINE__,__FUNCTION__,line);
    assert( id <= 0xffffff );   // see gf_gene_t.id. Ensembl IDs are never that big in practice
    return id;
}
static void gff_parse_ensid_fmt(const char *line, const char *needle, char *ss)
{
    ss = strstr(ss,needle);
    if ( !ss ) error("[%s:%d %s] Could not parse the line, \"%s\" not present: %s\n",__FILE__,__LINE__,__FUNCTION__,needle,line);
    ss += strlen(needle);
    char *se = ss;
    while ( *se && !isdigit(*se) ) se++;
    kstring_t str = {0,0,0};
    kputsn(ss,se-ss,&str);
    ss = se;
    while ( *se && isdigit(*se) ) se++;
    ksprintf(&str,"%%0%dd",(int)(se-ss));
    ENSID_FMT = str.s;
}
static inline int gff_parse_type(char *line)
{
    line = strstr(line,"ID=");
    if ( !line ) return -1;
    line += 3;
    if ( !strncmp(line,"transcript:",11) ) return GFF_TSCRIPT_LINE;
    else if ( !strncmp(line,"gene:",5) ) return GFF_GENE_LINE;
    return -1;
}
static inline int gff_parse_biotype(char *_line)
{
    char *line = strstr(_line,"biotype=");
    if ( !line ) return -1;

    line += 8;
    switch (*line)
    {
        case 'p': 
            if ( !strncmp(line,"protein_coding",14) ) return GF_PROTEIN_CODING;
            else if ( !strncmp(line,"polymorphic_pseudogene",22) ) return GF_POLYMORPHIC_PSEUDOGENE;
            else if ( !strncmp(line,"processed_transcript",20) ) return GF_PROCESSED_TRANSCRIPT;
            break;
        case 'a':
            if ( !strncmp(line,"antisense",9) ) return GF_ANTISENSE;
        case 'I':
            if ( !strncmp(line,"IG_C_gene",9) ) return GF_IG_C;
            else if ( !strncmp(line,"IG_D_gene",9) ) return GF_IG_D;
            else if ( !strncmp(line,"IG_J_gene",9) ) return GF_IG_J;
            else if ( !strncmp(line,"IG_LV_gene",10) ) return GF_IG_LV;
            else if ( !strncmp(line,"IG_V_gene",9) ) return GF_IG_V;
            break;
        case 'T':
            if ( !strncmp(line,"TR_C_gene",5) ) return GF_TR_C;
            else if ( !strncmp(line,"TR_D_gene",5) ) return GF_TR_D;
            else if ( !strncmp(line,"TR_J_gene",5) ) return GF_TR_J;
            else if ( !strncmp(line,"TR_V_gene",5) ) return GF_TR_V;
            break;
        case 'M':
            if ( !strncmp(line,"Mt_rRNA",7) ) return GF_MT_rRNA;
            else if ( !strncmp(line,"Mt_tRNA",7) ) return GF_MT_tRNA;
        case 'l':
            if ( !strncmp(line,"lincRNA",7) ) return GF_lincRNA;
        case 'm':
            if ( !strncmp(line,"miRNA",5) ) return GF_miRNA;
            else if ( !strncmp(line,"misc_RNA",8) ) return GF_MISC_RNA;
            else if ( !strncmp(line,"macro_lncRNA",12) ) return GF_macro_lncRNA;
        case 'r':
            if ( !strncmp(line,"rRNA",4) ) return GF_rRNA;
            if ( !strncmp(line,"ribozyme",8) ) return GF_ribozyme;
        case 's':
            if ( !strncmp(line,"snRNA",5) ) return GF_snRNA;
            else if ( !strncmp(line,"sRNA",4) ) return GF_sRNA;
            else if ( !strncmp(line,"scRNA",5) ) return GF_scRNA;
            else if ( !strncmp(line,"scaRNA",6) ) return GF_scaRNA;
            else if ( !strncmp(line,"snoRNA",6) ) return GF_snoRNA;
            else if ( !strncmp(line,"sense_intronic",14) ) return GF_SENSE_INTRONIC;
            else if ( !strncmp(line,"sense_overlapping",17) ) return GF_SENSE_OVERLAPPING;
    }
    return 0;
}
static inline int gff_ignored_biotype(char *ss)
{
    ss = strstr(ss,"biotype=") + 8;
    if ( !strncmp("non_stop_decay",ss,14) ) return 1;
    if ( !strncmp("IG_pseudogene",ss,13) ) return 1;
    if ( !strncmp("IG_D_pseudogene",ss,15) ) return 1;
    if ( !strncmp("IG_V_pseudogene",ss,15) ) return 1;
    if ( !strncmp("TR_J_pseudogene",ss,15) ) return 1;
    if ( !strncmp("translated_processed_pseudogene",ss,31) ) return 1;
    if ( !strncmp("processed_pseudogene",ss,20) ) return 1;
    if ( !strncmp("unitary_pseudogene",ss,18) ) return 1;
    if ( !strncmp("pseudogene",ss,10) ) return 1;
    if ( !strncmp("retained_intron",ss,15) ) return 1;
    if ( !strncmp("IG_C_pseudogene",ss,15) ) return 1;
    if ( !strncmp("IG_J_pseudogene",ss,15) ) return 1;
    if ( !strncmp("nonsense_mediated_decay",ss,23) ) return 1;
    if ( !strncmp("unprocessed_pseudogene",ss,22) ) return 1;
    if ( !strncmp("transcribed_unprocessed_pseudogene",ss,34) ) return 1;
    if ( !strncmp("transcribed_processed_pseudogene",ss,32) ) return 1;
    if ( !strncmp("3prime_overlapping_ncrna",ss,24) ) return 1;
    if ( !strncmp("3prime_overlapping_ncRNA",ss,24) ) return 1;
    if ( !strncmp("TEC",ss,3) ) return 1;
    if ( !strncmp("bidirectional_promoter_lncRNA",ss,29) ) return 1;
    if ( !strncmp("transcribed_unitary_pseudogene",ss,30) ) return 1;
    return 0;
}
gf_gene_t *gene_init(aux_t *aux, uint32_t gene_id)
{
    khint_t k = kh_get(int2gene, aux->gid2gene, (int)gene_id);
    gf_gene_t *gene = (k == kh_end(aux->gid2gene)) ? NULL : kh_val(aux->gid2gene, k);
    if ( !gene )
    {
        gene = (gf_gene_t*) calloc(1,sizeof(gf_gene_t));
        int ret;
        k = kh_put(int2gene, aux->gid2gene, (int)gene_id, &ret);
        kh_val(aux->gid2gene,k) = gene;
    }
    return gene;
}
void gff_parse_transcript(args_t *args, const char *line, char *ss, ftr_t *ftr)
{
    aux_t *aux = &args->init;
    int biotype = gff_parse_biotype(ss);
    if ( biotype <= 0 )
    {
        if ( !gff_ignored_biotype(ss) ) fprintf(stderr,"ignored transcript: %s\n",line);
        return;
    }

    // create a mapping from transcript_id to gene_id
    uint32_t trid = gff_parse_id(line, "ID=transcript:", ss);
    uint32_t gene_id = gff_parse_id(line, "Parent=gene:", ss);

    if ( !ENSID_FMT ) gff_parse_ensid_fmt(line, "ID=transcript:", ss);      // id prefix different across species

    tscript_t *tr = (tscript_t*) calloc(1,sizeof(tscript_t));
    tr->id     = trid;
    tr->strand = ftr->strand;
    tr->gene   = gene_init(aux, gene_id);
    tr->type   = biotype;
    tr->beg    = ftr->beg;
    tr->end    = ftr->end;

    khint_t k;
    int ret;
    k = kh_put(int2tscript, aux->id2tr, (int)trid, &ret);
    kh_val(aux->id2tr,k) = tr;
}
void gff_parse_gene(args_t *args, const char *line, char *ss, char *chr_beg, char *chr_end, ftr_t *ftr)
{
    int biotype = gff_parse_biotype(ss);
    if ( biotype <= 0 )
    {
        if ( !gff_ignored_biotype(ss) ) fprintf(stderr,"ignored gene: %s\n",line);
        return;
    }

    aux_t *aux = &args->init;

    // substring search for "ID=gene:ENSG00000437963"
    uint32_t gene_id = gff_parse_id(line, "ID=gene:", ss);
    gf_gene_t *gene = gene_init(aux, gene_id);
    assert( !gene->name );      // the gene_id should be unique

    gene->iseq = feature_set_seq(args, chr_beg,chr_end);

    // substring search for "Name=OR4F5"
    ss = strstr(chr_end+2,"Name=");
    if ( !ss ) error("Could not parse the line, \"Name=\" not present: %s\n", line);
    ss += 5;
    char *se = ss;
    while ( *se && *se!=';' && !isspace(*se) ) se++;
    gene->name = (char*) malloc(se-ss+1);
    memcpy(gene->name,ss,se-ss);
    gene->name[se-ss] = 0;
}
int gff_parse(args_t *args, char *line, ftr_t *ftr)
{
    // - skip empty lines and commented lines
    // - columns 
    //      1.      chr
    //      2.      <skip>
    //      3.      CDS, transcript, gene, ...
    //      4-5.    beg,end
    //      6.      <skip>
    //      7.      strand
    //      8.      phase
    //      9.      Parent=transcript:ENST(\d+);ID=... etc

    char *ss = line;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments

    char *chr_beg, *chr_end;
    gff_parse_chr(line, &chr_beg, &chr_end);
    ss = gff_skip(line, chr_end + 2);

    // 3. column: is this a CDS, transcript, gene, etc.
    if ( !strncmp("exon\t",ss,5) ) { ftr->type = GF_EXON; ss += 5; }
    else if ( !strncmp("CDS\t",ss,4) ) { ftr->type = GF_CDS; ss += 4; }
    else if ( !strncmp("three_prime_UTR\t",ss,16) ) { ftr->type = GF_UTR3; ss += 16; }
    else if ( !strncmp("five_prime_UTR\t",ss,15) ) { ftr->type = GF_UTR5; ss += 15; }
    else
    {
        ss = gff_skip(line, ss);
        ss = gff_parse_beg_end(line, ss, &ftr->beg,&ftr->end);
        ss = gff_skip(line, ss);
        int type = gff_parse_type(ss);
        if ( type!=GFF_TSCRIPT_LINE && type!=GFF_GENE_LINE ) 
        {
            // we ignore these, debug print to see new types:
            ss = strstr(ss,"ID=");
            if ( !ss ) return -1;   // no ID, ignore the line
            if ( !strncmp("chromosome",ss+3,10) ) return -1;
            if ( !strncmp("supercontig",ss+3,11) ) return -1;
            fprintf(stderr,"ignore: %s\n", line);
            return -1;
        }

        // 7. column: strand
        if ( *ss == '+' ) ftr->strand = STRAND_FWD;
        else if ( *ss == '-' ) ftr->strand = STRAND_REV;
        else error("Unknown strand: %c .. %s\n", *ss,ss);

        if ( type==GFF_TSCRIPT_LINE )
            gff_parse_transcript(args, line, ss, ftr);
        else
            gff_parse_gene(args, line, ss, chr_beg, chr_end, ftr);

        return -1;
    }
    ss = gff_parse_beg_end(line, ss, &ftr->beg,&ftr->end);
    ss = gff_skip(line, ss);

    // 7. column: strand
    if ( *ss == '+' ) ftr->strand = STRAND_FWD;
    else if ( *ss == '-' ) ftr->strand = STRAND_REV;
    else { fprintf(stderr,"Skipping unknown strand: %c\n", *ss); return -1; }
    ss += 2;

    // 8. column: phase (codon offset)
    if ( *ss == '0' ) ftr->phase = 0;
    else if ( *ss == '1' ) ftr->phase = 1;
    else if ( *ss == '2' ) ftr->phase = 2;
    else if ( *ss == '.' ) ftr->phase = 0;      // exons do not have phase
    else { fprintf(stderr,"Skipping unknown phase: %c, %s\n", *ss, line); return -1; }
    ss += 2;

    // substring search for "Parent=transcript:ENST00000437963"
    ftr->trid = gff_parse_id(line, "Parent=transcript:", ss);
    ftr->iseq = feature_set_seq(args, chr_beg,chr_end);
    return 0;
}

static int cmp_cds_ptr(const void *a, const void *b)
{
    // comparison function for qsort of transcripts's CDS
    if ( (*((gf_cds_t**)a))->beg < (*((gf_cds_t**)b))->beg ) return -1;
    if ( (*((gf_cds_t**)a))->beg > (*((gf_cds_t**)b))->beg ) return 1;
    return 0;
}

static inline void chr_beg_end(aux_t *aux, int iseq, char **chr_beg, char **chr_end)
{
    *chr_beg = *chr_end = aux->seq[iseq];
    while ( (*chr_end)[1] ) (*chr_end)++;
}
tscript_t *tscript_init(aux_t *aux, uint32_t trid)
{
    khint_t k = kh_get(int2tscript, aux->id2tr, (int)trid);
    tscript_t *tr = (k == kh_end(aux->id2tr)) ? NULL : kh_val(aux->id2tr, k);
    assert( tr );
    return tr;
}
void register_cds(args_t *args, ftr_t *ftr)
{
    // Make the CDS searchable via idx_cds. Note we do not malloc tr->cds just yet.
    //  ftr is the result of parsing a gff CDS line
    aux_t *aux = &args->init;

    tscript_t *tr = tscript_init(aux, ftr->trid);
    if ( tr->strand != ftr->strand ) error("Conflicting strand in transcript %"PRIu32" .. %d vs %d\n",ftr->trid,tr->strand,ftr->strand);
    
    gf_cds_t *cds = (gf_cds_t*) malloc(sizeof(gf_cds_t));
    cds->tr    = tr;
    cds->beg   = ftr->beg;
    cds->len   = ftr->end - ftr->beg + 1;
    cds->icds  = 0;     // to keep valgrind on mac happy
    cds->phase = ftr->phase;
    
    hts_expand(gf_cds_t*,tr->ncds+1,tr->mcds,tr->cds);
    tr->cds[tr->ncds++] = cds;
}
void register_utr(args_t *args, ftr_t *ftr)
{
    aux_t *aux = &args->init;
    gf_utr_t *utr = (gf_utr_t*) malloc(sizeof(gf_utr_t));
    utr->which = ftr->type==GF_UTR3 ? prime3 : prime5;
    utr->beg   = ftr->beg;
    utr->end   = ftr->end;
    utr->tr    = tscript_init(aux, ftr->trid);

    char *chr_beg, *chr_end;
    chr_beg_end(&args->init, utr->tr->gene->iseq, &chr_beg, &chr_end);
    regidx_push(args->idx_utr, chr_beg,chr_end, utr->beg,utr->end, &utr);
}
void register_exon(args_t *args, ftr_t *ftr)
{
    aux_t *aux = &args->init;
    gf_exon_t *exon = (gf_exon_t*) malloc(sizeof(gf_exon_t));
    exon->beg = ftr->beg;
    exon->end = ftr->end;
    exon->tr  = tscript_init(aux, ftr->trid);

    char *chr_beg, *chr_end;
    chr_beg_end(&args->init, exon->tr->gene->iseq, &chr_beg, &chr_end);
    regidx_push(args->idx_exon, chr_beg,chr_end, exon->beg - N_SPLICE_REGION_INTRON, exon->end + N_SPLICE_REGION_INTRON, &exon);
}

void tscript_init_cds(args_t *args)
{
    aux_t *aux = &args->init;

    // Sort CDS in all transcripts, set offsets, check their phase, length, create index (idx_cds)
    khint_t k;
    for (k=0; k<kh_end(aux->id2tr); k++)
    {
        if ( !kh_exist(aux->id2tr, k) ) continue;
        tscript_t *tr = (tscript_t*) kh_val(aux->id2tr, k);

        // position-to-tscript lookup
        char *chr_beg, *chr_end;
        chr_beg_end(aux, tr->gene->iseq, &chr_beg, &chr_end);
        regidx_push(args->idx_tscript, chr_beg, chr_end, tr->beg, tr->end, &tr);

        if ( !tr->ncds ) continue;      // transcript with no CDS

        // sort CDs
        qsort(tr->cds, tr->ncds, sizeof(gf_cds_t*), cmp_cds_ptr);

        // trim non-coding start; NB: we do not care about reading frame, just the CDS
        int i, len = 0;
        if ( tr->strand==STRAND_FWD )
        {
            tr->cds[0]->beg += tr->cds[0]->phase;
            tr->cds[0]->len -= tr->cds[0]->phase;
            tr->cds[0]->phase = 0;

            // sanity check phase
            for (i=0; i<tr->ncds; i++)
            {
                int phase = tr->cds[i]->phase ? 3 - tr->cds[i]->phase : 0;
                if ( phase!=len%3)
                    error("GFF3 assumption failed for transcript %s, CDS=%d: phase!=len%%3 (phase=%d, len=%d)\n",ENSID(tr->id),tr->cds[i]->beg+1,phase,len);
                assert( phase == len%3 );
                len += tr->cds[i]->len; 
            }
        }
        else
        {
            // Check that the phase is not bigger than CDS length. Curiously, this can really happen,
            // see Mus_musculus.GRCm38.85.gff3.gz, transcript:ENSMUST00000163141
            // todo: the same for the fwd strand
            i = tr->ncds - 1;
            int phase = tr->cds[i]->phase;
            while ( i>=0 && phase > tr->cds[i]->len )
            {
                phase -= tr->cds[i]->len;
                tr->cds[i]->phase = 0;
                tr->cds[i]->len   = 0;
                i--;
            }
            tr->cds[i]->len  -= tr->cds[i]->phase;
            tr->cds[i]->phase = 0;

            // sanity check phase
            for (i=tr->ncds-1; i>=0; i--)
            {
                int phase = tr->cds[i]->phase ? 3 - tr->cds[i]->phase : 0;
                if ( phase!=len%3)
                    error("GFF3 assumption failed for transcript %s, CDS=%d: phase!=len%%3 (phase=%d, len=%d)\n",ENSID(tr->id),tr->cds[i]->beg+1,phase,len);
                len += tr->cds[i]->len;
            }
        }

        // set len. At the same check that CDS within a transcript do not overlap
        len = 0;
        for (i=0; i<tr->ncds; i++)
        {
            tr->cds[i]->icds = i;
            len += tr->cds[i]->len; 
            if ( !i ) continue;

            gf_cds_t *a = tr->cds[i-1];
            gf_cds_t *b = tr->cds[i];
            if ( a->beg + a->len - 1 >= b->beg ) 
                error("Error: CDS overlap in the transcript %"PRIu32": %"PRIu32"-%"PRIu32" and %"PRIu32"-%"PRIu32"\n", 
                    kh_key(aux->id2tr, k), a->beg+1,a->beg+a->len, b->beg+1,b->beg+b->len);
        }
        if ( len%3 != 0 )
        {
            // There are 13k transcripts with incomplete 3' CDS. See for example ENST00000524289
            //  http://sep2015.archive.ensembl.org/Homo_sapiens/Transcript/Sequence_cDNA?db=core;g=ENSG00000155868;r=5:157138846-157159019;t=ENST00000524289
            // Also, the incomplete CDS can be too short (1 or 2bp), so it is not enough to trim the last one.

            if ( tr->strand==STRAND_FWD )
            {
                i = tr->ncds - 1;
                while ( i>=0 && len%3 )
                {
                    int dlen = tr->cds[i]->len >= len%3 ? len%3 : tr->cds[i]->len;
                    tr->cds[i]->len -= dlen;
                    len -= dlen;
                    i--;
                }
            }
            else
            {
                i = 0;
                while ( i<tr->ncds && len%3 )
                {
                    int dlen = tr->cds[i]->len >= len%3 ? len%3 : tr->cds[i]->len;
                    tr->cds[i]->len -= dlen;
                    tr->cds[i]->beg += dlen;
                    len -= dlen;
                    i++;
                }
            }
        }

        // set CDS offsets and insert into regidx
        len=0;
        for (i=0; i<tr->ncds; i++)
        {
            tr->cds[i]->pos = len;
            len += tr->cds[i]->len;
            regidx_push(args->idx_cds, chr_beg,chr_end, tr->cds[i]->beg,tr->cds[i]->beg+tr->cds[i]->len-1, &tr->cds[i]);
        }
    }
}

void regidx_free_gf(void *payload) { free(*((gf_cds_t**)payload)); }
void regidx_free_tscript(void *payload) { tscript_t *tr = *((tscript_t**)payload); free(tr->cds); free(tr); }

void init_gff(args_t *args)
{
    aux_t *aux = &args->init;
    aux->seq2int   = khash_str2int_init();   // chrom's numeric id
    aux->gid2gene  = kh_init(int2gene);      // gene id to gf_gene_t, for idx_gene
    aux->id2tr     = kh_init(int2tscript);   // transcript id to tscript_t
    args->idx_tscript = regidx_init(NULL, NULL, regidx_free_tscript, sizeof(tscript_t*), NULL);

    // parse gff
    kstring_t str = {0,0,0};
    htsFile *fp = hts_open(args->gff_fname,"r");
    if ( !fp ) error("Failed to read %s\n", args->gff_fname);
    while ( hts_getline(fp, KS_SEP_LINE, &str) > 0 )
    {
        hts_expand(ftr_t, aux->nftr+1, aux->mftr, aux->ftr);
        int ret = gff_parse(args, str.s, aux->ftr + aux->nftr);
        if ( !ret ) aux->nftr++;
    }
    free(str.s);
    if ( hts_close(fp)!=0 ) error("Close failed: %s\n", args->gff_fname);


    // process gff information: connect CDS and exons to transcripts
    args->idx_cds  = regidx_init(NULL, NULL, regidx_free_gf, sizeof(gf_cds_t*), NULL);
    args->idx_utr  = regidx_init(NULL, NULL, regidx_free_gf, sizeof(gf_utr_t*), NULL);
    args->idx_exon = regidx_init(NULL, NULL, regidx_free_gf, sizeof(gf_exon_t*), NULL);
    args->itr      = regitr_init(NULL);

    int i;
    for (i=0; i<aux->nftr; i++)
    {
        ftr_t *ftr = &aux->ftr[i];

        // check whether to keep this feature: is there a mapping trid -> gene_id -> gene?
        khint_t k = kh_get(int2tscript, aux->id2tr, (int)ftr->trid);
        if ( k==kh_end(aux->id2tr) ) continue;       // no such transcript

        tscript_t *tr = kh_val(aux->id2tr,k);
        if ( !tr->gene->name )
        {
            // not a supported biotype (e.g. gene:pseudogene, transcript:processed_transcript)
            regidx_free_tscript(&tr);
            kh_del(int2tscript, aux->id2tr,k);
            continue;
        }

        // populate regidx by category: 
        //      ftr->type   .. GF_CDS, GF_EXON, GF_UTR3, GF_UTR5
        //      gene->type  .. GF_PROTEIN_CODING, GF_MT_rRNA, GF_IG_C, ...
        if ( ftr->type==GF_CDS ) register_cds(args, ftr);
        else if ( ftr->type==GF_EXON ) register_exon(args, ftr);
        else if ( ftr->type==GF_UTR5 ) register_utr(args, ftr);
        else if ( ftr->type==GF_UTR3 ) register_utr(args, ftr);
        else
            error("something: %s\t%d\t%d\t%s\t%s\n", aux->seq[ftr->iseq],ftr->beg+1,ftr->end+1,ENSID(ftr->trid),gf_type2gff_string(ftr->type));
    }
    tscript_init_cds(args);

    if ( !args->quiet )
    {
        fprintf(stderr,"Indexed %d transcripts, %d exons, %d CDSs, %d UTRs\n", 
                regidx_nregs(args->idx_tscript),
                regidx_nregs(args->idx_exon),
                regidx_nregs(args->idx_cds),
                regidx_nregs(args->idx_utr));
    }

    free(aux->ftr);
    khash_str2int_destroy_free(aux->seq2int);
    // keeping only to destroy the genes at the end: kh_destroy(int2gene,aux->gid2gene);
    kh_destroy(int2tscript,aux->id2tr);
    free(aux->seq);
}

void init_data(args_t *args)
{
    if ( !args->quiet ) fprintf(stderr,"Parsing %s ...\n", args->gff_fname);
    init_gff(args);

    args->rid = -1;

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);

    args->fai = fai_load(args->fa_fname);
    if ( !args->fai ) error("Failed to load the fai index: %s\n", args->fa_fname);

    args->pos2vbuf  = kh_init(pos2vbuf);
    args->active_tr = khp_init(trhp);
    args->hap = (hap_t*) calloc(1,sizeof(hap_t));

    // init samples
    if ( args->sample_list && !strcmp("-",args->sample_list) )
    {
        // ignore all samples
        if ( args->output_type==FT_TAB_TEXT ) 
        {
            // significant speedup for plain VCFs
            bcf_hdr_set_samples(args->hdr,NULL,0);
        }
        args->phase = PHASE_DROP_GT;
    }
    else
        args->smpl = smpl_ilist_init(args->hdr, args->sample_list, args->sample_is_file, SMPL_STRICT);
    args->hdr_nsmpl = args->phase==PHASE_DROP_GT ? 0 : bcf_hdr_nsamples(args->hdr);

    if ( args->output_type==FT_TAB_TEXT )
    {
        args->out = args->output_fname ? fopen(args->output_fname,"w") : stdout;
        if ( !args->out ) error("Failed to open %s: %s\n", args->output_fname,strerror(errno));

        fprintf(args->out,"# This file was produced by: bcftools +csq(%s+htslib-%s)\n", bcftools_version(),hts_version());
        fprintf(args->out,"# The command line was:\tbcftools +%s", args->argv[0]);
        int i;
        for (i=1; i<args->argc; i++)
            fprintf(args->out," %s",args->argv[i]);
        fprintf(args->out,"\n");
        fprintf(args->out,"# LOG\t[2]Message\n");
        fprintf(args->out,"# CSQ"); i = 1;
        fprintf(args->out,"\t[%d]Sample", ++i);
        fprintf(args->out,"\t[%d]Haplotype", ++i);
        fprintf(args->out,"\t[%d]Chromosome", ++i);
        fprintf(args->out,"\t[%d]Position", ++i);
        fprintf(args->out,"\t[%d]Consequence", ++i);
        fprintf(args->out,"\n");
    }
    else
    {
        args->out_fh = hts_open(args->output_fname? args->output_fname : "-",hts_bcf_wmode(args->output_type));
        if ( args->out_fh == NULL ) error("Can't write to %s: %s\n", args->output_fname? args->output_fname : "standard output", strerror(errno));
        bcf_hdr_append_version(args->hdr,args->argc,args->argv,"bcftools/csq");
        bcf_hdr_append(args->hdr,"##INFO=<ID=BCSQ,Number=.,Type=String,Description=\"Consequence annotation from bcftools/csq. Format: ..\">");
        bcf_hdr_append(args->hdr,"##FORMAT=<ID=BCSQ,Number=1,Type=Integer,Description=\"Bitmask of indexes to INFO/BCSQ, with interleaved first/second haplotype. Use \\\"bcftools query -f'[%CHROM\\t%POS\\t%SAMPLE\\t%TBCSQ\\n]'\\\" to translate.\">");
        bcf_hdr_write(args->out_fh, args->hdr);
    }
    //if ( !args->quiet ) fprintf(stderr,"Read %d coding genes, %d transcripts. Ignoring NMDs, pseudogenes, introns, etc.\n", ngene-1,ntr);
    if ( !args->quiet ) fprintf(stderr,"Calling...\n");
}

void destroy_data(args_t *args)
{
    regidx_destroy(args->idx_cds);
    regidx_destroy(args->idx_utr);
    regidx_destroy(args->idx_exon);
    regidx_destroy(args->idx_tscript);
    regitr_destroy(args->itr);

    khint_t k,i,j;
    for (k=0; k<kh_end(args->init.gid2gene); k++)
    {
        if ( !kh_exist(args->init.gid2gene, k) ) continue;
        gf_gene_t *gene = (gf_gene_t*) kh_val(args->init.gid2gene, k);
        free(gene->name);
        free(gene);
    }
    kh_destroy(int2gene,args->init.gid2gene);

    if ( args->filter )
        filter_destroy(args->filter);

    khp_destroy(trhp,args->active_tr);
    kh_destroy(pos2vbuf,args->pos2vbuf);
    if ( args->smpl ) smpl_ilist_destroy(args->smpl);
    int ret;
    if ( args->out_fh )
        ret = hts_close(args->out_fh);
    else
        ret = fclose(args->out);
    if ( ret ) error("Error: close failed .. %s\n", args->output_fname?args->output_fname:"stdout");
    for (i=0; i<args->vcf_rbuf.m; i++)
    {
        vbuf_t *vbuf = args->vcf_buf[i];
        if ( !vbuf ) continue;
        for (j=0; j<vbuf->m; j++)
        {
            if ( vbuf->vrec[j]->line ) bcf_destroy(vbuf->vrec[j]->line);
            free(vbuf->vrec[j]->smpl);
            free(vbuf->vrec[j]->csq);
            free(vbuf->vrec[j]);
        }
        free(vbuf->vrec);
        free(vbuf);
    }
    free(args->vcf_buf);
    free(args->rm_tr);
    free(args->csq_buf);
    free(args->hap->stack);
    free(args->hap->sseq.s);
    free(args->hap->tseq.s);
    free(args->hap->tref.s);
    free(args->hap);
    fai_destroy(args->fai);
    free(args->gt_arr);
    free(args->str.s);
    free(args->str2.s);
    free(ENSID_FMT);
}

int hap_init(hap_node_t *parent, hap_node_t *child, gf_cds_t *cds, bcf1_t *rec, int ial)
{
    int i;
    kstring_t str = {0,0,0};
    tscript_t *tr = cds->tr;
    uint32_t tr_beg = tr->cds[0]->beg;
    child->icds = cds->icds;     // index of cds in the tscript's list of exons

    // splice variants
    if ( rec->pos < cds->beg || rec->pos + rec->rlen > cds->beg + cds->len )
    {
        child->seq  = NULL;
        child->sbeg = 0;
        child->rbeg = rec->pos;
        child->rlen = 0;
        child->dlen = 0;
        kputs(rec->d.allele[0],&str);
        kputc('>',&str);
        kputs(rec->d.allele[ial],&str);
        child->var  = str.s;
        child->type = HAP_SSS;
        // NB: not checking if this is the first/last exon, the prediction may be incorrect
        if ( rec->pos < cds->beg )
            child->csq  = tr->strand==STRAND_REV ? CSQ_SPLICE_DONOR : CSQ_SPLICE_ACCEPTOR;
        else
            child->csq  = tr->strand==STRAND_REV ? CSQ_SPLICE_ACCEPTOR : CSQ_SPLICE_DONOR;
        child->prev = parent->type==HAP_SSS ? parent->prev : parent;
        child->rec  = rec;
        return 0;
    }

    if ( parent->type==HAP_SSS ) parent = parent->prev;
    if ( parent->type==HAP_CDS )    
    {
        i = parent->icds;
        if ( i!=cds->icds )
        {
            // the variant is on a new exon, finish up the previous
            int len = tr->cds[i]->len - parent->rbeg - parent->rlen + tr->cds[i]->beg;
            if ( len > 0 )
                kputsn_(tr->ref + parent->rbeg + parent->rlen - tr_beg, len, &str);
        }

        // append any skipped non-variant exons
        while ( ++i < cds->icds )
            kputsn_(tr->ref + tr->cds[i]->beg - tr_beg, tr->cds[i]->len, &str);

        if ( parent->icds==child->icds )
        {
            int len = rec->pos - parent->rbeg - parent->rlen;
            if ( len < 0 )   // overlapping variants
            {
                free(str.s);
                return 1;
            }
            kputsn_(tr->ref + parent->rbeg + parent->rlen - tr_beg, len, &str);
        }
        else
            kputsn_(tr->ref + cds->beg - tr_beg, rec->pos - cds->beg, &str);
    }
    kputs(rec->d.allele[ial], &str);

    child->seq  = str.s;
    child->sbeg = cds->pos + (rec->pos - cds->beg);
    child->rbeg = rec->pos;    
    child->rlen = rec->rlen;
    child->type = HAP_CDS;
    child->prev = parent;
    child->rec  = rec;

    // check for splice_region or start_lost
    if ( tr->strand==STRAND_FWD )
    {
        if ( cds->icds==0 && cds->beg - rec->pos < 3 ) child->csq = CSQ_START_LOST;
    }
    else
    {
        if ( child->icds==tr->ncds - 1 && cds->beg + cds->len - 3 <= rec->pos  ) child->csq = CSQ_START_LOST;
    }
    if ( child->csq==0 )
    {
        if ( cds->beg - rec->pos < 3 || cds->beg + cds->len - 3 <= rec->pos ) child->csq = CSQ_SPLICE_REGION;
    }


    // set vlen and the "ref>alt" string
    {
        int rlen = strlen(rec->d.allele[0]);
        int alen = strlen(rec->d.allele[ial]);
        child->dlen = alen - rlen;
        child->var  = (char*) malloc(rlen+alen+2);
        memcpy(child->var,rec->d.allele[0],rlen);
        child->var[rlen] = '>';
        memcpy(child->var+rlen+1,rec->d.allele[ial],alen);
        child->var[rlen+alen+1] = 0;
    }

if ( !(child->rbeg + child->rlen - cds->beg <= cds->len) )
{
    fprintf(stderr,"\n%d .. %d + %d - %d <= %d, %d:%s\n",tr->id,child->rbeg, child->rlen, cds->beg, cds->len,rec->pos+1,child->var);
    fprintf(stderr,"parent_icds=%d  child_icds=%d\n",parent->icds,child->icds);
}
assert( child->rbeg + child->rlen - cds->beg <= cds->len );    // todo: overlapping variant
    return 0;
}
void hap_destroy(hap_node_t *hap)
{
    int i;
    for (i=0; i<hap->nchild; i++)
        if ( hap->child[i] ) hap_destroy(hap->child[i]);
    for (i=0; i<hap->ncsq_list; i++)
        free(hap->csq_list[i].str);
    free(hap->csq_list);
    free(hap->child);
    free(hap->cur_child);
    free(hap->seq);
    free(hap->var);
    free(hap);
}


/*
    ref:    spliced reference and its length (ref.l)
    seq:    part of the spliced query transcript on the reference strand to translate, its 
                length (seq.l) and the total length of the complete transcript (seq.m)
    sbeg:   offset of seq within the spliced query transcript
    rbeg:   offset of seq within ref, 0-based
    rend:   last base of seq within ref, plus one. If seq does not contain indels, it is rend=rbeg+seq->l
    strand: coding strand - 0:rev, 1:fwd
    tseq:   translated sequence (aa)
    fill:   frameshift, fill until the end (strand=fwd) or from the start (strand=rev)
 */
void cds_translate(kstring_t *_ref, kstring_t *_seq, uint32_t sbeg, uint32_t rbeg, uint32_t rend, int strand, kstring_t *tseq, int fill)
{
    char tmp[3], *codon, *end;
    int i, len, npad;

    kstring_t ref = *_ref;
    kstring_t seq = *_seq;

#define DBG 0
#if DBG
 fprintf(stderr,"translate:\n");
 fprintf(stderr,"    ref: l=%d %s\n", (int)ref.l,ref.s);
 fprintf(stderr,"    seq: l=%d m=%d %s\n", (int)seq.l,(int)seq.m,seq.s);
 fprintf(stderr,"    sbeg,rbeg,rend: %d,%d,%d\n", sbeg,rbeg,rend);
 fprintf(stderr,"    strand,fill: %d,%d\n", strand,fill);
#endif

    tseq->l = 0;
    if ( strand==STRAND_FWD )
    {
        // left padding
        npad = sbeg % 3;
#if DBG>1
        fprintf(stderr,"    npad: %d\n",npad);
#endif
        assert( npad<=rbeg );

        for (i=0; i<npad; i++)
            tmp[i] = ref.s[rbeg+i-npad];
        for (; i<3 && i-npad<seq.l; i++)
            tmp[i] = seq.s[i-npad];
        len = seq.l - i + npad;    // the remaining length of padded sseq
#if DBG>1
        fprintf(stderr,"\t i=%d\n", i);
#endif
        if ( i==3 )
        {
            kputc_(dna2aa(tmp), tseq);
#if DBG>1
            fprintf(stderr,"[1]%c%c%c\n",tmp[0],tmp[1],tmp[2]);
#endif
            codon = seq.s + 3 - npad;        // next codon
            end   = codon + len - 1 - (len % 3);    // last position of a valid codon
            while ( codon < end )
            {
                kputc_(dna2aa(codon), tseq);
#if DBG>1
                fprintf(stderr,"[2]%c%c%c\n",codon[0],codon[1],codon[2]);
#endif
                codon += 3;
            }
            end = seq.s + seq.l - 1;
            for (i=0; codon+i<=end; i++) tmp[i] = codon[i];
        }

        // right padding
        codon = ref.s + rend;
        if ( i>0 )
        {
#if DBG>1
            if(i==1)fprintf(stderr,"[3]%c\n",tmp[0]);
            if(i==2)fprintf(stderr,"[3]%c%c\n",tmp[0],tmp[1]);
#endif
            for (; i<3; i++)
            {
                tmp[i] = *codon;
                codon++;
            }
            kputc_(dna2aa(tmp), tseq);
#if DBG>1
            fprintf(stderr,"[4]%c%c%c\n",tmp[0],tmp[1],tmp[2]);
#endif
        }
        if ( fill!=0 )
        {
            end = ref.s + ref.l;
            while ( codon+3 <= end )
            {
                kputc_(dna2aa(codon), tseq);
#if DBG>1
                fprintf(stderr,"[5]%c%c%c\t%c\n",codon[0],codon[1],codon[2],dna2aa(codon));
#endif
                codon += 3;
            }
        }
    }
    else
    {
        // right padding - number of bases to take from ref
        npad = (seq.m - (sbeg + seq.l)) % 3; 
#if DBG>1
        fprintf(stderr,"    npad: %d\n",npad);
#endif
        assert( npad>=0 && sbeg+seq.l+npad<=seq.m );  // todo: first codon on the rev strand

        if ( npad==2 )
        {
            tmp[1] = ref.s[rend];
            tmp[2] = ref.s[rend+1];
            i = 0;
        }
        else if ( npad==1 )
        {
            tmp[2] = ref.s[rend];
            i = 1;
        }
        else
            i = 2;

        end = seq.s + seq.l;
        for (; i>=0 && end>seq.s; i--) tmp[i] = *(--end);
#if DBG>1
        fprintf(stderr,"\t i=%d\n", i);
        if(i==1)fprintf(stderr,"[0]    %c\n",tmp[2]);
        if(i==0)fprintf(stderr,"[0]  %c%c\n",tmp[1],tmp[2]);
#endif
        if ( i==-1 )
        {
#if DBG>1
            fprintf(stderr,"[1]%c%c%c\t%c\n",tmp[0],tmp[1],tmp[2], cdna2aa(tmp));
#endif
            kputc_(cdna2aa(tmp), tseq);
            codon = end - 3;
            while ( codon >= seq.s )
            {
                kputc_(cdna2aa(codon), tseq);
#if DBG>1
                fprintf(stderr,"[2]%c%c%c\t%c\n",codon[0],codon[1],codon[2], cdna2aa(codon));
#endif
                codon -= 3;
            }
            if ( seq.s-codon==2 )
            {
                tmp[2] = seq.s[0]; 
                i = 1;
            }
            else if ( seq.s-codon==1 )
            {
                tmp[1] = seq.s[0]; 
                tmp[2] = seq.s[1];
                i = 0;
            }
            else
                i = -1;
#if DBG>1
            if(i==1)fprintf(stderr,"[3]   %c\n",tmp[2]);
            if(i==0)fprintf(stderr,"[3] %c%c\n",tmp[1],tmp[2]);
#endif
        }
        // left padding
        end = ref.s + rbeg;
        if ( i>=0 )
        {
            for (; i>=0 && end>=ref.s; i--) tmp[i] = *(--end);
            kputc_(cdna2aa(tmp), tseq);
#if DBG>1
            fprintf(stderr,"[4]%c%c%c\t%c\n",tmp[0],tmp[1],tmp[2],cdna2aa(tmp));
#endif
        }
        if ( fill!=0 )
        {
            codon = end - 3;
            while ( codon >= ref.s )
            {
                kputc_(cdna2aa(codon), tseq);
#if DBG>1
                fprintf(stderr,"[5]%c%c%c\t%c\n",codon[0],codon[1],codon[2],cdna2aa(codon));
#endif
                codon -= 3;
            }
        }
    }
    kputc_(0, tseq);
#if DBG
 fprintf(stderr,"    tseq: %s\n", tseq->s);
#endif
}

void tscript_splice_ref(tscript_t *tr)
{
    int i, len = 0;
    for (i=0; i<tr->ncds; i++) 
        len += tr->cds[i]->len;

    tr->nsref = len;
    tr->sref  = (char*) malloc(len+1);
    len = 0;

    for (i=0; i<tr->ncds; i++)
    {
        memcpy(tr->sref+len, tr->ref + tr->cds[i]->beg - tr->cds[0]->beg, tr->cds[i]->len);
        len += tr->cds[i]->len;
    }
    tr->sref[len] = 0;
}

int csq_push(args_t *args, csq_t *csq, bcf1_t *rec)
{
    khint_t k = kh_get(pos2vbuf, args->pos2vbuf, (int)csq->pos);
    vbuf_t *vbuf = (k == kh_end(args->pos2vbuf)) ? NULL : kh_val(args->pos2vbuf, k);
    if ( !vbuf ) error("This should not happen. %s:%d  %s\n",bcf_seqname(args->hdr,rec),csq->pos+1,csq->str);
    
    int i;
    for (i=0; i<vbuf->n; i++)
        if ( vbuf->vrec[i]->line==rec ) break;
    if ( i==vbuf->n ) error("This should not happen.. %s:%d  %s\n", bcf_seqname(args->hdr,rec),csq->pos+1,csq->str);
    vrec_t *vrec = vbuf->vrec[i];

    for (i=0; i<vrec->ncsq; i++)
        if ( !strcmp(csq->str,vrec->csq[i]) ) break;

    csq->vrec = vrec;
    csq->idx  = i;
    if ( i<vrec->ncsq ) return 1; // the consequence already exists

    vrec->ncsq++;
    hts_expand(char*, vrec->ncsq, vrec->mcsq, vrec->csq);
    vrec->csq[vrec->ncsq-1] = csq->str;

    return 0;
}

//  soff .. position of the variant within the trimmed query transcript
//  sbeg .. position of the variant within the query transcript
//  rbeg .. position on the reference transcript (if there are no indels, then rbeg=send)
//  rpos .. VCF position
#define node2soff(i) (hap->stack[i].slen - (hap->stack[i].node->rlen + hap->stack[i].node->dlen))
#define node2sbeg(i) (hap->sbeg + node2soff(i))
#define node2send(i) (hap->sbeg + hap->stack[i].slen)
#define node2rbeg(i) (hap->stack[i].node->sbeg)
#define node2rend(i) (hap->stack[i].node->sbeg + hap->stack[i].node->rlen)
#define node2rpos(i) (hap->stack[i].node->rbeg)

void kput_csq(uint32_t csq, kstring_t *str)
{
    int i, n = sizeof(csq_strings)/sizeof(char*);
    for (i=1; i<n; i++)
        if ( csq_strings[i] && csq&(1<<i) ) { kputs(csq_strings[i],str); break; }
    i++;
    for (; i<n; i++)
        if ( csq_strings[i] && csq&(1<<i) ) { kputc_('&',str); kputs(csq_strings[i],str); }
}

void hap_add_csq(args_t *args, hap_t *hap, hap_node_t *node, int tlen, int ibeg, int iend, int dlen)
{
    int i;
    tscript_t *tr = hap->tr;
    int ref_node = tr->strand==STRAND_FWD ? ibeg : iend;

    int inode = node->ncsq_list++;
    hts_expand(csq_t,node->ncsq_list,node->mcsq_list,node->csq_list);
    csq_t *csq = &node->csq_list[inode];
    csq->pos  = node2rpos(ref_node);
    csq->type = 0;
    for (i=ibeg; i<=iend; i++)
        csq->type |= hap->stack[i].node->csq;

    if ( hap->stack[ibeg].node->type != HAP_SSS  )
    {
        if ( dlen!=0 )
        {
            if ( dlen%3 )
                csq->type |= CSQ_FRAMESHIFT_VARIANT;
            else if ( dlen<0 )
                csq->type |= CSQ_INFRAME_DELETION;
            else
                csq->type |= CSQ_INFRAME_INSERTION;
        }
        else
        {
            for (i=0; i<hap->tref.l; i++) 
                if ( hap->tref.s[i] != hap->tseq.s[i] ) break;
            if ( i==hap->tref.l )
                csq->type |= CSQ_SYNONYMOUS_VARIANT;
            else if ( hap->tref.s[i] ==  '*' )
                csq->type |= CSQ_STOP_LOST;
            else if ( hap->tseq.s[i] ==  '*' )
                csq->type |= CSQ_STOP_GAINED;
            else
                csq->type |= CSQ_MISSENSE_VARIANT;
        }
    }

    kstring_t str = {0,0,0};

    // gene, transcript, strand
    if ( hap->upstream_stop ) kputc_('*', &str);
    kput_csq(csq->type, &str);
    ksprintf(&str, "|%s|%s|%c",tr->gene->name,ENSID(tr->id),tr->strand==STRAND_FWD ? '+' : '-');

    if ( hap->stack[ibeg].node->type == HAP_SSS  )
    {
        node->csq_list[inode].str = str.s;
        csq_push(args, node->csq_list+inode, hap->stack[ibeg].node->rec);
        return;
    }

    // create the aa variant string
    int aa_rbeg = tr->strand==STRAND_REV ? (hap->tr->nsref - node2rend(iend))/3+1 : node2rbeg(ibeg)/3+1;
    int aa_sbeg = tr->strand==STRAND_REV ? (tlen - node2send(iend))/3+1 : node2sbeg(ibeg)/3+1;
    for (i=0; i<hap->tref.l; i++)
        if ( hap->tref.s[i]=='*' ) break;
    if ( i!=hap->tref.l )
    {
        hap->tref.l = i+1;
        hap->tref.s[i+1] = 0;
    }
    for (i=0; i<hap->tseq.l; i++)
        if ( hap->tseq.s[i]=='*' ) break;
    if ( i!=hap->tseq.l )
    {
        hap->tseq.l = i+1;
        hap->tseq.s[i+1] = 0;
        hap->upstream_stop = 1;
    }
    kputc_('|', &str);
    kputw(aa_rbeg, &str);
    kputs(hap->tref.s, &str);
    kputc_('>', &str);
    kputw(aa_sbeg, &str);
    kputs(hap->tseq.s, &str);
    kputc_('|', &str);

    // create the dna variant string and, in case of combined variants,
    // insert silent CSQ_PRINTED_UPSTREAM variants
    for (i=ibeg; i<=iend; i++)
    {
        if ( i>ibeg ) kputc_('+', &str);
        kputw(node2rpos(i)+1, &str);
        kputs(hap->stack[i].node->var, &str);

        // csq are printed at one position only for combined variants, the rest is
        // silent and references the first
        if ( i!=ref_node )
        {
            kstring_t tmp_str = {0,0,0};
            kputc('@', &tmp_str);
            kputw(node2rpos(ref_node)+1, &tmp_str);

            node->ncsq_list++;
            hts_expand(csq_t,node->ncsq_list,node->mcsq_list,node->csq_list);
            csq_t *ctmp = &node->csq_list[node->ncsq_list-1];
            ctmp->type = CSQ_PRINTED_UPSTREAM;
            ctmp->pos  = node2rpos(i);
            ctmp->str  = tmp_str.s;

            csq_push(args, ctmp, hap->stack[i].node->rec);
        }
    }

    node->csq_list[inode].str = str.s;
assert(hap->stack[ibeg].node->rec);
    csq_push(args, node->csq_list+inode, hap->stack[ref_node].node->rec);
}

void hap_finalize(args_t *args, hap_t *hap)
{
    tscript_t *tr = hap->tr;
    if ( !tr->sref )
        tscript_splice_ref(tr);

    kstring_t sref;
    sref.s = tr->sref;
    sref.l = tr->nsref;
    sref.m = sref.l;

    int istack = 0;
    hts_expand(hstack_t,1,hap->mstack,hap->stack);

    hap->sseq.l = 0;
    hap->tseq.l = 0;
    hap->stack[0].node = tr->root;
    hap->stack[0].ichild = -1;
    hap->stack[0].slen = 0;
    hap->stack[0].dlen = 0;

    while ( istack>=0 )
    {
        hstack_t *stack  = &hap->stack[istack];
        hap_node_t *node = hap->stack[istack].node;
        while ( ++hap->stack[istack].ichild < node->nchild )
        {
            if ( node->child[stack->ichild] ) break;
        }
        if ( stack->ichild == node->nchild ) { istack--; continue; }

        node = node->child[stack->ichild];

        istack++;
        hts_expand(hstack_t,istack+1,hap->mstack,hap->stack);
        stack = &hap->stack[istack-1];

        hap->stack[istack].node = node;
        hap->stack[istack].ichild = -1;

        hap->sseq.l = stack->slen;
        if ( node->type==HAP_CDS ) kputs(node->seq, &hap->sseq);
        hap->stack[istack].slen = hap->sseq.l;
        hap->stack[istack].dlen = hap->stack[istack-1].dlen + node->dlen;

        if ( !node->nend ) continue;    // not a leaf node

        // The spliced sequence has been built for the current haplotype and stored
        // in hap->sseq. Now we break it and output as independent parts
        
        kstring_t sseq;
        sseq.m = sref.m + hap->stack[istack].dlen;  // total length of the spliced query transcript
        hap->upstream_stop = 0;

        int i = 1, dlen, ibeg;
        while ( i<istack && hap->stack[i].node->type == HAP_SSS ) i++;
        hap->sbeg = hap->stack[i].node->sbeg;

        if ( tr->strand==STRAND_FWD )
        {
            i = 0, dlen = 0, ibeg = -1;
            while ( ++i <= istack )
            {
                if ( hap->stack[i].node->type == HAP_SSS )
                {
                    // start/stop/splice site overlap: don't know how to build the haplotypes correctly, skipping
                    hap_add_csq(args,hap,node,0,i,i,0);
                    continue;
                }
                dlen += hap->stack[i].node->dlen;
                if ( i<istack )
                {
                    if ( dlen%3 )   // frameshift
                    {
                        if ( ibeg==-1 ) ibeg = i;
                        continue;
                    }
                    int icur  = node2sbeg(i);
                    int inext = node2sbeg(i+1);
                    if ( icur/3 == inext/3 )    // in the same codon, can't be flushed yet
                    {
                        if ( ibeg==-1 ) ibeg = i;
                        continue;
                    }
                }
                if ( ibeg<0 ) ibeg = i;
                int ioff = node2soff(ibeg);
                int icur = node2sbeg(ibeg);
                int rbeg = node2rbeg(ibeg);
                int rend = node2rend(i);

                sseq.l = hap->stack[i].slen - ioff;
                sseq.s = hap->sseq.s + ioff;
                cds_translate(&sref, &sseq, icur,rbeg,rend, tr->strand, &hap->tseq, dlen%3);

                sseq.s = sref.s + rbeg;
                sseq.l = node2rend(i) - rbeg;
                sseq.m = sref.m;
                cds_translate(&sref, &sseq, rbeg,rbeg,rend, tr->strand, &hap->tref, dlen%3);
                sseq.m = sref.m + hap->stack[istack].dlen;

                hap_add_csq(args,hap,node,0, ibeg,i,dlen);
                ibeg = -1;
                dlen = 0;
            }
        }
        else
        {
            i = istack + 1, dlen = 0, ibeg = -1;
            while ( --i > 0 )
            {
                if ( hap->stack[i].node->type == HAP_SSS )
                {
                    hap_add_csq(args,hap,node,0,i,i,0);
                    continue;
                }
                dlen += hap->stack[i].node->dlen;
                if ( i>1 && hap->stack[i-1].node->type != HAP_SSS )
                {
                    if ( dlen%3 )
                    {
                        if ( ibeg==-1 ) ibeg = i;
                        continue;
                    }
                    int icur  = sseq.m - 1 - node2sbeg(i);
                    int inext = sseq.m - 1 - node2sbeg(i-1);
                    if ( icur/3 == inext/3 )
                    {
                        if ( ibeg==-1 ) ibeg = i;
                        continue;
                    }
                }
                if ( ibeg<0 ) ibeg = i;
                int ioff = node2soff(i);
                int icur = node2sbeg(i);
                int rbeg = node2rbeg(i);
                int rend = node2rend(ibeg);

                sseq.l = hap->stack[ibeg].slen - ioff;
                sseq.s = hap->sseq.s + ioff;
                cds_translate(&sref, &sseq, icur,rbeg,rend, tr->strand, &hap->tseq, dlen%3);

                sseq.s = sref.s + rbeg;
                sseq.l = node2rend(ibeg) - rbeg;
                sseq.m = sref.m;
                cds_translate(&sref, &sseq, rbeg,rbeg,rend, tr->strand, &hap->tref, dlen%3);
                sseq.m = sref.m + hap->stack[istack].dlen;

                hap_add_csq(args,hap,node,sseq.m, i,ibeg,dlen);
                ibeg = -1;
                dlen = 0;
            }
        }
    }
}

static inline void csq_print_text(args_t *args, csq_t *csq, int ismpl, int ihap)
{
    if ( csq->type==CSQ_PRINTED_UPSTREAM ) return;

    char *smpl = ismpl >= 0 ? args->hdr->samples[ismpl] : "-";
    const char *chr = bcf_hdr_id2name(args->hdr,args->rid);

    fprintf(args->out,"CSQ\t%s\t", smpl);
    if ( ihap>0 )
        fprintf(args->out,"%d", ihap);
    else
        fprintf(args->out,"-");

    fprintf(args->out,"\t%s\t%d\t%s\n",chr,csq->pos+1,csq->str);
}
static inline void hap_print_text(args_t *args, tscript_t *tr, int ismpl, int ihap, hap_node_t *node)
{
    if ( !node || !node->ncsq_list ) return;

    char *smpl = ismpl >= 0 ? args->hdr->samples[ismpl] : "-";
    const char *chr = bcf_hdr_id2name(args->hdr,args->rid);

    int i;
    for (i=0; i<node->ncsq_list; i++)
    {
        csq_t *csq = node->csq_list + i;
        if ( csq->type==CSQ_PRINTED_UPSTREAM ) continue;
        assert( csq->str );

        fprintf(args->out,"CSQ\t%s\t", smpl);
        if ( ihap>0 )
            fprintf(args->out,"%d", ihap);
        else
            fprintf(args->out,"-");

        fprintf(args->out,"\t%s\t%d\t%s\n",chr,csq->pos+1,csq->str);
    }
}

static inline void hap_stage_vcf(args_t *args, tscript_t *tr, int ismpl, int ihap, hap_node_t *node)
{
    if ( !node || !node->ncsq_list || ismpl<0 ) return;

    int i;
    for (i=0; i<node->ncsq_list; i++)
    {
        csq_t *csq = node->csq_list + i;
        vrec_t *vrec = csq->vrec;
        // more than 16 consequences, ignoring the rest...
        if ( 2*csq->idx+ihap >= 32 ) break;
        vrec->smpl[ismpl] |= 1<<(2*csq->idx + ihap);
    }
}

void hap_flush(args_t *args, uint32_t pos)
{
    int i,j;
    tr_heap_t *heap = args->active_tr;

    while ( heap->ndat && heap->dat[0]->end<=pos )
    {
        tscript_t *tr = heap->dat[0];
        khp_delete(trhp, heap);

        args->hap->tr = tr;
        hap_finalize(args, args->hap);

        if ( args->output_type==FT_TAB_TEXT )   // plain text output, not a vcf
        {
            if ( args->phase==PHASE_DROP_GT )
                hap_print_text(args, tr, -1,0, tr->hap[0]);
            else
            {
                for (i=0; i<args->smpl->n; i++)
                {
                    for (j=0; j<2; j++)
                        hap_print_text(args, tr, args->smpl->idx[i],j+1, tr->hap[i*2+j]);
                }
            }
        }
        else if ( args->phase!=PHASE_DROP_GT )
        {
            for (i=0; i<args->smpl->n; i++)
            {
                for (j=0; j<2; j++)
                    hap_stage_vcf(args, tr, args->smpl->idx[i],j, tr->hap[i*2+j]);
            }
        }

        // mark the transcript for deletion. Cannot delete it immediately because
        // by-position VCF output will need them when flushed by vcf_buf_push
        args->nrm_tr++;
        hts_expand(tscript_t*,args->nrm_tr,args->mrm_tr,args->rm_tr);
        args->rm_tr[args->nrm_tr-1] = tr;
    }
}

#define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }

void vbuf_push(args_t *args, bcf1_t **rec_ptr)
{
    int i;

    assert(rec_ptr);
    bcf1_t *rec = *rec_ptr;

    // check for duplicate records
    i = args->vcf_rbuf.n ? rbuf_last(&args->vcf_rbuf) : -1;
    if ( i<0 || args->vcf_buf[i]->vrec[0]->line->pos!=rec->pos ) 
    {
        // vcf record with a new pos
        rbuf_expand0(&args->vcf_rbuf, vbuf_t*, args->vcf_rbuf.n+1, args->vcf_buf);
        i = rbuf_append(&args->vcf_rbuf);
        if ( !args->vcf_buf[i] ) args->vcf_buf[i] = (vbuf_t*) calloc(1,sizeof(vbuf_t));
        args->vcf_buf[i]->n = 0;
    }
    vbuf_t *vbuf = args->vcf_buf[i];
    vbuf->n++;
    hts_expand0(vrec_t*, vbuf->n, vbuf->m, vbuf->vrec);
    if ( !vbuf->vrec[vbuf->n - 1] )
        vbuf->vrec[vbuf->n - 1] = (vrec_t*) calloc(1,sizeof(vrec_t));

    vrec_t *vrec = vbuf->vrec[vbuf->n - 1];
    if ( args->phase!=PHASE_DROP_GT && args->smpl->n )
    {
        if ( !vrec->smpl ) vrec->smpl = (uint32_t*) calloc(args->hdr_nsmpl,sizeof(*vrec->smpl));
        else memset(vrec->smpl,0,args->hdr_nsmpl*sizeof(*vrec->smpl));
    }
    if ( !vrec->line ) vrec->line = bcf_init1();
    SWAP(bcf1_t*, (*rec_ptr), vrec->line);

    int ret;
    khint_t k = kh_put(pos2vbuf, args->pos2vbuf, (int)rec->pos, &ret);
    kh_val(args->pos2vbuf,k) = vbuf;

    return;
}

void vbuf_flush(args_t *args)
{
    if ( args->active_tr->ndat ) return; // cannot output buffered VCF lines (args.vbuf) until all active transcripts are gone

    int i,j;
    while ( (i=rbuf_shift(&args->vcf_rbuf))>=0 )
    {
        vbuf_t *vbuf = args->vcf_buf[i];
        for (i=0; i<vbuf->n; i++)
        {
            vrec_t *vrec = vbuf->vrec[i];
            if ( !args->out_fh ) // not a VCF output
            {
                vrec->ncsq = 0;
                continue;
            }
            if ( !vrec->ncsq )
            {
                bcf_write(args->out_fh, args->hdr, vrec->line);
                continue;
            }
            
            args->str.l = 0;
            kputs(vrec->csq[0], &args->str);
            for (j=1; j<vrec->ncsq; j++)
            {
                kputc_(',', &args->str);
                kputs(vrec->csq[j], &args->str);
            }
            bcf_update_info_string(args->hdr, vrec->line, "BCSQ", args->str.s);
            if ( args->hdr_nsmpl )
                bcf_update_format_int32(args->hdr, vrec->line, "BCSQ", vrec->smpl, args->hdr_nsmpl);
            vrec->ncsq = 0;

            bcf_write(args->out_fh, args->hdr, vrec->line);
        }
        vbuf->n = 0;
    }

    for (i=0; i<args->nrm_tr; i++)
    {
        tscript_t *tr = args->rm_tr[i];
        hap_destroy(tr->root);
        tr->root = NULL;
        free(tr->hap);
        free(tr->ref);
        free(tr->sref);
    }
    args->nrm_tr = 0;

    for (i=0; i<args->ncsq_buf; i++) free(args->csq_buf[i].str);
    args->ncsq_buf = 0;
}

int test_cds(args_t *args, bcf1_t *rec)
{
    int i;
    const char *chr = bcf_seqname(args->hdr,rec);
    if ( !regidx_overlap(args->idx_cds,chr,rec->pos,rec->pos+rec->rlen-1, args->itr) ) return 0;
    while ( regitr_overlap(args->itr) )
    {
        gf_cds_t *cds = regitr_payload(args->itr,gf_cds_t*);
        tscript_t *tr = cds->tr;

        if ( !tr->root )
        {
            // initialize the transcript and its haplotype tree, fetch the reference sequence
            int len, cds_end = tr->cds[tr->ncds-1]->beg + tr->cds[tr->ncds-1]->len - 1;
            tr->ref = faidx_fetch_seq(args->fai, chr, tr->cds[0]->beg, cds_end, &len);
            if ( !tr->ref )
                error("faidx_fetch_seq failed %s:%d-%d\n", chr,tr->cds[0]->beg+1,cds_end+1);
            if ( len != cds_end - tr->cds[0]->beg + 1 )
                error("incorrect length from faidx_fetch_seq %s:%d-%d .. %d\n", chr,tr->cds[0]->beg+1,cds_end+1,len);

            tr->root = (hap_node_t*) calloc(1,sizeof(hap_node_t));
            tr->nhap = args->phase==PHASE_DROP_GT ? 1 : 2*args->smpl->n;     // maximum ploidy = diploid
            tr->hap  = (hap_node_t**) malloc(tr->nhap*sizeof(hap_node_t*));
            for (i=0; i<tr->nhap; i++) tr->hap[i] = NULL;
            tr->root->nend = tr->nhap;
            tr->root->type = HAP_ROOT;

            khp_insert(trhp, args->active_tr, &tr);
        }

        if ( args->phase==PHASE_DROP_GT )
        {
            if ( rec->d.allele[1][0]=='<' || rec->d.allele[1][0]=='*' ) { continue; }
            hap_node_t *parent = tr->hap[0] ? tr->hap[0] : tr->root;
            hap_node_t *child  = (hap_node_t*)calloc(1,sizeof(hap_node_t));
            if ( hap_init(parent, child, cds, rec, 1)!=0 )
            {
                // overlapping variant, cannot apply
                if ( !args->quiet )
                    fprintf(stderr,"Warning: Skipping overlapping variants at %s:%d\t%s>%s\n", chr,rec->pos+1,rec->d.allele[0],rec->d.allele[1]);
                if ( args->out ) 
                    fprintf(args->out,"LOG\tWarning: Skipping overlapping variants at %s:%d\t%s>%s\n", chr,rec->pos+1,rec->d.allele[0],rec->d.allele[1]);
                free(child);
                continue;
            }
            parent->nend--;
            parent->nchild = 1;
            parent->mchild = 1;
            parent->child  = (hap_node_t**) malloc(sizeof(hap_node_t*));
            parent->child[0] = child;
            tr->hap[0] = child;
            tr->hap[0]->nend = 1;
            continue;
        }

        // apply the VCF variants and extend the haplotype tree
        int j, ismpl, ihap, ngts = bcf_get_genotypes(args->hdr, rec, &args->gt_arr, &args->mgt_arr);
        ngts /= bcf_hdr_nsamples(args->hdr);
        for (ismpl=0; ismpl<args->smpl->n; ismpl++)
        {
            int32_t *gt = args->gt_arr + args->smpl->idx[ismpl]*ngts;
            if ( gt[0]==bcf_gt_missing ) continue;

            if ( ngts>1 && gt[0]!=gt[1] && gt[1]!=bcf_gt_missing && gt[1]!=bcf_int32_vector_end )
            {
                if ( args->phase==PHASE_MERGE )
                {
                    if ( !bcf_gt_allele(gt[0]) ) gt[0] = gt[1];
                }
                else if ( !bcf_gt_is_phased(gt[0]) && !bcf_gt_is_phased(gt[1]) )
                {
                    if ( args->phase==PHASE_REQUIRE )
                        error("Unphased genotype at %s:%d, sample %s. See the --phase option.\n", chr,rec->pos+1,args->hdr->samples[args->smpl->idx[ismpl]]);
                    if ( args->phase==PHASE_SKIP )
                        continue;
                }
            }

            for (ihap=0; ihap<ngts; ihap++)
            {
                if ( gt[ihap]==bcf_gt_missing || gt[ihap]==bcf_int32_vector_end ) continue;

                i = 2*ismpl + ihap;

                int ial = bcf_gt_allele(gt[ihap]);
                if ( !ial ) continue;
                assert( ial < rec->n_allele );
                if ( rec->d.allele[ial][0]=='<' || rec->d.allele[ial][0]=='*' ) { continue; }

                hap_node_t *parent = tr->hap[i] ? tr->hap[i] : tr->root;
                if ( parent->cur_rec==rec && parent->cur_child[ial]>=0 )
                {
                    // this haplotype has been seen in another sample
                    tr->hap[i] = parent->child[ parent->cur_child[ial] ];
                    tr->hap[i]->nend++;
                    parent->nend--;
                    continue;
                }

                hap_node_t *child = (hap_node_t*)calloc(1,sizeof(hap_node_t));
                if ( hap_init(parent, child, cds, rec, ial)!=0 )
                {
                    // overlapping variant, cannot apply
                    if ( !args->quiet )
                        fprintf(stderr,"Warning: Skipping overlapping variants at %s:%d, sample %s\t%s>%s\n",
                                chr,rec->pos+1,args->hdr->samples[args->smpl->idx[ismpl]],rec->d.allele[0],rec->d.allele[ial]);
                    if ( args->out  )
                        fprintf(args->out,"LOG\tWarning: Skipping overlapping variants at %s:%d, sample %s\t%s>%s\n",
                                chr,rec->pos+1,args->hdr->samples[args->smpl->idx[ismpl]],rec->d.allele[0],rec->d.allele[ial]);

                    free(child);
                    continue;
                }

                if ( parent->cur_rec!=rec )
                {
                    hts_expand(int,rec->n_allele,parent->mcur_child,parent->cur_child);
                    for (j=0; j<rec->n_allele; j++) parent->cur_child[j] = -1;
                    parent->cur_rec = rec;
                }

                j = parent->nchild++;
                hts_expand0(hap_node_t*,parent->nchild,parent->mchild,parent->child);
                parent->cur_child[ial] = j;
                parent->child[j] = child;
                tr->hap[i] = child;
                tr->hap[i]->nend++;
                parent->nend--;
            }
        }
    }
    return 1;
}

void csq_stage(args_t *args, csq_t *csq, bcf1_t *rec)
{
    // known issues: tab output leads to unsorted output. This is because
    // coding haplotypes are printed in one go and buffering is not used
    // with tab output. VCF output is OK though.

    if ( csq_push(args, csq, rec)!=0 ) return;    // the consequence already exists

    int i,j,ngt = 0;
    if ( args->phase!=PHASE_DROP_GT )
    {
        ngt = bcf_get_genotypes(args->hdr, rec, &args->gt_arr, &args->mgt_arr);
        if ( ngt>0 ) ngt /= bcf_hdr_nsamples(args->hdr);
    }
    if ( ngt<=0 )
    {
        if ( args->output_type==FT_TAB_TEXT )
            csq_print_text(args, csq, -1,0);
        return;
    }
    assert( ngt<=2 );

    if ( args->output_type==FT_TAB_TEXT )
    {
        for (i=0; i<args->smpl->n; i++)
        {
            int32_t *gt = args->gt_arr + args->smpl->idx[i]*ngt;
            for (j=0; j<ngt; j++)
            {
                if ( gt[j]==bcf_gt_missing || gt[j]==bcf_int32_vector_end || !bcf_gt_allele(gt[j]) ) continue;
                csq_print_text(args, csq, args->smpl->idx[i],j+1);
            }
        }
        return;
    }

    vrec_t *vrec = csq->vrec;
    for (i=0; i<args->smpl->n; i++)
    {
        int32_t *gt = args->gt_arr + args->smpl->idx[i]*ngt;
        for (j=0; j<ngt; j++)
        {
            if ( gt[j]==bcf_gt_missing || gt[j]==bcf_int32_vector_end || !bcf_gt_allele(gt[j]) ) continue;

            // more than 16 consequences, ignoring the rest...
            if ( 2*csq->idx+j >= 32 ) break;
            vrec->smpl[i] |= 1<<(2*csq->idx + j);
        }
    }
}
int test_utr(args_t *args, bcf1_t *rec)
{
    const char *chr = bcf_seqname(args->hdr,rec);
    if ( !regidx_overlap(args->idx_utr,chr,rec->pos,rec->pos+rec->rlen-1, args->itr) ) return 0;
    while ( regitr_overlap(args->itr) )
    {
        gf_utr_t *utr = regitr_payload(args->itr, gf_utr_t*);
        tscript_t *tr = utr->tr;

        hts_expand(csq_t, args->ncsq_buf+1, args->mcsq_buf, args->csq_buf);
        csq_t *csq = &args->csq_buf[args->ncsq_buf++];
        csq->pos  = rec->pos;
        csq->type = utr->which==prime5 ? CSQ_UTR5 : CSQ_UTR3;

        args->str2.l = 0;
        kput_csq(csq->type, &args->str2);
        ksprintf(&args->str2, "|%s|%s|%c",tr->gene->name,ENSID(tr->id),tr->strand==STRAND_FWD ? '+' : '-');
        csq->str = strdup(args->str2.s);

        csq_stage(args, csq, rec);
    }
    return 1;
}
int test_splice(args_t *args, bcf1_t *rec)
{
    const char *chr = bcf_seqname(args->hdr,rec);

    int rec_end = rec->pos + rec->rlen - 1;
    if ( !regidx_overlap(args->idx_exon,chr,rec->pos,rec_end, args->itr) ) return 0;

    int ret = 0;
    while ( regitr_overlap(args->itr) )
    {
        gf_exon_t *exon = regitr_payload(args->itr, gf_exon_t*);
        tscript_t *tr = exon->tr;
        if ( !tr->ncds ) continue;  // not a coding transcript, no interest in splice sites

        // check whether the record falls inside (dpos<0) or outside (dpos>0) the exon
        int dpos = exon->beg - rec->pos, side = exon->beg;
        if ( abs(dpos) > abs(exon->beg - rec_end) ) dpos = exon->beg - rec_end, side = exon->beg;
        if ( abs(dpos) > abs(exon->end - rec->pos) ) dpos = rec->pos - exon->end, side = exon->end;
        if ( abs(dpos) > abs(exon->end - rec_end) ) dpos = rec_end - exon->end, side = exon->end;

        // skip the variant if it is at begining of the first exon or at the end of the last exon
        if ( side==tr->beg || side==tr->end ) continue;

        int type = 0;
        if ( dpos < 0 )   // inside the exon
        {
            if ( abs(dpos) < N_SPLICE_REGION_EXON ) type = CSQ_SPLICE_REGION;
        }
        else if ( dpos <= N_SPLICE_DONOR )
        {
            if ( tr->strand==STRAND_FWD )
                type = side==exon->beg ? CSQ_SPLICE_ACCEPTOR : CSQ_SPLICE_DONOR;
            else
                type = side==exon->beg ? CSQ_SPLICE_DONOR : CSQ_SPLICE_ACCEPTOR;
        }
        else if ( dpos <= N_SPLICE_REGION_INTRON ) type = CSQ_SPLICE_REGION;

        if ( !type ) continue;
        ret = 1;

        hts_expand(csq_t, args->ncsq_buf+1, args->mcsq_buf, args->csq_buf);
        csq_t *csq = &args->csq_buf[args->ncsq_buf++];
        csq->pos  = rec->pos;
        csq->type = type;

        args->str2.l = 0;
        kput_csq(csq->type, &args->str2);
        ksprintf(&args->str2, "|%s|%s|%c",tr->gene->name,ENSID(tr->id),tr->strand==STRAND_FWD ? '+' : '-');
        csq->str = strdup(args->str2.s);

        csq_stage(args, csq, rec);
    }
    return ret;
}
int test_tscript(args_t *args, bcf1_t *rec)
{
    const char *chr = bcf_seqname(args->hdr,rec);
    if ( !regidx_overlap(args->idx_tscript,chr,rec->pos,rec->pos+rec->rlen-1, args->itr) ) return 0;
    while ( regitr_overlap(args->itr) )
    {
        tscript_t *tr = regitr_payload(args->itr, tscript_t*);

        hts_expand(csq_t, args->ncsq_buf+1, args->mcsq_buf, args->csq_buf);
        csq_t *csq = &args->csq_buf[args->ncsq_buf++];
        csq->pos  = rec->pos;
        csq->type = GF_is_coding(tr->type) ? CSQ_INTRON : CSQ_NON_CODING;

        args->str2.l = 0;
        kput_csq(csq->type, &args->str2);
        kputc_('|', &args->str2);
        kputs(tr->gene->name, &args->str2);
        if ( csq->type==CSQ_NON_CODING )
        {
            kputc_('|', &args->str2);
            kputs(gf_type2gff_string(tr->type), &args->str2);
        }
        csq->str = strdup(args->str2.s);

        csq_stage(args, csq, rec);
    }
    return 1;
}

void process(args_t *args, bcf1_t **rec_ptr)
{
    if ( !rec_ptr )
    {
        hap_flush(args, REGIDX_MAX);
        vbuf_flush(args);
        return;
    }

    bcf1_t *rec = *rec_ptr;

    int call_csq = 1;
    if ( !rec->n_allele ) call_csq = 0;   // no alternate allele
    else if ( rec->n_allele==2 && (rec->d.allele[1][0]=='<' || rec->d.allele[1][0]=='*') ) call_csq = 0;     // gVCF, no alt allele
    else if ( args->filter )
    {
        call_csq = filter_test(args->filter, rec, NULL);
        if ( args->filter_logic==FLT_EXCLUDE ) call_csq = call_csq ? 0 : 1;
    }
    if ( !call_csq )
    {
        if ( !args->out_fh ) return;    // not a VCF output
        vbuf_push(args, rec_ptr);
        vbuf_flush(args);
        return;
    }

    if ( args->rid != rec->rid ) 
    {
        hap_flush(args, REGIDX_MAX);
        vbuf_flush(args);
    }
    args->rid = rec->rid;
    vbuf_push(args, rec_ptr);

    int hit = test_cds(args, rec);
    if ( !hit ) hit = test_utr(args, rec);
    if ( !hit ) hit = test_splice(args, rec);
    if ( !hit ) hit = test_tscript(args, rec);

    hap_flush(args, rec->pos-1);
    vbuf_flush(args);

    return;
}

const char *usage(void)
{
    return 
        "\n"
        "About: Haplotype-aware consequence caller. GFF3 annotation file can be downloaded e.g. from\n"
        "           ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/\n"
        "           ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/\n"
        "Usage: bcftools csq [options] in.vcf\n"
        "\n"
        "Options:\n"
        "   -e, --exclude <expr>            exclude sites for which the expression is true\n"
        "   -f, --fasta-ref <file>          reference file in fasta format\n"
        "   -g, --gff-annot <file>          gff3 annotation file\n"
        "   -i, --include <expr>            select sites for which the expression is true\n"
        "   -o, --output <file>             write output to a file [standard output]\n"
        "   -O, --output-type <b|u|z|v|t>   t:plain tab-delimited text output; b,u,z,v: as in \"bcftools view\"\n"
        "   -p, --phase <a|m|r|s>           how to construct haplotypes and how to deal with unphased data: [r]\n"
        "                                     a: take GTs as is, create haplotypes regardless of phase (0/1 -> 0|1)\n"
        "                                     m: merge all GTs into a single haplotype (0/1 -> 1, 1/2 -> 1)\n"
        "                                     r: require phased GTs, throw an error on unphased het GTs\n"
        "                                     s: skip unphased GTs\n"
        "   -q, --quiet                     suppress warning messages\n"
        "   -r, --regions <region>          restrict to comma-separated list of regions\n"
        "   -R, --regions-file <file>       restrict to regions listed in a file\n"
        "   -s, --samples <-|list>          samples to include or \"-\" to apply all variants and ignore samples\n"
        "   -S, --samples-file <file>       samples to include\n"
        "   -t, --targets <region>          similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file <file>       similar to -R but streams rather than index-jumps\n"
        "\n"
        "Example:\n"
        "   bcftools csq -f hs37d5.fa -g Homo_sapiens.GRCh37.82.gff3.gz in.vcf\n"
        "\n";
}

int main_csq(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc = argc; args->argv = argv;
    args->output_type = FT_VCF;

    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"gff-annot",1,0,'g'},
        {"fasta-ref",1,0,'f'},
        {"include",1,0,'i'},
        {"exclude",1,0,'e'},
        {"output",1,0,'o'},
        {"output-type",1,NULL,'O'},
        {"phase",1,0,'p'},
        {"quiet",0,0,'q'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {0,0,0,0}
    };
    int c, targets_is_file = 0, regions_is_file = 0; 
    char *targets_list = NULL, *regions_list = NULL;
    while ((c = getopt_long(argc, argv, "?hr:R:t:T:i:e:f:o:O:g:s:S:p:q",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'q': args->quiet = 1; break;
            case 'p':
                switch (optarg[0]) 
                {
                    case 'a': args->phase = PHASE_AS_IS; break;
                    case 'm': args->phase = PHASE_MERGE; break;
                    case 'r': args->phase = PHASE_REQUIRE; break;
                    case 's': args->phase = PHASE_SKIP; break;
                    default: error("The -p code \"%s\" not recognised\n", optarg);
                }
                break;
            case 'f': args->fa_fname = optarg; break;
            case 'g': args->gff_fname = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 't': args->output_type = FT_TAB_TEXT; break;
                          case 'b': args->output_type = FT_BCF_GZ; break;
                          case 'u': args->output_type = FT_BCF; break;
                          case 'z': args->output_type = FT_VCF_GZ; break;
                          case 'v': args->output_type = FT_VCF; break;
                          default: error("The output type \"%s\" not recognised\n", optarg);
                      }
                      break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'r': regions_list = optarg; break;
            case 'R': regions_list = optarg; regions_is_file = 1; break;
            case 's': args->sample_list = optarg; break;
            case 'S': args->sample_list = optarg; args->sample_is_file = 1; break;
            case 't': targets_list = optarg; break;
            case 'T': targets_list = optarg; targets_is_file = 1; break;
            case 'h':
            case '?': error("%s",usage());
            default: error("The option not recognised: %s\n\n", optarg); break;
        }
    }
    char *fname = NULL;
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else error("%s", usage());
    }
    else fname = argv[optind];
    if ( argc - optind>1 ) error("%s", usage());
    if ( !args->fa_fname ) error("Missing the --fa-ref option\n");
    if ( !args->gff_fname ) error("Missing the --gff option\n");
    args->sr = bcf_sr_init();
    if ( targets_list && bcf_sr_set_targets(args->sr, targets_list, targets_is_file, 0)<0 )
        error("Failed to read the targets: %s\n", targets_list);
    if ( regions_list && bcf_sr_set_regions(args->sr, regions_list, regions_is_file)<0 )
        error("Failed to read the regions: %s\n", regions_list);
    if ( !bcf_sr_add_reader(args->sr, fname) )
        error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->sr->errnum));
    args->hdr = bcf_sr_get_header(args->sr,0);

    init_data(args);
    while ( bcf_sr_next_line(args->sr) )
    {
        process(args, &args->sr->readers[0].buffer[0]);
    }
    process(args,NULL);

    destroy_data(args);
    bcf_sr_destroy(args->sr);
    free(args);

    return 0;
}

