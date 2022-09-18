/* The MIT License

   Copyright (c) 2016-2022 Genome Research Ltd.

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
    Illumina TOP/BOT strand convention causes a lot of pain. This tool
    attempts to determine the strand convention and convert it to the
    forward reference strand.

    On TOP strand, we can encounter
        unambiguous SNPs:
            A/G
            A/C
        ambiguous (context-dependent) SNPs:
            C/G
            A/T

    On BOT strand:
        unambiguous SNPs:
            T/G
            T/C
        ambiguous (context-dependent) SNPs:
            T/A
            G/C


    For unambiguous pairs (A/C, A/G, T/C, T/G), the knowledge of reference base
    at the SNP position is enough to determine the strand:

         TOP      REF   ->  ALLELES   TOP_ON_STRAND
         -------------------------------------------
         A/C     A or C      A/C         1
          "      T or G      T/G        -1
         A/G     A or G      A/G         1
          "      T or C      T/C        -1


    For ambiguous pairs (A/T, C/G), a sequence walking must be performed
    (simultaneously upstream and downstream) until the first unambiguous pair
    is encountered. The 5' base determines the strand:

         TOP    5'REF_BASE   ->  ALLELES   TOP_ON_STRAND
         ------------------------------------------------
         A/T    A or T             A/T          1
          "     C or G             T/A         -1
         C/G    A or T             C/G          1
          "     C or G             G/C         -1

 */


#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/kfunc.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/synced_bcf_reader.h>
#include "bcftools.h"

#define MODE_STATS          1
#define MODE_TOP2FWD        2
#define MODE_FLIP2FWD       3
#define MODE_USE_ID         4
#define MODE_REF_ALT        5
#define MODE_FLIP_ALL       6
#define MODE_SWAP_REF_ALT   7

#define FIX_ERR  (1<<0)
#define FIX_SKIP (1<<1)
#define FIX_NONE (1<<2)
#define FIX_FLIP (1<<3)
#define FIX_SWAP (1<<4)
#define FIX_GT   (1<<5)
const char *info_annots[] = {"err","skip","none","flip","swap","GT"};

typedef struct
{
    uint32_t pos;
    uint8_t ref;
}
marker_t;

KHASH_MAP_INIT_STR(i2m, marker_t)
typedef khash_t(i2m) i2m_t;

typedef struct
{
    char *dbsnp_fname;
    int mode, discard;
    bcf_hdr_t *hdr_in, *hdr_out;
    faidx_t *fai;
    int rid, skip_rid;
    i2m_t *i2m;
    int32_t *gts, ngts, pos;
    uint32_t nsite,nok,nflip,nunresolved,nswap,nflip_swap,nonSNP,nonACGT,nonbiallelic,nerr;
    uint32_t count[4][4], npos_err, unsorted;
    uint32_t dirty;
    kstring_t str;
    char *info_tag;
}
args_t;

args_t args;

const char *about(void)
{
    return "Fix reference strand orientation, e.g. from Illumina/TOP to fwd.\n";
}

const char *usage(void)
{
    return
        "\n"
        "About: This tool helps to determine and fix strand orientation.\n"
        "       Currently the following modes are recognised:\n"
        "           flip     .. swap or flip REF/ALT columns and GTs for non-ambiguous SNPs and ignore the rest\n"
        "           flip-all .. swap or flip REF/ALT columns and GTs for all SNPs, including ambiguous (A/T, C/G) sites\n"
        "           id       .. swap or flip REF/ALT columns and GTs using the ID column to determine the REF allele\n"
        "           ref-alt  .. swap or flip REF/ALT columns to match the reference, do not modify the genotypes\n"
        "           stats    .. collect and print stats\n"
        "           swap     .. swap REF/ALT columns to match the reference, do not modify the genotypes\n"
        "           top      .. convert from Illumina TOP strand to fwd\n"
        "       Each record is annotated with the INFO/FIXREF tag which records what change was made:\n"
        "           err      .. the alleles are invalid and could not be repaired\n"
        "           skip     .. the record was not considered at all (e.g. multiallelic or indel)\n"
        "           none     .. the alleles were valid and no change was required\n"
        "           flip     .. the REF and ALT alleles were flipped (e.g. A,C -> T,G)\n"
        "           swap     .. the REF and ALT columns were swapped (e.g. A,C -> C,A)\n"
        "           GT       .. FORMAT/GT was modified to reflect REF,ALT swap (e.g. 1/1 -> 0/0)\n"
        "\n"
        "       WARNING: Do not use the program blindly, make an effort to\n"
        "       understand what strand convention your data uses! Make sure\n"
        "       the reason for mismatching REF alleles is not a different\n"
        "       reference build!!\n"
        "\n"
        "       Please check this page before messing up your VCF even more\n"
        "           http://samtools.github.io/bcftools/howtos/plugin.fixref.html\n"
        "\n"
        "Usage: bcftools +fixref [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -d, --discard             Discard sites which could not be resolved\n"
        "   -f, --fasta-ref FILE.fa   Reference sequence\n"
        "   -i, --use-id FILE.vcf     Swap REF/ALT using the ID column to determine the REF allele, implies -m id.\n"
        "                             Download the dbSNP file from\n"
        "                               https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf\n"
        "   -m, --mode STRING         Collect statistics (\"stats\") or convert (see the list of modes above) [stats]\n"
        "   -t, --tag-name STRING     The name of the INFO annotation to record what change was made [FIXREF]\n"
        "\n"
        "Examples:\n"
        "   # run stats\n"
        "   bcftools +fixref file.bcf -- -f ref.fa\n"
        "\n"
        "   # convert from TOP to fwd\n"
        "   bcftools +fixref file.bcf -Ob -o out.bcf -- -f ref.fa -m top\n"
        "\n"
        "   # match the REF/ALT alleles based on the ID column, discard unknown sites\n"
        "   bcftools +fixref file.bcf -Ob -o out.bcf -- -d -f ref.fa -i All_20151104.vcf.gz\n"
        "\n"
        "   # assuming the reference build is correct, just flip to fwd, discarding the rest\n"
        "   bcftools +fixref file.bcf -Ob -o out.bcf -- -d -f ref.fa -m flip\n"
        "\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    memset(&args,0,sizeof(args_t));
    args.skip_rid = -1;
    args.hdr_in  = in;
    args.hdr_out = out;
    args.mode = MODE_STATS;
    args.info_tag = "FIXREF";
    char *ref_fname = NULL;
    static struct option loptions[] =
    {
        {"mode",required_argument,NULL,'m'},
        {"tag-name",required_argument,NULL,'t'},
        {"discard",no_argument,NULL,'d'},
        {"fasta-ref",required_argument,NULL,'f'},
        {"use-id",required_argument,NULL,'i'},
        {NULL,0,NULL,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "?hf:m:di:t:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'm':
                if ( !strcasecmp(optarg,"top") ) args.mode = MODE_TOP2FWD;
                else if ( !strcasecmp(optarg,"flip") ) args.mode = MODE_FLIP2FWD;
                else if ( !strcasecmp(optarg,"flip-all") ) args.mode = MODE_FLIP_ALL;
                else if ( !strcasecmp(optarg,"id") ) args.mode = MODE_USE_ID;
                else if ( !strcasecmp(optarg,"ref-alt") ) args.mode = MODE_REF_ALT;
                else if ( !strcasecmp(optarg,"swap") ) args.mode = MODE_SWAP_REF_ALT;
                else if ( !strcasecmp(optarg,"stats") ) args.mode = MODE_STATS;
                else error("The source strand convention not recognised: %s\n", optarg);
                break;
            case 'i': args.dbsnp_fname = optarg; args.mode = MODE_USE_ID; break;
            case 'd': args.discard = 1; break;
            case 'f': ref_fname = optarg; break;
            case 't': args.info_tag = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( !ref_fname ) error("Expected the -f option\n");
    args.fai = fai_load(ref_fname);
    if ( !args.fai ) error("Failed to load the fai index: %s\n", ref_fname);

    if ( bcf_hdr_printf(args.hdr_out,"##INFO=<ID=%s,Number=.,Type=String,Description=\"The change made by bcftools/fixref\">",args.info_tag) ) return -1;

    if ( args.mode==MODE_STATS ) return 1;
    return 0;
}

static bcf1_t *set_ref_alt(args_t *args, bcf1_t *rec, const char ref, const char alt, int swap_gt)
{
    rec->d.allele[0][0] = ref;
    rec->d.allele[1][0] = alt;
    rec->d.shared_dirty |= BCF1_DIRTY_ALS;

    if ( !swap_gt ) return rec;    // only fix the alleles, leaving GTs unchanged

    int ngts = bcf_get_genotypes(args->hdr_in, rec, &args->gts, &args->ngts);
    if ( ngts<=0 ) return rec;  // no samples, no genotypes

    int i, j, nsmpl = bcf_hdr_nsamples(args->hdr_in);
    ngts /= nsmpl;
    for (i=0; i<nsmpl; i++)
    {
        int32_t *ptr = args->gts + i*ngts;
        for (j=0; j<ngts; j++)
        {
            if ( ptr[j]==bcf_gt_unphased(0) ) ptr[j] = bcf_gt_unphased(1);
            else if ( ptr[j]==bcf_gt_phased(0) ) ptr[j] = bcf_gt_phased(1);
            else if ( ptr[j]==bcf_gt_unphased(1) ) ptr[j] = bcf_gt_unphased(0);
            else if ( ptr[j]==bcf_gt_phased(1) ) ptr[j] = bcf_gt_phased(0);
        }
    }
    bcf_update_genotypes(args->hdr_out,rec,args->gts,args->ngts);

    return rec;
}

static inline int nt2int(char nt)
{
    nt = toupper(nt);
    if ( nt=='A' ) return 0;
    if ( nt=='C' ) return 1;
    if ( nt=='G' ) return 2;
    if ( nt=='T' ) return 3;
    return -1;
}
#define int2nt(x) "ACGT"[x]
#define revint(x) ("3210"[x]-'0')

static int fetch_ref(args_t *args, bcf1_t *rec)
{
    // Get the reference allele
    int len;
    char *ref = faidx_fetch_seq(args->fai, (char*)bcf_seqname(args->hdr_in,rec), rec->pos, rec->pos, &len);
    if ( !ref )
    {
        if ( faidx_has_seq(args->fai, bcf_seqname(args->hdr_in,rec))==0 )
        {
            fprintf(stderr,"Ignoring sequence \"%s\"\n", bcf_seqname(args->hdr_in,rec));
            args->skip_rid = rec->rid;
            return -2;
        }
        error("faidx_fetch_seq failed at %s:%"PRId64"\n", bcf_seqname(args->hdr_in,rec),(int64_t) rec->pos+1);
    }
    int ir = nt2int(*ref);
    free(ref);
    return ir;
}

static void dbsnp_destroy(args_t *args)
{
    if ( !args->i2m ) return;
    khint_t k;
    for (k = 0; k < kh_end(args->i2m); ++k)
        if (kh_exist(args->i2m, k)) free((char*)kh_key(args->i2m, k));
    kh_destroy(i2m, args->i2m);
}
static void dbsnp_init(args_t *args, const char *chr)
{
    dbsnp_destroy(args);
    args->i2m = kh_init(i2m);
    bcf_srs_t *sr = bcf_sr_init();
    if ( bcf_sr_set_regions(sr, chr, 0) != 0 ) goto done;
    if ( !args->dbsnp_fname ) error("No ID file specified, use -i/--use-id\n");
    if ( !bcf_sr_add_reader(sr,args->dbsnp_fname) ) error("Failed to open %s: %s\n", args->dbsnp_fname,bcf_sr_strerror(sr->errnum));
    while ( bcf_sr_next_line(sr) )
    {
        bcf1_t *rec = bcf_sr_get_line(sr, 0);
        if ( rec->d.allele[0][1]!=0 || rec->d.allele[1][1]!=0 ) continue;   // skip non-snps

        int ref = nt2int(rec->d.allele[0][0]);
        if ( ref<0 ) continue;     // non-[ACGT] base

        if ( !rec->d.id || (rec->d.id[0]=='.' && !rec->d.id[1]) ) continue;
        char *id = strdup(rec->d.id);

        int ret, k;
        k = kh_put(i2m, args->i2m, id, &ret);
        if ( ret<0 ) error("An error occurred while inserting the key \"%s\"\n", id);
        if ( ret==0 ) { free(id); continue; }       // skip ambiguous id
        kh_val(args->i2m, k).pos = (uint32_t)rec->pos;
        kh_val(args->i2m, k).ref = ref;
    }
done:
    bcf_sr_destroy(sr);
}

static bcf1_t *dbsnp_check(args_t *args, bcf1_t *rec, int ir, int ia, int ib)
{
    int k, ref,pos;
    char *id = rec->d.id;
    if ( !id || (id[0]=='.' && !id[1]) ) goto no_info;

    k = kh_get(i2m, args->i2m, id);
    if ( k==kh_end(args->i2m) ) goto no_info;

    pos = (int)kh_val(args->i2m, k).pos;
    if ( pos != rec->pos )
    {
        rec->pos = pos;
        ir = fetch_ref(args, rec);
        args->npos_err++;
    }

    ref = kh_val(args->i2m, k).ref;
	if ( ref!=ir )
        error("Reference base mismatch at %s:%"PRId64" .. %c vs %c\n",bcf_seqname(args->hdr_in,rec),(int64_t) rec->pos+1,int2nt(ref),int2nt(ir));

    if ( ia==ref ) { args->dirty = FIX_NONE; return rec; }
    if ( ib==ref ) { args->dirty = FIX_SWAP; args->nswap++; return set_ref_alt(args,rec,int2nt(ib),int2nt(ia),1); }

no_info:
    args->nunresolved++;
    return args->discard ? NULL : rec;
}

static bcf1_t *process_record(bcf1_t *rec)
{
    if ( rec->rid == args.skip_rid ) return NULL;

    bcf1_t *ret = args.mode==MODE_STATS ? NULL : rec;
    args.nsite++;

    args.dirty = FIX_SKIP;

    // Skip non-SNPs
    if ( bcf_get_variant_types(rec)!=VCF_SNP )
    {
        args.nonSNP++;
        return args.discard ? NULL : ret;
    }

    // Get the reference allele
    int ir = fetch_ref(&args, rec);
    if ( ir==-2 ) return NULL;
    if ( ir==-1 )
    {
        args.nonACGT++;
        return args.discard ? NULL : ret;     // not A,C,G,T
    }

    if ( rec->n_allele!=2 )
    {
        // not a biallelic site
        args.nonbiallelic++;
        return args.discard ? NULL : ret;
    }

    int ia = nt2int(rec->d.allele[0][0]);
    if ( ia<0 )
    {
        // not A,C,G,T
        args.nonACGT++;
        return args.discard ? NULL : ret;
    }

    int ib = nt2int(rec->d.allele[1][0]);
    if ( ib<0 )
    {
        // not A,C,G,T
        args.nonACGT++;
        return args.discard ? NULL : ret;
    }

    if ( ia==ib )
    {
        // should not happen in well-formed VCF
        args.nonSNP++;
        return args.discard ? NULL : ret;
    }
    args.count[ia][ib]++;

    if ( ir==ia ) args.nok++;

    if ( args.mode==MODE_USE_ID )
    {
        if ( !args.i2m || args.rid!=rec->rid )
        {
            args.pos = 0;
            args.rid = rec->rid;
            dbsnp_init(&args,bcf_seqname(args.hdr_in,rec));
        }
        ret = dbsnp_check(&args, rec, ir,ia,ib);
        if ( !args.unsorted && args.pos > rec->pos )
        {
            fprintf(stderr,
                "Warning: corrected position(s) results in unsorted VCF, for example %s:%"PRId64" comes after %s:%d\n"
                "         The command `bcftools sort` can be used to fix the order.\n",
                bcf_seqname(args.hdr_in,rec),(int64_t) rec->pos+1,bcf_seqname(args.hdr_in,rec),args.pos);
            args.unsorted = 1;
        }
        args.pos = rec->pos;
        return ret;
    }
    if ( args.mode==MODE_REF_ALT ) // only change the REF/ALT column, leave the genotypes as is
    {
        if ( ir==ia ) { args.dirty = FIX_NONE; return ret; }
        if ( ir==ib ) { args.dirty = FIX_SWAP; args.nswap++; return set_ref_alt(&args,rec,int2nt(ib),int2nt(ia),0); }
        if ( ir==revint(ia) ) { args.dirty = FIX_FLIP; args.nflip++; return set_ref_alt(&args,rec,int2nt(revint(ia)),int2nt(revint(ib)),0); }
        if ( ir==revint(ib) ) { args.dirty = FIX_FLIP|FIX_SWAP; args.nflip_swap++; return set_ref_alt(&args,rec,int2nt(revint(ib)),int2nt(revint(ia)),0); }
        args.dirty = FIX_ERR;
        args.nerr++;
        return ret;
    }
    if ( args.mode==MODE_SWAP_REF_ALT ) // only swap the REF/ALT column but never flip, leave the genotypes as is
    {
        if ( ir==ia ) { args.dirty = FIX_NONE; return ret; }
        if ( ir==ib ) { args.dirty = FIX_SWAP; args.nswap++; return set_ref_alt(&args,rec,int2nt(ib),int2nt(ia),0); }
        args.dirty = FIX_ERR;
        args.nerr++;
        return ret;
    }
    if ( args.mode==MODE_FLIP2FWD || args.mode==MODE_FLIP_ALL )
    {
        int pair = 1 << ia | 1 << ib;
        if ( args.mode==MODE_FLIP2FWD && (pair==0x9 || pair==0x6) )   // skip ambiguous pairs: A/T or C/G
        {
            args.nunresolved++;
            return args.discard ? NULL : ret;
        }
        if ( ir==ia ) { args.dirty = FIX_NONE; return ret; }
        if ( ir==ib ) { args.dirty = FIX_SWAP|FIX_GT; args.nswap++; return set_ref_alt(&args,rec,int2nt(ib),int2nt(ia),1); }
        if ( ir==revint(ia) ) { args.dirty = FIX_FLIP; args.nflip++; return set_ref_alt(&args,rec,int2nt(revint(ia)),int2nt(revint(ib)),0); }
        if ( ir==revint(ib) ) { args.dirty = FIX_FLIP|FIX_SWAP|FIX_GT; args.nflip_swap++; return set_ref_alt(&args,rec,int2nt(revint(ib)),int2nt(revint(ia)),1); }
        args.dirty = FIX_ERR;
        args.nerr++;
        return ret;
    }
    if ( args.mode==MODE_TOP2FWD )
    {
        int pair = 1 << ia | 1 << ib;
        if ( pair != 0x9 && pair != 0x6 )    // unambiguous pair: A/C or A/G
        {
            if ( ir==ia ) return ret;

            int ia_rev = revint(ia);
            if ( ir==ia_rev )               // vcfref is A, faref is T, flip
            {
                args.dirty = FIX_FLIP;
                args.nflip++;
                return set_ref_alt(&args,rec,int2nt(ia_rev),int2nt(revint(ib)),0);
            }
            if ( ir==ib )                   // vcfalt is faref, swap
            {
                args.dirty = FIX_SWAP|FIX_GT;
                args.nswap++;
                return set_ref_alt(&args,rec,int2nt(ib),int2nt(ia),1);
            }
            if ( ib!=revint(ir) ) { args.dirty = FIX_ERR; args.nerr++; return ret; }

            args.dirty = FIX_FLIP|FIX_SWAP|FIX_GT;
            args.nflip_swap++;
            return set_ref_alt(&args,rec,int2nt(revint(ib)),int2nt(revint(ia)),1);
        }
        else    // ambiguous pair, sequence walking must be performed
        {
            int len, win = rec->pos > 100 ? 100 : rec->pos, beg = rec->pos - win, end = rec->pos + win;
            char *ref = faidx_fetch_seq(args.fai, (char*)bcf_seqname(args.hdr_in,rec), beg,end, &len);
            if ( !ref ) error("faidx_fetch_seq failed at %s:%"PRId64"\n", bcf_seqname(args.hdr_in,rec),(int64_t) rec->pos+1);
            if ( end - beg + 1 != len ) error("FIXME: check win=%d,len=%d at %s:%"PRId64"  (%d %d)\n", win,len, bcf_seqname(args.hdr_in,rec),(int64_t) rec->pos+1, end,beg);

            int i, mid = rec->pos - beg, strand = 0;
            for (i=1; i<=win; i++)
            {
                int ra = nt2int(ref[mid-i]);
                int rb = nt2int(ref[mid+i]);
                if ( ra<0 || rb<0 || ra==rb ) continue;     // skip N's and non-infomative pairs: A/A, C/C, G/G, T/T
                pair = 1 << ra | 1 << rb;
                if ( pair==0x9 || pair==0x6 ) continue;     // skip ambiguous pairs: A/T or C/G
                strand = 1 << ra & 0x9 ? 1 : -1;
                break;
            }
            free(ref);

            if ( strand==1 )
            {
                if ( ir==ia ) { args.dirty = FIX_NONE; return ret; }
                if ( ir==ib )
                {
                    args.dirty = FIX_SWAP|FIX_GT;
                    args.nswap++;
                    return set_ref_alt(&args,rec,int2nt(ib),int2nt(ia),1);
                }
            }
            else if ( strand==-1 )
            {
                int ia_rev = revint(ia);
                int ib_rev = revint(ib);
                if ( ir==ia_rev )
                {
                    args.dirty = FIX_FLIP;
                    args.nflip++;
                    return set_ref_alt(&args,rec,int2nt(ia_rev),int2nt(ib_rev),0);
                }
                if ( ir==ib_rev )
                {
                    args.dirty = FIX_FLIP|FIX_SWAP|FIX_GT;
                    args.nflip_swap++;
                    return set_ref_alt(&args,rec,int2nt(ib_rev),int2nt(ia_rev),1);
                }
            }

            args.nunresolved++;
            return args.discard ? NULL : ret;
        }
    }
    return ret;
}

bcf1_t *process(bcf1_t *rec)
{
    args.dirty = 0;
    bcf1_t *ret = process_record(rec);
    if ( ret && args.dirty )
    {
        args.str.l = 0;
        int i;
        for (i=0; i<6; i++)
        {
            if ( !(args.dirty & (1<<i)) ) continue;
            if ( args.str.l ) kputc(',', &args.str);
            kputs(info_annots[i], &args.str);
        }
        bcf_update_info_string(args.hdr_out,rec,args.info_tag,args.str.s);
    }
    return ret;
}


int top_mask[4][4] =
{
    {0,1,1,1},
    {0,0,1,0},
    {0,0,0,0},
    {0,0,0,0},
};
int bot_mask[4][4] =
{
    {0,0,0,0},
    {0,0,0,0},
    {0,1,0,0},
    {1,1,1,0},
};

void destroy(void)
{
    uint32_t i,j,tot = 0;
    uint32_t top_err = 0, bot_err = 0;
    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            tot += args.count[i][j];
            if ( !top_mask[i][j] && args.count[i][j] ) top_err++;
            if ( !bot_mask[i][j] && args.count[i][j] ) bot_err++;
        }
    }
    uint32_t nskip = args.nonACGT+args.nonSNP+args.nonbiallelic;
    uint32_t ncmp  = args.nsite - nskip;

    fprintf(stderr,"# SC, guessed strand convention\n");
    fprintf(stderr,"SC\tTOP-compatible\t%d\n",top_err?0:1);
    fprintf(stderr,"SC\tBOT-compatible\t%d\n",bot_err?0:1);

    fprintf(stderr,"# ST, substitution types\n");
    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            if ( i==j ) continue;
            fprintf(stderr,"ST\t%c>%c\t%u\t%.1f%%\n", int2nt(i),int2nt(j),args.count[i][j], args.count[i][j]*100./tot);
        }
    }
    fprintf(stderr,"# NS, Number of sites:\n");
    fprintf(stderr,"NS\ttotal        \t%u\n", args.nsite);
    fprintf(stderr,"NS\tref match    \t%u\t%.1f%%\n", args.nok,100.*args.nok/ncmp);
    fprintf(stderr,"NS\tref mismatch \t%u\t%.1f%%\n", ncmp-args.nok,100.*(ncmp-args.nok)/ncmp);
    if ( args.mode!=MODE_STATS )
    {
        fprintf(stderr,"NS\tflipped      \t%u\t%.1f%%\n", args.nflip,100.*args.nflip/(args.nsite-nskip));
        fprintf(stderr,"NS\tswapped      \t%u\t%.1f%%\n", args.nswap,100.*args.nswap/(args.nsite-nskip));
        fprintf(stderr,"NS\tflip+swap    \t%u\t%.1f%%\n", args.nflip_swap,100.*args.nflip_swap/(args.nsite-nskip));
        fprintf(stderr,"NS\tunresolved   \t%u\t%.1f%%\n", args.nunresolved,100.*args.nunresolved/(args.nsite-nskip));
        fprintf(stderr,"NS\tfixed pos    \t%u\t%.1f%%\n", args.npos_err,100.*args.npos_err/(args.nsite-nskip));
    }
    fprintf(stderr,"NS\terrors       \t%u\n", args.nerr);
    fprintf(stderr,"NS\tskipped      \t%u\n", nskip);
    fprintf(stderr,"NS\tnon-ACGT     \t%u\n", args.nonACGT);
    fprintf(stderr,"NS\tnon-SNP      \t%u\n", args.nonSNP);
    fprintf(stderr,"NS\tnon-biallelic\t%u\n", args.nonbiallelic);

    free(args.str.s);
    free(args.gts);
    if ( args.fai ) fai_destroy(args.fai);
    dbsnp_destroy(&args);
}
