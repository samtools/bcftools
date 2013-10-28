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

typedef struct
{
    aln_aux_t aln;
    char *tseq, *seq;
    int mseq;
    bcf1_t **lines;
    rbuf_t rbuf;
    int buf_win;            // maximum distance between two records to consider
    int aln_win;            // the realignment window size (maximum repeat size)
    bcf_srs_t *files;       // using the synced reader only for -r option
    bcf_hdr_t *hdr;
    faidx_t *fai;
	char **argv, *ref_fname, *vcf_fname, *region;
	int argc, rmdup, output_type;
}
args_t;

void _vcfnorm_debug_print(aln_aux_t *aux)
{
    cell_t *mat = aux->mat;
    char *ref = aux->ref;
    char *seq = aux->seq; 
    int nref  = aux->nref;
    int nseq  = aux->nseq;
    int nlen  = nref>nseq ? nref : nseq;
    int k     = (nref+1)*(nseq+1)-1;
    int kd    = nref+2;
    int i = k/(nref+1);
    int j = k - i*(nref+1);
    assert(i>0 && j>0);
    int l = k, ialn = 0, nout_ref = 0, nout_seq = 0, ipos = 0;
    char *aln_ref = (char*) malloc(sizeof(char)*(nlen+1));
    char *aln_seq = (char*) malloc(sizeof(char)*(nlen+1));
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
        else     // m
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
    while (k>kmin && mat[k].dir==DM) { k -= kd; ipos++; }

    i = k/(nref+1);             // seq[nseq-i]
    j = k - i*(nref+1);         // ref[nref-j]

    if ( !i && !j ) return -1;  // this is a legitimate case, consider MNPs
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

    int i, ref_len = strlen(line->d.allele[0]), len = ref_len;
    for (i=1; i<line->n_allele; i++)
    {
        int l = strlen(line->d.allele[i]);
        if ( len < l ) len = l;
    }
    if ( len==1 ) return 0;    // SNP

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
    int ref_winlen;
    char *ref = faidx_fetch_seq(args->fai, (char*)args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos-win, line->pos+ref_len, &ref_winlen);
    if ( !ref ) error("faidx_fetch_seq failed at %s:%d\n", args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos-win);
    assert( ref_winlen==ref_len+win+1 );

    // Sanity check: the reference sequence must match the REF allele
    if ( strncasecmp(&ref[win],line->d.allele[0],ref_len) )
    {
        for (i=0; i<ref_len; i++)
            if ( toupper(ref[win+i])!=toupper(line->d.allele[0][i]) ) break;
        error("\nSanity check failed, the reference sequence differs at %s:%d[%d] .. '%c' vs '%c'\n", 
            args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1, i+1,toupper(ref[win+i]),toupper(line->d.allele[0][i]));
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
        args->aln.nref = ref_winlen;
        args->aln.nseq = i;

        if ( align(args, &args->aln)<0 ) 
        {
            // something went wrong - output the original line
            free(ref);
            return 0;
        }

        // fprintf(stderr, "%s  \t nref=%d\n", ref, ref_winlen);
        // fprintf(stderr, "%s  \t nseq=%d\n", args->tseq, i);
        // fprintf(stderr, "pos=%d win=%d  ipos=%d lref=%d lseq=%d\n", line->pos+1, win, args->aln.ipos, args->aln.lref, args->aln.lseq);
        // fprintf(stderr, "-> "); for (k=args->aln.ipos; k<args->aln.lref; k++) fprintf(stderr, "%c", ref[k]); fprintf(stderr, "\n");
        // fprintf(stderr, "-> "); for (k=args->aln.ipos; k<args->aln.lseq; k++) fprintf(stderr, "%c", args->tseq[k]); fprintf(stderr, "\n");
        // fprintf(stderr, "\n"); 

        ipos[j] = args->aln.ipos;   // position before the first difference (w.r.t. the window)
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

    // todo: 
    //      - modify the alleles only if needed. For now redoing always to catch errors
    //      - in some cases the realignment does not really improve things, see the case at 2:114 in test/norm.vcf

    // REF
    kstring_t str = {0,0,0};
    kputsn_(&ref[min_pos], max_ref-min_pos, &str); 
    kputc_(0, &str);
    // ALTs
    for (k=1; k<line->n_allele; k++)
    {
        // prefix the sequence with REF bases if the other alleles were aligned more to the left
        int nprefix = ipos[k] - min_pos;
        if ( nprefix ) kputsn_(&ref[min_pos], nprefix, &str);

        // the ALT sequence 
        int nseq = iseq[k] - ipos[k];
        if ( nseq )
        {
            char *alt = args->seq + (k-1)*len + ipos[k];
            kputsn_(alt, nseq, &str);
        }

        // suffix invoked by other deletions which must be added to match the REF
        int nsuffix = max_ref - iref[k];
        if ( nsuffix ) kputsn_(&ref[iref[k]], nsuffix, &str);
        kputc_(0, &str);
    }
    // create new block of alleles 
    char *rmme = line->d.als;
    line->d.allele[0] = line->d.als = str.s;
    line->d.m_als = str.m;
    char *t = str.s;
    for (k=1; k<line->n_allele; k++)
    {
        while (*t) t++;
        line->d.allele[k] = ++t;
    }
    free(rmme);
    free(ref);
    return 1;
}

void flush_buffer(args_t *args, htsFile *file, int n)
{
    int i, k, prev_rid = -1, prev_pos = 0, prev_type = 0;
    for (i=0; i<n; i++)
    {
        k = rbuf_shift(&args->rbuf);
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
}

static void init_data(args_t *args)
{
    args->hdr = args->files->readers[0].header;
    rbuf_init(&args->rbuf, 100);
    args->lines = (bcf1_t**) calloc(args->rbuf.m, sizeof(bcf1_t*));
    args->fai = fai_load(args->ref_fname);
    if ( !args->fai ) error("Failed to load the fai index: %s\n", args->ref_fname);
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->rbuf.m; i++)
        if ( args->lines[i] ) bcf_destroy1(args->lines[i]);
    free(args->lines);
    fai_destroy(args->fai);
    if ( args->mseq ) free(args->seq);
    if ( args->aln.nmat ) free(args->aln.mat);
    if ( args->aln.m_arr ) { free(args->aln.ipos_arr); free(args->aln.lref_arr); free(args->aln.lseq_arr); }
}


#define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }
static void normalize_vcf(args_t *args)
{
    htsFile *out = hts_open("-", hts_bcf_wmode(args->output_type));
    bcf_hdr_append_version(args->hdr, args->argc, args->argv, "bcftools_norm");
    bcf_hdr_write(out, args->hdr);

    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = args->files->readers[0].buffer[0];
        if ( realign(args, line)<0 ) continue;   // exclude broken VCF lines

        // still on the same chromosome?
        int i, j, ilast = rbuf_last(&args->rbuf); 
        if ( ilast>=0 && line->rid != args->lines[ilast]->rid ) flush_buffer(args, out, args->rbuf.n); // new chromosome

        // insert into sorted buffer
        i = j = ilast = rbuf_add(&args->rbuf);
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
}

static void usage(void)
{
	fprintf(stderr, "About:   Left-align and normalize indels.\n");
	fprintf(stderr, "Usage:   bcftools norm [options] -f ref.fa <file.vcf.gz>\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "    -D, --remove-duplicates           remove duplicate lines of the same type. [Todo: merge genotypes, don't just throw away.]\n");
	fprintf(stderr, "    -f, --fasta-ref <file>            reference sequence\n");
    fprintf(stderr, "    -o, --output-type <type>          'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n");
	fprintf(stderr, "    -r, --region <file|reg>           comma-separated list of regions or regions listed in tab-delimited indexed file\n");
	fprintf(stderr, "    -w, --win <int,int>               alignment window and buffer window [50,1000]\n");
	fprintf(stderr, "\n");
	exit(1);
}

int main_vcfnorm(int argc, char *argv[])
{
	int c;
	args_t *args  = (args_t*) calloc(1,sizeof(args_t));
	args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->aln_win = 50;
    args->buf_win = 1000;

	static struct option loptions[] = 
	{
		{"help",0,0,'h'},
		{"fasta-ref",1,0,'f'},
		{"region",1,0,'r'},
		{"win",1,0,'w'},
		{"remove-duplicates",0,0,'D'},
        {"output-type",1,0,'o'},
		{0,0,0,0}
	};
	while ((c = getopt_long(argc, argv, "hr:f:w:Do:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'o': 
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                }
                break;
			case 'D': args->rmdup = 1; break;
			case 'f': args->ref_fname = optarg; break;
			case 'r': args->region = optarg; break;
            case 'w': { if (sscanf(optarg,"%d,%d",&args->aln_win,&args->buf_win)!=2) error("Could not parse --win %s\n", optarg); break; }
			case 'h': 
			case '?': usage();
			default: error("Unknown argument: %s\n", optarg);
		}
	}
	if ( argc!=optind+1 || !args->ref_fname ) usage();   // none or too many files given
    if ( args->region )
    {
        if ( bcf_sr_set_targets(args->files, args->region,0)<0 ) error("Failed to read the targets: %s\n", args->region);
    }

    if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open or the file not indexed: %s\n", argv[optind]);
    init_data(args);
    normalize_vcf(args);
    destroy_data(args);
    bcf_sr_destroy(args->files);
	free(args);
	return 0;
}

