#include <stdio.h>
#include <unistd.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include "bcftools.h"

static void usage(void)
{
    fprintf(stderr, "About:     Convert between VCF/BCF.\n");   // Most of the original functionality was moved into vcfsubset
    fprintf(stderr, "Usage:     bcftools view [options] <in.bcf>|<in.vcf>|<in.vcf.gz>|-\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -l INT       compression level [%d]\n", -1);
    fprintf(stderr, "   -n FILE      output file name [stdout]\n");
    fprintf(stderr, "   -o TYPE      output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n");
    fprintf(stderr, "   -v           stdin input is VCF\n");
    fprintf(stderr, "\n");
}

int main_vcfview(int argc, char *argv[])
{
	int c, clevel = -1, in_type = FT_BCF, out_type = FT_VCF;
	char *fname_out = NULL, moder[8], modew[8];

	while ((c = getopt(argc, argv, "l:bvo:n:z?hu")) >= 0) {
		switch (c) {
            case 'o': 
                switch (optarg[0]) {
                    case 'b': out_type = FT_BCF_GZ; break;
                    case 'u': out_type = FT_BCF; break;
                    case 'z': out_type = FT_VCF_GZ; break;
                    case 'v': out_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                }
                break;
            case 'l': clevel = atoi(optarg); out_type |= FT_GZ; break;
            case 'v': in_type  = FT_VCF; break;
            case 'b': out_type = FT_BCF_GZ; break;
            case 'u': out_type = FT_BCF; break;
            case 'z': out_type = FT_VCF_GZ; break;
            case 'n': fname_out = optarg; break;
            case '?':
            case 'h': usage(); return 1; break;
        }
    }
	if (argc!=optind+1) { usage(); return 1; }

    // Init reader
	strcpy(moder, "r");
	if ( (!strcmp("-",argv[optind]) && (in_type & FT_BCF)) || (hts_file_type(argv[optind]) & FT_BCF)) strcat(moder, "b");
	htsFile *fp_in = hts_open(argv[optind], moder, NULL);
	bcf_hdr_t *hdr = vcf_hdr_read(fp_in);
    if ( !hdr ) error("Fail to read VCF/BCF header: %s\n", argv[optind]); 
	bcf1_t *rec = bcf_init1();

    // Init writer
    strcpy(modew, "w");
    if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
    if (out_type & FT_GZ) strcat(modew,"z");
    if (out_type & FT_BCF) strcat(modew, "b");
    if (out_type == FT_BCF) strcat(modew, "u"); // uncompressed BCF output
    htsFile *fp_out = hts_open(fname_out ? fname_out : "-", modew, NULL);

    vcf_hdr_write(fp_out, hdr);
    while ( vcf_read1(fp_in, hdr, rec) >= 0) vcf_write1(fp_out, hdr, rec);

	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);
	hts_close(fp_in);
    hts_close(fp_out);

	return 0;
}

