#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

bcf_hdr_t *hdr;
int *arr = NULL, marr = 0;

const char *about(void)
{
    return "Fill INFO fields AN and AC.\n";
}

int init(bcf_hdr_t *h)
{
    hdr = h;
    bcf_hdr_append(hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
    bcf_hdr_append(hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    return 0;
}

int process(bcf1_t *rec)
{
    hts_expand(int,rec->n_allele,marr,arr);
    int ret = bcf_calc_ac(hdr,rec,arr,BCF_UN_FMT);
    if ( ret>0 )
    {
        int i, an = 0;
        for (i=0; i<rec->n_allele; i++) an += arr[i];
        bcf_update_info_int32(hdr, rec, "AN", &an, 1);
        bcf_update_info_int32(hdr, rec, "AC", arr+1, rec->n_allele-1);
    }
    return 0;
}

void destroy(void) 
{
    free(arr);
}


