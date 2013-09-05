#ifndef __FILTER_H__
#define __FILTER_H__

#include <htslib/vcf.h>

typedef struct _filter_t filter_t;

/**
  *  @hdr:  BCF header file
  *  @str:  see the bcftools filter command help for description
  */
filter_t *filter_init(bcf_hdr_t *hdr, const char *str);

void filter_destroy(filter_t *filter);

/**
  *  filter_test() - test whether the BCF record passes the test
  *  Returns 1 if the expression is true and 0 if false. 
  */
int filter_test(filter_t *filter, bcf1_t *rec);

#endif
