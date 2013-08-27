#ifndef __BCFTOOLS_H__
#define __BCFTOOLS_H__

#include <stdarg.h>
#include <htslib/vcf.h>
#include "version.h"

char *bcftools_version(void);
void error(const char *format, ...);
void bcf_hdr_append_version(bcf_hdr_t *hdr, int argc, char **argv, const char *cmd);

#endif
