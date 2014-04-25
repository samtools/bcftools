#include <htslib/hts.h>
#include "version.h"

void version(const char **bcftools_version, const char **htslib_version)
{
    *bcftools_version = BCFTOOLS_VERSION;
    *htslib_version = hts_version();
}

