#!/usr/bin/env bash

cd "$(dirname "${BASH_SOURCE[0]}")"

BCFTOOLS=${1:-bcftools}
export LC_ALL=C

for x in *.vcf; do
    "$BCFTOOLS" +trio-dnm3 -p 1X:child,dad,mom $x | "$BCFTOOLS" query -uf$x'\t%CHROM:%POS\t%EXP\t[ %DNM]\t[ %VAF]\t[ %AD]\t%REF,%ALT\t[ %VA]\t[ %SP]';
done | sort | sed 's/fp/-/;s/tp/o/'

