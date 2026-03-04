#!/usr/bin/env bash

cd "$(dirname "${BASH_SOURCE[0]}")"

for x in *.vcf; do
    bcftools +trio-dnm3 -p 1X:child,dad,mom $x | bcftools query -uf$x'\t%CHROM:%POS\t%EXP\t[ %DNM]\t[ %VAF]\t[ %AD]\t%REF,%ALT\t[ %VA]\t[ %SP]';
done | sort -k3,3g -k4,4g -k3,3d | sed 's,fp,-, ; s,tp,o,'

