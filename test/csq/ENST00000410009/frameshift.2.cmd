{bin}/bcftools csq -lb -f ENST00000410009.fa -g ENST00000410009.gff frameshift.vcf | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
