{bin}/bcftools csq -f ENST00000000001.fa -g ENST00000000001.gff test.vcf -p a | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
