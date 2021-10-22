{bin}/bcftools csq -p a -f ENST00000000000.fa -g ENST00000000000.gff upstream-stop.vcf  | {bin}/bcftools query -f'[%POS\t%REF\t%ALT\t%SAMPLE\t%GT\t%TBCSQ\n]' | ../sort-csq -q
