{bin}/bcftools csq -f ENST00000422159.fa -g ENST00000422159.gff multiallelic.vcf  | {bin}/bcftools query -f'[%POS\t%REF\t%ALT\t%GT\t%TBCSQ\n]' | ../sort-csq -q
echo
{bin}/bcftools csq -f ENST00000422159.fa -g short.gff multiallelic.vcf | {bin}/bcftools query -f'[%POS\t%REF\t%ALT\t%GT\t%TBCSQ\n]' | ../sort-csq -q
