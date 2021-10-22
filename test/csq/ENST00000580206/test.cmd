{bin}/bcftools csq -f ENST00000580206.fa -g ENST00000580206.gff compound.ins-snv.vcf  | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
echo
{bin}/bcftools csq -f ENST00000580206.fa -g short.gff compound.ins-snv.vcf  | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
echo
{bin}/bcftools csq -f ENST00000580206.fa -g ENST00000580206.gff single.ins-snv.vcf  | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
echo
{bin}/bcftools csq -f ENST00000580206.fa -g short.gff single.ins-snv.vcf  | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
echo
{bin}/bcftools csq -f ENST00000580206.fa -g ENST00000580206.gff compound.del-snv.vcf  | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
echo
{bin}/bcftools csq -f ENST00000580206.fa -g short.gff compound.del-snv.vcf  | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
echo
{bin}/bcftools csq -f ENST00000580206.fa -g short.gff single.del-snv.vcf  | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
echo
{bin}/bcftools csq -f ENST00000580206.fa -g ENST00000580206.gff compound.del-ins.vcf  | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
echo
{bin}/bcftools csq -f ENST00000580206.fa -g short.gff compound.del-ins.vcf  | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
echo
{bin}/bcftools csq -f ENST00000580206.fa -g short.gff single.del-ins.vcf  | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
