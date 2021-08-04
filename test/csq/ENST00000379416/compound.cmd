{bin}/bcftools csq -B 3 -f ENST00000379416.fa -g ENST00000379416.gff compound.del-snp.vcf       | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
echo
{bin}/bcftools csq -B 3 -f ENST00000379416.fa -g ENST00000379416.gff compound.ins-del.vcf       | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
echo
{bin}/bcftools csq -B 3 -f ENST00000379416.fa -g ENST00000379416.gff compound.ins-ins.vcf       | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
echo
{bin}/bcftools csq -B 3 -f ENST00000379416.fa -g ENST00000379416.gff compound.ins-snp.vcf       | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
echo
{bin}/bcftools csq -B 3 -f ENST00000379416.fa -g ENST00000379416.gff single.vcf                 | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
