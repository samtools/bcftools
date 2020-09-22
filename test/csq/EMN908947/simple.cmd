{bin}/bcftools csq -f EMN908947.fa -g simple.gff simple.vcf | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%EXP\n%POS\t%REF\t%ALT\t%BCSQ\n\n'
