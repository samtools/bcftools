{bin}/bcftools csq -f EMN908947.fa -g ribosomal-slippage.gff ribosomal-slippage.vcf --force | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%EXP\n%POS\t%REF\t%ALT\t%BCSQ\n\n'
