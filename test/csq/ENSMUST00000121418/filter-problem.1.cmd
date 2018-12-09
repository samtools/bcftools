{bin}/bcftools csq -f ENSMUST00000121418.fa -g ENSMUST00000121418.gff filter-problem.vcf -i 'POS!=25105' --phase m | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
