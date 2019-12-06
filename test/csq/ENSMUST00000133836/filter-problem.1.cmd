{bin}/bcftools csq -f ENSMUST00000133836.fa -g ENSMUST00000133836.gff filter-problem.vcf -i 'FILTER="PASS"' --phase m | ../sort-csq | {bin}/bcftools query -f'%POS\t%REF\t%ALT\t%BCSQ\n'
