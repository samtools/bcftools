{bin}/bcftools csq -f ENST00000368801.2.fa -g ENST00000368801.2.gff alt-gencode.vcf -v 0 -C 1 | {bin}/bcftools query -f'%BCSQ\n' | ../sort-csq -q
