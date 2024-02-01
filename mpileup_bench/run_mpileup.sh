#!/bin/sh

bcftools=${BCFTOOLS:-bcftools}

if [ $# -lt 3 ]
then
    echo Usage: run_mpileup file.bam region outdir [bcftools-arguments]
    exit 1
fi

file=$1
region=$2
dir=$3
shift 3
args=${@+"$@"}
echo $file / $region / $dir / $args

mkdir -p $dir

# Standard files, downloaded by the get_data.sh script.
# They can be overridden on the command line before running the script
# eg "BED=hard_regions.bed ./run_mpileup.sh [opts]"
HREF38=${HREF38:-/nfs/srpipe_references/references/Human/GRCh38_full_analysis_set_plus_decoy_hla/all/fasta/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa}
TRUTH=${TRUTH:-HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}
BED=${BED:-HG005_GRCh38_1_22_v4.2.1_benchmark.bed}


# -- Do the actual variant calling
echo "=== Running $bcftools mpileup $args -a AD --fasta-ref $HREF38 -r $region $file | bcftools call -vm -"
eval $bcftools mpileup $args -a AD --fasta-ref $HREF38 -r $region $file 2>$dir/bcftools.mpileup.out | $bcftools call -vm - > $dir/bcftools.vcf 2>$dir/bcftools.call.out

# A primary evaluation.
# They key thing here is it leaves behind the .isec directory with the
# intersection of the truth and call sets.
echo "=== ./compare_vcf_simple.sh $TRUTH $dir/bcftools.vcf "" $BED $region"
QUAL=30 NORM=1 ./compare_vcf_simple.sh $TRUTH $dir/bcftools.vcf "" $BED $region

# Produce a .plot file for use in gnuplot, along with a basic summary too.
echo "=== ./plot_isec.pl $dir/bcftools.vcf.isec indel > $dir/plot++"
./plot_isec.pl $dir/bcftools.vcf.isec indel > $dir/plot++
grep ALL $dir/plot++ > $dir/plot
grep INS $dir/plot++ > $dir/plot.ins
grep DEL $dir/plot++ > $dir/plot.del

awk 'BEGIN {n=0} $6 >= n {print;n=50*(1+int($6/50))}' $dir/plot | cut -c 1-28|head -20


# Example gnuplot:
# set xlabel "FPr"
# set ylabel "FNr"
# set yrange [0:10]
# set title "PacBio 50x indels"

# FP vs FN
# a=8;b=10;t="ALL";plot "<grep ".t." dev/plot" using a:b with linespoints title "dev" lw 2 ps .6, "<grep ".t." cns/plot" using a:b with linespoints title "indels-cns" lw 2 ps .6

# GT vs FN
# a=9;b=10; ... as above
