#!/bin/sh

# Execute bcftools on HG002, HG003 and HG004 samples.  This could be done as
# a trio with pedegree information, but for now we just do them without
# additional data so we can evaluate calling in the presence of other samples.
#
# We have the GIAB benchmark files, but the HG002-4 BAMs need to be passed in
# as inputs as this may vary by instrument.

# Standard files, downloaded by the get_data.sh script.
# They can be overridden on the command line before running the script
# eg "BED=hard_regions.bed ./run_mpileup.sh [opts]"

# Note if in the shell you do:  === () { /bin/true; }
# then you can cut and paste the === lines as the first part becomes a
# semi-colon capable comment.

HREF38=${HREF38:-/nfs/srpipe_references/references/Human/GRCh38_full_analysis_set_plus_decoy_hla/all/fasta/Homo_sapiens.GRCh38_full_analysis_set_plus_decoy_hla.fa}
TRUTH002=${TRUTH002:-HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}
TRUTH003=${TRUTH003:-HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}
TRUTH004=${TRUTH004:-HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}
BED002=${BED002:-HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed}
BED003=${BED003:-HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed}
BED004=${BED004:-HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed}
bcftools=${BCFTOOLS:-bcftools}
BENCHDIR=${BENCHDIR:-.}

if [ $# -lt 5 ]
then
    echo Usage: run_multi.sh HG002.bam HG003.bam HG004.bam region outdir [bcftools-arguments]
    echo
    echo Internal variables may be overridden on the command line, eg
    echo BCFTOOLS=~/bcftools/bcftools.myPR ./run_multi.sh [args]
    echo
    echo Current settings:
    echo BCFTOOLS=${bcftools}
    echo HREF38=${HREF38}
    echo TRUTH002=${TRUTH002}
    echo TRUTH003=${TRUTH003}
    echo TRUTH004=${TRUTH004}
    echo BED002=${BED002}
    echo BED003=${BED003}
    echo BED004=${BED004}
    echo BENCHDIR=${BENCHDIR}
    exit 1
fi


file2=$1
file3=$2
file4=$3
region=$4
dir=$5
shift 5
args=${@+"$@"}

mkdir -p $dir

# Change to "if false" if we wish to adjust the plots only without recalling.
if true
then

# Do the actual variant calling
echo "===; Running $bcftools mpileup $args -a AD --fasta-ref $HREF38 -r $region $file2 $file3 $file4 | bcftools call -vm -"
eval $bcftools mpileup $args -a AD --fasta-ref $HREF38 -r $region $file2 $file3 $file4 2>$dir/bcftools.mpileup.out | $bcftools call -vm - > $dir/bcftools.vcf 2>$dir/bcftools.call.out

# Split into the input samples
echo
echo "===; $bcftools +split $dir/bcftools.vcf -o $dir"
eval $bcftools +split $dir/bcftools.vcf -o $dir

set -- `$bcftools view -h $dir/bcftools.vcf|tail -1|cut -f 10-`
vcf2=$1.vcf
vcf3=$2.vcf
vcf4=$3.vcf

# Remove non-sample data
echo
echo "===; $bcftools norm -m -both -f $HREF38 $dir/$vcf2 | grep -v '[.0]/[.0]:' > $dir/_$vcf2"
eval $bcftools norm -m -both -f $HREF38 $dir/$vcf2 | grep -v '[.0]/[.0]:' > $dir/_$vcf2
echo "===; $bcftools norm -m -both -f $HREF38 $dir/$vcf3 | grep -v '[.0]/[.0]:' > $dir/_$vcf3"
eval $bcftools norm -m -both -f $HREF38 $dir/$vcf3 | grep -v '[.0]/[.0]:' > $dir/_$vcf3
echo "===; $bcftools norm -m -both -f $HREF38 $dir/$vcf4 | grep -v '[.0]/[.0]:' > $dir/_$vcf4"
eval $bcftools norm -m -both -f $HREF38 $dir/$vcf4 | grep -v '[.0]/[.0]:' > $dir/_$vcf4

# A primary evaluation.
# They key thing here is it leaves behind the .isec directory with the
# intersection of the truth and call sets.
echo
echo "===; QUAL=30 NORM=1 ${BENCHDIR}/compare_vcf_simple.sh $TRUTH002 $dir/_$vcf2 '' $BED002 $region"
QUAL=30 NORM=1 ${BENCHDIR}/compare_vcf_simple.sh $TRUTH002 $dir/_$vcf2 "" $BED002 $region
echo "===; QUAL=30 NORM=1 ${BENCHDIR}/compare_vcf_simple.sh $TRUTH003 $dir/_$vcf3 '' $BED003 $region"
QUAL=30 NORM=1 ${BENCHDIR}/compare_vcf_simple.sh $TRUTH003 $dir/_$vcf3 "" $BED003 $region
echo "===; QUAL=30 NORM=1 ${BENCHDIR}/compare_vcf_simple.sh $TRUTH004 $dir/_$vcf4 '' $BED004 $region"
QUAL=30 NORM=1 ${BENCHDIR}/compare_vcf_simple.sh $TRUTH004 $dir/_$vcf4 "" $BED004 $region


# Produce a .plot file for use in gnuplot, along with a basic summary too.
echo "===; ${BENCHDIR}/plot_isec.pl $dir/_$vcf2.isec indel > $dir/2.all"
${BENCHDIR}/plot_isec.pl $dir/_$vcf2.isec indel > $dir/2.out
grep ALL $dir/2.out > $dir/2.all
grep INS $dir/2.out > $dir/2.ins
grep DEL $dir/2.out > $dir/2.del
awk 'BEGIN {n=0} $6 >= n {print;n=50*(1+int($6/50))}' $dir/2.all | cut -c 1-28|head -20

echo "===; ${BENCHDIR}/plot_isec.pl $dir/_$vcf3.isec indel > $dir/3.all"
${BENCHDIR}/plot_isec.pl $dir/_$vcf3.isec indel > $dir/3.out
grep ALL $dir/3.out > $dir/3.all
grep INS $dir/3.out > $dir/3.ins
grep DEL $dir/3.out > $dir/3.del
awk 'BEGIN {n=0} $6 >= n {print;n=50*(1+int($6/50))}' $dir/3.all | cut -c 1-28|head -20

echo "===; ${BENCHDIR}/plot_isec.pl $dir/_$vcf4.isec indel > $dir/4.all"
${BENCHDIR}/plot_isec.pl $dir/_$vcf4.isec indel > $dir/4.out
grep ALL $dir/4.out > $dir/4.all
grep INS $dir/4.out > $dir/4.ins
grep DEL $dir/4.out > $dir/4.del
awk 'BEGIN {n=0} $6 >= n {print;n=50*(1+int($6/50))}' $dir/4.all | cut -c 1-28|head -20

else

# Short cut to replot results without rerunning the calling
set -- `$bcftools view -h $dir/bcftools.vcf|tail -1|cut -f 10-`
vcf2=$1.vcf
vcf3=$2.vcf
vcf4=$3.vcf

fi

# Generate GNU plots
echo === Running gnuplot to create "HG00[234]*.png"
gnuplot <<EOF
set xlabel "FPr"
set ylabel "FNr"
set yrange [0:10]
set title "indels"

set terminal png size 800,600

set output "$dir/HG002_FP.png"
a=8;b=10;plot "$dir/2.all" using a:b with linespoints title "HG002-FP vs FN" lw 2 ps .6,

set xlabel "GTr"
set output "$dir/HG002_GT.png"
a=9;b=10;plot "$dir/2.all" using a:b with linespoints title "HG002-GT vs FN" lw 2 ps .6,
EOF

gnuplot <<EOF
set xlabel "FPr"
set ylabel "FNr"
set yrange [0:10]
set title "indels"

set terminal png size 800,600
set output "$dir/HG003_FP.png"

a=8;b=10;plot "$dir/3.all" using a:b with linespoints title "HG003-FP vs FN" lw 2 ps .6,

set xlabel "GTr"
set output "$dir/HG003_GT.png"
a=9;b=10;plot "$dir/3.all" using a:b with linespoints title "HG003-GT vs FN" lw 2 ps .6,
EOF

gnuplot <<EOF
set xlabel "FPr"
set ylabel "FNr"
set yrange [0:10]
set title "indels"

set terminal png size 800,600
set output "$dir/HG004_FP.png"

a=8;b=10;plot "$dir/4.all" using a:b with linespoints title "HG004-FP vs FN" lw 2 ps .6,

set xlabel "GTr"
set output "$dir/HG004_GT.png"
a=9;b=10;plot "$dir/4.all" using a:b with linespoints title "HG004-GT vs FN" lw 2 ps .6,
EOF

# Example gnuplot:
# set xlabel "FPr"
# set ylabel "FNr"
# set yrange [0:10]
# set title "PacBio 50x indels"

# FP vs FN
# a=8;b=10;t="ALL";plot "<grep ".t." dev/plot" using a:b with linespoints title "dev" lw 2 ps .6, "<grep ".t." cns/plot" using a:b with linespoints title "indels-cns" lw 2 ps .6

# GT vs FN
# a=9;b=10; ... as above
