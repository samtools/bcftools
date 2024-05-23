#!/bin/bash
#
# Usage compare_vcf A.vcf B.vcf [region-skip.bed]

bcftools=${BCFTOOLS:-bcftools}
vt=${VT:-/nfs/users/nfs_j/jkb/lustre/vt/vt}

if [ $# -lt 2 ]
then
    echo Usage: compare_vcf A.vcf B.vcf [region-skip.bed] [region-include.bed] [region]
    exit 1
fi

v1=$1
v2=$2
exclude=$3
include=$4
region=$5

qual=${QUAL:-30}

#href=${HREF:-/nfs/srpipe_references/references/Human/1000Genomes_hs37d5/all/fasta/hs37d5.fa} 
href=${REF:-$HREF38}

pp() {
    awk "END {printf(\"%f\", 100 *$1 / $2)}" < /dev/null
}

norm() {
    # Also consider norm "-m +both" without a "-d both" step.
    # This produces more differences, but is also more correct.

    v=$1
    if [ x"$region" != x ]
    then
	#$bcftools norm -t $region -f $href $v 2>/dev/null | $bcftools norm -d both -N | $bcftools view -T ^$exclude | $bcftools view -T $include > $v.norm.vcf
	#$bcftools norm -m -both -t $region -f $href $v 2>/dev/null | $bcftools view -T ^$exclude | $bcftools view -T $include > $v.norm.vcf
	if [ x"$exclude" != x ]
	then
	    $bcftools norm -m -both -t $region -f $href $v 2>/dev/null | $vt decompose_blocksub - 2>/dev/null | $bcftools view -T ^$exclude | $bcftools view -T $include > $v.norm.vcf
	else
	    $bcftools norm -m -both -t $region -f $href $v 2>/dev/null | $vt decompose_blocksub - 2>/dev/null | $bcftools view -T $include > $v.norm.vcf
	fi
    elif [ x"$include" != x ]
    then
	$bcftools norm -m -both -f $href $v 2>/dev/null | $vt decompose_blocksub - 2>/dev/null | $bcftools view -T ^$exclude | $bcftools view -T $include > $v.norm.vcf
    elif [ x"$exclude" != x ]
    then
	$bcftools norm -m -both -f $href $v 2>/dev/null | $vt decompose_blocksub - 2>/dev/null | $bcftools view -T ^$exclude > $v.norm.vcf
    else
	$bcftools norm -m -both -f $href $v 2>/dev/null | $vt decompose_blocksub - 2>/dev/null > $v.norm.vcf
    fi
    rm $v.norm.vcf.gz 2>/dev/null
    bgzip $v.norm.vcf
    $bcftools index $v.norm.vcf.gz
}

if [ "$NORM" != "" ]
then
    norm $v1
fi
norm $v2

# NB: consider CHM13_1 chr20:2241427
# The Illumina data is all mapped with high mqual and shows homozygous deletion
# of AAAC in an AAAC STR.  The truth set claims AAAC del is heterozygous.
#
# $Bcftools isec therefore claims it is shared between both sets (in 0002 and 0003.vcf)
# but also only in the query set (0001) due to the extra allele.

#$bcftools isec -p $v1.isec $v1.norm.vcf.gz $v2.norm.vcf.gz
$bcftools isec -c both -p $v2.isec $v1.norm.vcf.gz $v2.norm.vcf.gz
# Or "-c any"? Works better possibly


# 0000 is private to v1 => FN
# 0001 is private to v2 => FP
# 0002/3 are records common to v1/v2 (from v1 or v2 respectively => TP)

# Depth filtering is beneficial to all tools, so we use it in this evaluation.
case $v2 in
    *15x*)  dp=15;DP=30;;
    *30x*)  dp=30;DP=60;;
    *53X*)  dp=53;DP=106;;
    *60x*)  dp=60;DP=120;;
    *100x*) dp=100;DP=200;;
    *150x*) dp=150;DP=300;;
    *300x*) dp=300;DP=600;;
    *)     dp=90;DP=120;;
#    *)     DP=90;;
esac

case $v2 in
    *gatk*)
	prog=GATK
        #gatk: https://software.broadinstitute.org/gatk/documentation/article.php?id=3225
	# s_filt_exp="QUAL < $qual || QD < 2 || FS > 60 || MQ < 40 || SOR > 3 || MQRankSum < -12.5 || ReadPosRankSum < -8"
	# i_filt_exp="QUAL < $qual || QD < 2 || FS > 200 || ReadPosRankSum < -20"
	s_filt_exp="QUAL < $qual || QD < 2 || FS > 60 || MQ < 40 || SOR > 3 || MQRankSum < -2.5 || ReadPosRankSum < -8 || INFO/DP>$DP"
	i_filt_exp="QUAL < $qual || QD < 2 || FS > 60 || ReadPosRankSum < -8 || INFO/DP>$DP || SOR > 3"

	# better on GIAB.  Also on others?  Unknown.  SOR 3->6 and MQ 40->20
	s_filt_exp="QUAL < $qual || QD < 2 || FS > 60 || MQ < 20 || SOR > 6 || MQRankSum < -4 || ReadPosRankSum < -8 || INFO/DP>$DP"

	# Depth aware variant of above
	s_filt_exp="QUAL < $qual || QD < 2 || FS > 55+INFO/DP/6 || MQ < 20 || SOR > 6+INFO/DP/25 || MQRankSum < -(3+INFO/DP/50) || ReadPosRankSum < -(3.5+INFO/DP/50) || INFO/DP>$DP"

	;;
    *freebayes*)
        prog=Freebayes
        #freebayes: https://wiki.uiowa.edu/download/attachments/145192256/erik%20garrison%20-%20iowa%20talk%202.pdf?api=v2
	#s_filt_exp="QUAL < $qual || SAF <= 0 || SAR <= 0 || RPR <= 1 || RPL <= 1 || INFO/DP > $DP"
	s_filt_exp="QUAL < $qual || SAF <= 0 || SAR <= 0 || RPR <= 0 || RPL <= 0 || INFO/DP > $DP"

	#i_filt_exp=$s_filt_exp
	#i_filt_exp="QUAL < $qual || INFO/DP > $DP"
	i_filt_exp="QUAL < $qual || RPR <= 0 || RPL <= 0 || INFO/DP > $DP"
	;;
    *bcftools*)
        prog=Bcftools

	# Simple filters; QUAL, DP and for indels IDV/IMF
	s_filt_exp="QUAL < $qual || DP>$DP"
	i_filt_exp="IDV < 3 || IMF < 0.1 || DP>$DP || QUAL < $qual"

	# Suggest:
	#s_filt_exp="QUAL < $qual || DP>$dp*2 || MQBZ < -(4+$dp/30) || RPBZ > (3.5+$dp/60) || RPBZ < -(3.5+$dp/60) || FORMAT/SP > 30+$dp/2 || SCBZ > 3+$dp/50"
	#q_filt_exp="IDV<(2+DP*$qual/2000) || IMF < 0.02+(($qual+1)/($qual+31))*(($qual+1)/($qual+31))/3 || DP>$dp*2 || MQBZ < -(4+$dp/30) || FORMAT/SP > 30+$dp/2 || RPBZ+SCBZ > 5.5"

	# I_FILT_EXP='IDV<(2+DP*$qual/2000) || IMF < 0.02+(($qual+1)/($qual+31))*(($qual+1)/($qual+31))/3 || DP>$dp*2 || MQBZ < -(4+$dp/30) || FORMAT/SP > 30+$dp/2 || RPBZ+SCBZ > 5.5' S_FILT_EXP='QUAL < $qual || DP>$dp*2 || MQBZ < -(4+$dp/30) || RPBZ > (3.5+$dp/60) || RPBZ < -(3.5+$dp/60) || FORMAT/SP > 30+$dp/2 || SCBZ > 3+$dp/50' ./compare_vcf_simple.sh

	;;
    *octopus*)
        prog=octopus
	s_filt_exp="FILTER!=\"PASS\" || QUAL < $qual || INFO/DP>$DP"
	i_filt_exp="FILTER!=\"PASS\" || QUAL < $qual || INFO/DP>$DP"
	;;
    *)
        prog=Unknown
	echo Unrecognised type, no specific filters
	s_filt_exp="QUAL<$qual || INFO/DP>$DP"
	i_filt_exp="QUAL<$qual || INFO/DP>$DP"
	;;
esac

S_FILT_EXP=`eval echo \"${S_FILT_EXP:-$s_filt_exp}\"`
I_FILT_EXP=`eval echo \"${I_FILT_EXP:-$i_filt_exp}\"`

#echo "SNP:   $S_FILT_EXP" 1>&2
#echo "Indel: $I_FILT_EXP" 1>&2

# We classify wrong genotypes as both a true and untrue call.  We don't
# label it as a false negative as the call exists and it's not as bad
# a situation as missing it entirely.

# Hacky addition of GTFAIL filter if 0002.vcf and 0003.vcf have a different
# GT call.
perl -e 'BEGIN {$"="\t"} open(F1,"<'$v2'.isec/0002.vcf");open(F2,"<'$v2'.isec/0003.vcf");$_=<F2>;print;print "##FILTER=<ID=GTFAIL,Description=\"GT mismatch\">\n";while(<F1>) {next if /^#/;chomp($_);@F=split("\t",$_);$F[-1]=~/^(\d+)[|\/](\d+)/;$gt=($1<$2)?"$1/$2":"$2/$1";while (<F2>) {if (/^#/) {print;next} else {last}};chomp($_);@G=split("\t",$_);$G[-1]=~/^(\d+)[|\/](\d+)/;$gt2=($1<$2)?"$1/$2":"$2/$1";if ($gt ne $gt2) {$G[6]="GTFAIL"} print "@G\n"}' > $v2.isec/0003b.vcf

# Produce isec/filtered.vcf as a filtered copy of the call set, so
# we can pass this to rtg vcfeval for an alternative way of evaluating
# data sets.
$bcftools view -e "QUAL < $qual || (TYPE='snp' && ($S_FILT_EXP)) || (TYPE='indel' && ($I_FILT_EXP))" $v2 > $v2.isec/filtered.vcf

# QUAL 1 is recommended minimum for freebayes to remove detritus.
v1_snp=`    $bcftools view -H -i "TYPE='snp'" $v2.isec/0000.vcf|wc -l`
v2_snp=`    $bcftools view -H -i "TYPE='snp' && QUAL >= 1" $v2.isec/0001.vcf|wc -l`
v2_snp_hq=` $bcftools view -H -i "TYPE='snp' && QUAL >= $qual" $v2.isec/0001.vcf|wc -l`
v2_snp_fi=` $bcftools view    -i "TYPE='snp'" $v2.isec/0001.vcf | bcftools view -H -e "$S_FILT_EXP" -|wc -l`
v12_snp=`   $bcftools view -H -i "TYPE='snp'" $v2.isec/0002.vcf|wc -l`
v12_snp_hq=`$bcftools view -H -i "TYPE='snp' && QUAL >= $qual" $v2.isec/0003.vcf|wc -l`
v12_snp_fi=`$bcftools view    -i "TYPE='snp'" $v2.isec/0003.vcf | bcftools view -H -e "$S_FILT_EXP" -|wc -l`

v12_snp_gt=`   $bcftools view -H -f GTFAIL -i "TYPE='snp' && QUAL >= 1" $v2.isec/0003b.vcf|wc -l`
v12_snp_hq_gt=`$bcftools view -H -f GTFAIL -i "TYPE='snp' && QUAL >= $qual" $v2.isec/0003b.vcf|wc -l`
v12_snp_fi_gt=`$bcftools view    -f GTFAIL -i "TYPE='snp'" $v2.isec/0003b.vcf | bcftools view -H -e "$I_FILT_EXP" -|wc -l`

#v2_snp=`expr $v2_snp + $v12_snp_gt`
#v2_snp_hq=`expr $v2_snp_hq + $v12_snp_hq_gt`
#v2_snp_fi=`expr $v2_snp_fi + $v12_snp_fi_gt`


v1_indel=`    $bcftools view -H -i "TYPE='indel'" $v2.isec/0000.vcf|wc -l`
v2_indel=`    $bcftools view -H -i "TYPE='indel' && QUAL >= 1" $v2.isec/0001.vcf|wc -l`
v2_indel_hq=` $bcftools view -H -i "TYPE='indel' && QUAL >= $qual" $v2.isec/0001.vcf|wc -l`
v2_indel_fi=` $bcftools view    -i "TYPE='indel'" $v2.isec/0001.vcf | bcftools view -H -e "$I_FILT_EXP" -|wc -l`
v12_indel=`   $bcftools view -H -i "TYPE='indel'" $v2.isec/0002.vcf|wc -l`
v12_indel_hq=`$bcftools view -H -i "TYPE='indel' && QUAL >= $qual" $v2.isec/0003.vcf|wc -l`
v12_indel_fi=`$bcftools view    -i "TYPE='indel'" $v2.isec/0003.vcf | bcftools view -H -e "$I_FILT_EXP" -|wc -l`

v12_indel_gt=`   $bcftools view -H -f GTFAIL -i "TYPE='indel'" $v2.isec/0003b.vcf|wc -l`
v12_indel_hq_gt=`$bcftools view -H -f GTFAIL -i "TYPE='indel' && QUAL >= $qual" $v2.isec/0003b.vcf|wc -l`
v12_indel_fi_gt=`$bcftools view    -f GTFAIL -i "TYPE='indel'" $v2.isec/0003b.vcf | bcftools view -H -e "$I_FILT_EXP" -|wc -l`

#v2_indel=`expr $v2_indel + $v12_indel_gt`
#v2_indel_hq=`expr $v2_indel_hq + $v12_indel_hq_gt`
#v2_indel_fi=`expr $v2_indel_fi + $v12_indel_fi_gt`

# quality trimmed FN aren't the records private to v1 above QUAL, but the
# total number of records not in v12 after filtering.  Thus as we increase
# acceptance threshold to reduce FP we increase FN.
v1_snp_hq=`expr $v1_snp + $v12_snp - $v12_snp_hq`
v1_snp_fi=`expr $v1_snp + $v12_snp - $v12_snp_fi`
v1_indel_hq=`expr $v1_indel + $v12_indel - $v12_indel_hq`
v1_indel_fi=`expr $v1_indel + $v12_indel - $v12_indel_fi`

# Assumption A.vcf is truth set and B.vcf is test set
if [ "$FORMAT" = "tex" ]
then
	printf '\\bigskip\n'
	printf '\\begin{minipage}{\\linewidth}\n'
	printf '\\centering\n'
	printf '\\captionof{table}{%s: FIXME}\n' $prog
	printf '{\\begin{tabular}{ll|r|rr}\n'
	printf 'Variants & & \\textbf{Q>0} & \\textbf{Q>=%d} & \\textbf{Filtered} \\\\ \\midrule\n' $qual
	printf 'SNP   & TP & %7d & %7d & %7d \\\\\n' $v12_snp  $v12_snp_hq  $v12_snp_fi
	printf 'SNP   & FP & %7d & %7d & %7d \\\\\n' $v2_snp   $v2_snp_hq   $v2_snp_fi
	printf 'SNP   & GT & %7d & %7d & %7d \\\\\n' $v12_snp_gt $v12_snp_hq_gt $v12_snp_fi_gt
	printf 'SNP   & FN & %7d & %7d & %7d \\\\\n' $v1_snp   $v1_snp_hq   $v1_snp_fi
	printf '\\midrule\n';
	printf 'InDel & TP & %7d & %7d & %7d \\\\\n' $v12_indel $v12_indel_hq $v12_indel_fi
	printf 'InDel & FP & %7d & %7d & %7d \\\\\n' $v2_indel  $v2_indel_hq  $v2_indel_fi
	printf 'InDel & GT & %7d & %7d & %7d \\\\\n' $v12_indel_gt $v12_indel_hq_gt $v12_indel_fi_gt
	printf 'InDel & FN & %7d & %7d & %7d \\\\\n' $v1_indel  $v1_indel_hq  $v1_indel_fi
	printf '\\end{tabular}}\n'
	printf '\\par\n'
	printf '\\textbf{CRAM qual size x,x}\n'
	printf '\\bigskip\n'
	printf '\\end{minipage}\n'

elif [ "$FORMAT" = "gnuplot" ]
then
    # v1_snp  = total number of FN SNPs (at any qual; from truth set)
    # v2_snp  = total number of FP SNPs (at any qual; our extra calls)
    # v12_snp = total number of TP SNPs
    # => v12_snp + v1_snp = total true SNPs
    #
    # Columns:
    # Filename
    # Qual $q
    # No. TP SNP in this 10-qual bin
    # No. FP SNP in this 10-qual bin
    # No. TP Indel in this 10-qual bin
    # No. FP Indel in this 10-qual bin
    # Total no. SNPs
    # Number FN SNP at QUAL >= $q
    # Number FP SNP at QUAL >= $q
    # Total no. INDELs
    # Number FN INDEL at QUAL >= $q
    # Number FP INDEL at QUAL >= $q
    # //Number GT SNP errs at QUAL >= $q
    tot_snp=`expr $v1_snp + $v12_snp`
    tot_indel=`expr $v1_indel + $v12_indel`

    # For gnuplot
    for Q in `seq 0 10 250`
    do
	qual=$Q
	qual_max=`expr $Q + 10`
	v2_snp_hq=` $bcftools view -H -i "TYPE='snp' && QUAL >= $qual && QUAL < $qual_max" $v2.isec/0001.vcf|wc -l`
	v2_snp_hq2=` $bcftools view -H -i "TYPE='snp' && QUAL >= $qual" $v2.isec/0001.vcf|wc -l`
	v12_snp_hq=`$bcftools view -H -i "TYPE='snp' && QUAL >= $qual && QUAL < $qual_max" $v2.isec/0003.vcf|wc -l`
	v12_snp_hq2=`$bcftools view -H -i "TYPE='snp' && QUAL >= $qual" $v2.isec/0003.vcf|wc -l`
	v2_indel_hq=` $bcftools view -H -i "TYPE='indel' && QUAL >= $qual && QUAL < $qual_max" $v2.isec/0001.vcf|wc -l`
	v2_indel_hq2=` $bcftools view -H -i "TYPE='indel' && QUAL >= $qual" $v2.isec/0001.vcf|wc -l`
	v12_indel_hq=`$bcftools view -H -i "TYPE='indel' && QUAL >= $qual && QUAL < $qual_max" $v2.isec/0003.vcf|wc -l`
	v12_indel_hq2=`$bcftools view -H -i "TYPE='indel' && QUAL >= $qual" $v2.isec/0003.vcf|wc -l`

	# Total number of SNPs minus true HQ snps we call
	v12_snp_hq2=`expr $tot_snp - $v12_snp_hq2`
	v12_indel_hq2=`expr $tot_indel - $v12_indel_hq2`
	printf "$v2\t$qual\t$v12_snp_hq $v2_snp_hq\t$v12_indel_hq $v2_indel_hq\t$tot_snp $v12_snp_hq2 $v2_snp_hq2\t$tot_indel $v12_indel_hq2 $v2_indel_hq2\n"
    done

elif [ "$FORMAT" = "percent" ]
then
	printf "SNP          Q>0 /   Q>=$qual / Filtered\n"
	x=`expr $v12_snp + $v1_snp`
	printf "SNP   TP %7.2f / %7.2f / %7.2f\n" \
	       `pp $v12_snp $x` `pp $v12_snp_hq $x` `pp $v12_snp_fi $x`
	x=`expr $v12_snp + $v2_snp`
	printf "SNP   FP %7.2f / %7.2f / %7.2f\n" \
	       `pp $v2_snp $x` `pp $v2_snp_hq $x` `pp $v2_snp_fi $x`
	x=`expr $v12_snp`
	printf 'SNP   GT %7.2f / %7.2f / %7.2f\n' \
	       `pp $v12_snp_gt $x` `pp $v12_snp_hq_gt $x` `pp $v12_snp_fi_gt $x`
	x=`expr $v12_snp + $v1_snp`
	printf "SNP   FN %7.2f / %7.2f / %7.2f\n" \
	       `pp $v1_snp $x` `pp $v1_snp_hq $x` `pp $v1_snp_fi $x`	
#printf "SNP   %4.1f%% prec, %4.1f%% rec\n" 100.0*$v12_snp_hq/($v12_snp_hq+$v2_snp_hq) 100.0*$v12_snp_hq/($v12_snp_hq+$v1_snp_hq);
	printf "\n";
	x=`expr $v12_indel + $v1_indel`
	printf "InDel TP %7.2f / %7.2f / %7.2f\n" \
	       `pp $v12_indel $x` `pp $v12_indel_hq $x` `pp $v12_indel_fi $x`
	x=`expr $v12_indel + $v2_indel`
	printf "InDel FP %7.2f / %7.2f / %7.2f\n" \
	       `pp $v2_indel $x` `pp $v2_indel_hq $x` `pp $v2_indel_fi $x`
	x=`expr $v12_indel`
	printf 'InDel GT %7.2f / %7.2f / %7.2f\n' \
	       `pp $v12_indel_gt $x` `pp $v12_indel_hq_gt $x` `pp $v12_indel_fi_gt $x`
	x=`expr $v12_indel + $v1_indel`
	printf "InDel FN %7.2f / %7.2f / %7.2f\n" \
	       `pp $v1_indel $x` `pp $v1_indel_hq $x` `pp $v1_indel_fi $x`

else
	printf "SNP          Q>0 /   Q>=$qual / Filtered\n"
	printf "SNP   TP %7d / %7d / %7d\n" $v12_snp  $v12_snp_hq  $v12_snp_fi
	printf "SNP   FP %7d / %7d / %7d\n" $v2_snp   $v2_snp_hq   $v2_snp_fi
	printf 'SNP   GT %7d / %7d / %7d\n' $v12_snp_gt $v12_snp_hq_gt $v12_snp_fi_gt
	printf "SNP   FN %7d / %7d / %7d\n" $v1_snp   $v1_snp_hq   $v1_snp_fi
	#printf "SNP   %4.1f%% prec, %4.1f%% rec\n" 100.0*$v12_snp_hq/($v12_snp_hq+$v2_snp_hq) 100.0*$v12_snp_hq/($v12_snp_hq+$v1_snp_hq);
	printf "\n";
	printf "InDel TP %7d / %7d / %7d\n" $v12_indel $v12_indel_hq $v12_indel_fi
	printf "InDel FP %7d / %7d / %7d\n" $v2_indel  $v2_indel_hq  $v2_indel_fi
	printf 'InDel GT %7d / %7d / %7d\n' $v12_indel_gt $v12_indel_hq_gt $v12_indel_fi_gt
	printf "InDel FN %7d / %7d / %7d\n" $v1_indel  $v1_indel_hq  $v1_indel_fi
	#printf "InDel %4.1f%% prec, %4.1f%% rec\n" 100.0*$v12_indel_hq/($v12_indel_hq+$v2_indel_hq) 100.0*$v12_indel_hq/($v12_indel_hq+$v1_indel_hq);
fi

#rm $v1.norm* $v2.norm*
#rm -rf $v2.isec
