#!/usr/bin/perl -w

# Reports the cumulative total of false positives, genotype assignment errors,
# and false negatives.  These are using an assumed filter of QUAL>=x, so as we
# increase QUAL filtering we're removing errors (FP and GT) and increasing
# missing calls (FN).
#
# We also report these in percentage terms.
# FP is fraction of called variants that are false
# GT is fraction of true variants with the wrong genotype
# FN is fraction of true variants not called

use strict;

my $dir=shift(@ARGV);
my $type=shift(@ARGV);

my ($fn_count, $I_fn_count, $D_fn_count) = (0,0,0);
my ($fp_count, $I_fp_count, $D_fp_count) = (0,0,0);
my ($gt_count, $I_gt_count, $D_gt_count) = (0,0,0);
my (@fp, @I_fp, @D_fp);
my (@tp, @I_tp, @D_tp);
my (@gt, @I_gt, @D_gt);
my ($total_true, $I_total_true, $D_total_true) = (0,0,0);

# False negatives.  No score, so just count them as base number
open(my $fn, "bcftools query -i 'TYPE=\"$type\"' -f '%QUAL %REF %ALT\n' $dir/0000.vcf|");
while (<$fn>) {
    my ($q,$r,$a) = split(/\s+/,$_);
    $fn_count++;
    $I_fn_count++ if (length($r) < length($a));
    $D_fn_count++ if (length($r) > length($a));
}
close($fn);

# False positives, bin by QUAL
open(my $fp, "bcftools query -i 'TYPE=\"$type\"' -f '%QUAL %REF %ALT\n' $dir/0001.vcf|");
while (<$fp>) {
    my ($q,$r,$a) = split(/\s+/,$_);
    @fp[int($q)]++;
    $fp_count++;
    if (length($r) < length($a)) {
	@I_fp[int($q)]++;
	$I_fp_count++;
    }
    if (length($r) > length($a)) {
	@D_fp[int($q)]++;
	$D_fp_count++;
    }
}
close($fp);

# True positives, bin by QUAL
my $total_tp_call = 0;
my $I_total_tp_call = 0;
my $D_total_tp_call = 0;
open(my $tp, "bcftools query -i 'TYPE=\"$type\"' -f '%QUAL %REF %ALT %FILTER\n' $dir/0003b.vcf|");
while (<$tp>) {
    my ($q,$r,$a,$f) = split(/\s+/,$_);
    $total_tp_call++;
    @tp[int($q)]++;
    if (length($r) < length($a)) {
	$I_total_tp_call++;
	@I_tp[int($q)]++;
    }
    if (length($r) > length($a)) {
	$D_total_tp_call++;
	@D_tp[int($q)]++;
    }

    if (/GTFAIL/) {
	@gt[int($q)]++;
	$gt_count++;
	if (length($r) < length($a)) {
	    @I_gt[int($q)]++;
	    $I_gt_count++;
	}
	if (length($r) > length($a)) {
	    @D_gt[int($q)]++;
	    $D_gt_count++;
	}
    }
}
close($tp);
$total_true = $fn_count + $total_tp_call;
$I_total_true = $I_fn_count + $I_total_tp_call;
$D_total_true = $D_fn_count + $D_total_tp_call;

print STDERR "Total true indel: $total_true\n";
print STDERR "Total true ins:   $I_total_true\n";
print STDERR "Total true del::  $D_total_true\n";

print STDERR "    QUAL :    FP    GT    FN\n";

for (my $qual = 0; $qual < 5000; $qual++) {
    $fp_count -= $fp[$qual] if (defined($fp[$qual]));
    $fn_count += $tp[$qual] if (defined($tp[$qual]));
    $gt_count -= $gt[$qual] if (defined($gt[$qual]));

    last if $fp_count == 0;

    $I_fp_count -= $I_fp[$qual] if (defined($I_fp[$qual]));
    $I_fn_count += $I_tp[$qual] if (defined($I_tp[$qual]));
    $I_gt_count -= $I_gt[$qual] if (defined($I_gt[$qual]));

    $D_fp_count -= $D_fp[$qual] if (defined($D_fp[$qual]));
    $D_fn_count += $D_tp[$qual] if (defined($D_tp[$qual]));
    $D_gt_count -= $D_gt[$qual] if (defined($D_gt[$qual]));

    my $total_call = $total_true - $fn_count + $fp_count;
    my $I_total_call = $I_total_true - $I_fn_count + $I_fp_count;
    my $D_total_call = $D_total_true - $D_fn_count + $D_fp_count;

    printf("ALL %4d : %5d %5d %5d", $qual, $fp_count, $gt_count, $fn_count);
    printf(" : %7.4f %7.4f %7.4f   $total_call\n",
	100*$fp_count/$total_call,
	100*$gt_count/$total_true,
	100*$fn_count/$total_true);

    printf("INS %4d : %5d %5d %5d",
	   $qual, $I_fp_count, $I_gt_count, $I_fn_count);
    printf(" : %7.4f %7.4f %7.4f\n",
	100*$I_fp_count/$I_total_call,
	100*$I_gt_count/$I_total_true,
	100*$I_fn_count/$I_total_true);

    printf("DEL %4d : %5d %5d %5d",
	   $qual, $D_fp_count, $D_gt_count, $D_fn_count);
    printf(" : %7.4f %7.4f %7.4f\n",
	100*$D_fp_count/$D_total_call,
	100*$D_gt_count/$D_total_true,
	100*$D_fn_count/$D_total_true);
}

__END__

/nfs/users/nfs_j/jkb/work/samtools_master/bcftools/plot_isec.pl HG002.GRCh38.60x.RG.bcftools-10h.vcf.isec indel > _2
Total true indel: 45604
Total true ins:   22194
Total true del::  23410

# Plot FP (x) vs FN (y)
a=8;b=10;plot "<grep ALL _1" using a:b with linespoints

# Plot GT (x) vs FN (y)
a=9;b=10;plot "<grep ALL _1" using a:b with linespoints

eg
a=8;b=10;t="ALL";plot "<grep ".t." _1" using a:b with linespoints title "dev" lw 2 ps .6, "<grep ".t." _2" using a:b with linespoints title "new" lw 2 ps .6, "<grep ".t." _3" using a:b with linespoints title "indels-2.0" lw 2 ps .6

