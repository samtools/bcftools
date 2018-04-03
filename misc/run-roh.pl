#!/usr/bin/env perl
#
#  The MIT License
#  
#  Copyright (c) 2017 Genome Research Ltd.
#  
#  Author: Petr Danecek <pd3@sanger.ac.uk>
#  
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#  
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#  
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.

use strict;
use warnings;
use Carp;
use Data::Dumper;

my $opts = parse_params();
run_roh($opts);
eval_roh($opts);

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg,"\n"; }
    print 
        "About: This is a convenience wrapper for \"bcftools roh\" which takes multiple VCF/BCF files\n",
        "       and locates regions private to a sample or shared across multiple samples. On input it\n",
        "       expects a directory with .vcf, .vcf.gz or .bcf files, a file with allele frequencies\n",
        "       and optionally a genetic map. See http://samtools.github.io/bcftools/howtos/roh-calling.html\n",
        "       for details\n",
        "Usage: run-roh.pl [OPTIONS]\n",
        "Options:\n",
        "   -a, --af-annots <file>      Allele frequency annotations [1000GP-AFs/AFs.tab.gz]\n",
        "   -i, --indir <dir>           Input directory with VCF files\n",
        "       --include <expr>        Select sites for which the expression is true\n",
        "       --exclude <expr>        Exclude sites for which the epxression is true\n",
        "   -l, --min-length <num>      Filter input regions shorter than this [1e6]\n",
        "   -m, --genmap <dir>          Directory with genetic map in IMPUTE2 format (optional)\n",
        "   -M, --rec-rate <float>      constant recombination rate per bp (optional)\n",
        "   -n, --min-markers <num>     Filter input regions with fewer marker than this [100]\n",
        "   -o, --outdir <dir>          Output directory\n",
        "   -q, --min-qual <num>        Filter input regions with quality smaller than this [10]\n",
        "       --roh-args <string>     Extra arguments to pass to bcftools roh\n",
        "   -s, --silent                Quiet output, do not print commands\n",
        "   -h, -?, --help              This help message\n",
        "\n";
    exit -1;
}
sub parse_params
{
    my $opts =
    {
        af_annots   => '1000GP-AFs/AFs.tab.gz', 
        verbose     => 1,
        min_length  => 1e6,
        min_markers => 100,
        min_qual    => 10,
        roh_args    => '',
    };
    while (defined(my $arg=shift(@ARGV)))
    {
        if (                 $arg eq '--roh-args' ) { $$opts{roh_args}=shift(@ARGV); next }
        if (                 $arg eq '--include' ) { $$opts{include_expr}=shift(@ARGV); next }
        if (                 $arg eq '--exclude' ) { $$opts{exclude_expr}=shift(@ARGV); next }
        if ( $arg eq '-q' || $arg eq '--min-qual' ) { $$opts{min_qual}=shift(@ARGV); next }
        if ( $arg eq '-l' || $arg eq '--min-length' ) { $$opts{min_length}=shift(@ARGV); next }
        if ( $arg eq '-n' || $arg eq '--min-markers' ) { $$opts{min_markers}=shift(@ARGV); next }
        if ( $arg eq '-s' || $arg eq '--silent' ) { $$opts{verbose}=0; next }
        if ( $arg eq '-a' || $arg eq '--af-annots' ) { $$opts{af_annots}=shift(@ARGV); next }
        if ( $arg eq '-m' || $arg eq '--genmap' ) { $$opts{genmap}=shift(@ARGV); next }
        if ( $arg eq '-M' || $arg eq '--rec-rate' ) { $$opts{rec_rate}=shift(@ARGV); next }
        if ( $arg eq '-o' || $arg eq '--outdir' ) { $$opts{outdir}=shift(@ARGV); next }
        if ( $arg eq '-i' || $arg eq '--indir' ) { $$opts{indir}=shift(@ARGV); next }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    if ( !exists($$opts{outdir}) ) { error("Missing the -o, --outdir option.\n") }
    if ( !exists($$opts{indir}) ) { error("Missing the -i, --indir option.\n") }
    if ( ! -e $$opts{af_annots} ) { error("The annotation file does not exist: $$opts{af_annots}\n"); }
    if ( ! -e "$$opts{af_annots}.tbi" ) { error("The annotation file is not indexed: $$opts{af_annots}.tbi\n"); }
    if ( ! -e "$$opts{af_annots}.hdr" ) { error("The annotation file has no header: $$opts{af_annots}.hdr\n"); }
    if ( exists($$opts{genmap}) && ! -d "$$opts{genmap}" ) { error("The directory with genetic maps does not exist: $$opts{genmap}\n"); }
    if ( exists($$opts{include_expr}) ) { $$opts{include_expr} =~ s/\'/\'\\\'\'/g; $$opts{inc_exc} .= qq[ -i '$$opts{include_expr}']; }
    if ( exists($$opts{exclude_expr}) ) { $$opts{exclude_expr} =~ s/\'/\'\\\'\'/g; $$opts{inc_exc} .= qq[ -e '$$opts{exclude_expr}']; }
    return $opts;
}

sub cmd
{
    my ($cmd,%args) = @_;

    if ( $args{verbose} ) { print STDERR $cmd,"\n"; }

    # Why not to use backticks? Perl calls /bin/sh, which is often bash. To get the correct
    #   status of failing pipes, it must be called with the pipefail option.

    my $kid_io;
    my $pid = open($kid_io, "-|");
    if ( !defined $pid ) { error("Cannot fork: $!"); }

    my @out;
    if ($pid) 
    {
        # parent
        @out = <$kid_io>;
        close($kid_io);
    } 
    else 
    {      
        # child
        exec('/bin/bash', '-o','pipefail','-c', $cmd) or error("Failed to run the command [/bin/sh -o pipefail -c $cmd]: $!");
    }

    if ( exists($args{exit_on_error}) && !$args{exit_on_error} ) { return @out; }

    my $exit_status = $?;
    my $status = exists($args{require_status}) ? $args{require_status} : 0;
    if ( $status ne $exit_status ) 
    {
        my $msg;
        if ( $? & 0xff )
        {
            $msg = "The command died with signal ".($? & 0xff);
        }
        else
        {
            $msg = "The command exited with status ".($? >> 8)." (expected $status)";
        }
        $msg .= ":\n\t$cmd\n\n";
        if ( @out ) {  $msg .= join('',@out,"\n\n"); }
        error($msg); 
    }
    return @out;
}

# determine the common prefix/suffix in file names like:
#   genetic_map_chr12_combined_b37.txt
#   genetic_map_chr13_combined_b37.txt
sub parse_genmap_path
{
    my ($opts) = @_;
    if ( !exists($$opts{genmap}) or !-d $$opts{genmap} ) { return ''; }
    my @files  = glob("$$opts{genmap}/*");
    my $prefix = $files[0];
    my $suffix = $files[0];
    for my $file (@files)
    {
        while ( length($prefix) && index($file,$prefix)==-1 )
        {
            substr($prefix,length($prefix)-1,1,'');
        }
        if ( !length($prefix) ) { last; }
    }
    for my $file (@files)
    {
        while ( length($suffix) && rindex($file,$suffix)==-1 )
        {
            substr($suffix,0,1,'');
        }
        if ( !length($suffix) ) { last; }
    }
    my @test = glob("$prefix*$suffix");
    if ( @test != @files ) 
    { 
        print STDERR "Warning: Could not determine the genetic map files [$prefix][$suffix]\n";
        return '';
    }
    return "-m $prefix\{CHROM}$suffix";
}

sub run_roh
{
    my ($opts) = @_;
    cmd("mkdir -p $$opts{outdir}",%$opts) unless -d $$opts{outdir};

    my $chr_fname = "$$opts{outdir}/chr-names.txt";
    open(my $fh,'>',$chr_fname) or error("$chr_fname: $!");
    for my $chr (1..22,'X') { print $fh "chr$chr\t$chr\n"; }
    close($fh) or error("close failed: $chr_fname");

    my @files = ();
    opendir(my $dh,$$opts{indir}) or error("$$opts{indir}: $!");
    while (my $file = readdir($dh))
    {
        if ( !($file=~/\.vcf$/i) && !($file=~/\.vcf\.gz$/i) && !($file=~/\.bcf$/i) ) { next; }
        my $outfile = "$$opts{outdir}/$`.bcf";
        push @files,$outfile;
        if ( -e $outfile ) { next; }

        my $cmd = 
            "bcftools annotate --rename-chrs $chr_fname '$$opts{indir}/$file' -Ou | " .
            "bcftools annotate -c CHROM,POS,REF,ALT,AF1KG -h $$opts{af_annots}.hdr -a $$opts{af_annots} ";

        if ( exists($$opts{inc_exc}) )
        {
            $cmd .= " -Ou | bcftools view $$opts{inc_exc} ";
        }
        $cmd .= "-Ob -o $outfile.part && ";
        $cmd .= "mv $outfile.part $outfile";

        cmd($cmd, %$opts);
    }
    closedir($dh) or error("close failed: $$opts{indir}");

    my $genmap = parse_genmap_path($opts);
    if ( exists($$opts{rec_rate}) ) { $genmap .= " -M $$opts{rec_rate}"; }

    for my $file (@files)
    {
        if ( -e "$file.txt.gz" ) { next; }
        my @out = cmd("bcftools roh $$opts{roh_args} --AF-tag AF1KG $genmap $file -Orz -o $file.txt.gz.part 2>&1 | tee -a $file.log",%$opts);
        for my $line (@out)
        {
            if ( !($line=~m{total/processed:\s+(\d+)/(\d+)}) ) { next; }
            my $total = $1;
            my $used  = $2;
            if ( !$total or $used/$total < 0.3 )
            {
                print STDERR @out;
                print STDERR "WARNING: Less than 30% of sites was used!\n\n"; 
            }
        }
        cmd(qq[bcftools query -f'GT\\t%CHROM\\t%POS[\\t%SAMPLE\\t%GT]\\n' $file | gzip -c >> $file.txt.gz.part && mv $file.txt.gz.part $file.txt.gz],%$opts);
    }
    $$opts{files} = \@files;
}

sub next_region
{
    my ($regions) = @_;
    my $chr = undef;
    for my $chrom (sort keys %$regions)
    {
        if ( %{$$regions{$chrom}} ) { $chr = $chrom; last; }
        delete($$regions{$chrom});
    }
    if ( !defined $chr ) { return undef; }

    my %min = ();
    for my $smpl (keys %{$$regions{$chr}})
    {
        if ( !exists($$regions{$chr}{$smpl}) ) { next; }
        my $reg = $$regions{$chr}{$smpl}[0];
        if ( !%min ) 
        {
            $min{chr} = $chr;
            $min{beg} = $$reg{beg};
            $min{end} = $$reg{end};
            next;
        }
        if ( $min{beg} > $$reg{beg} ) { $min{beg} = $$reg{beg}; }
    }
    if ( !%min ) { return undef; }
    for my $smpl (keys %{$$regions{$chr}})
    {
        if ( !exists($$regions{$chr}{$smpl}) ) { next; }
        my $reg = $$regions{$chr}{$smpl}[0];
        if ( $min{end} > $$reg{end} ) { $min{end} = $$reg{end}; }
        if ( $min{end} > $$reg{beg} - 1 && $min{beg} != $$reg{beg} ) { $min{end} = $$reg{beg} - 1; }
    }
    return \%min;
}

sub eval_roh
{
    my ($opts) = @_;
    my %regions = ();
    my %samples = ();
    my %lengths = ();
    for my $file (@{$$opts{files}})
    {
        open(my $fh,"gunzip -c $file.txt.gz |") or error("gunzip -c $file.txt.gz: $!");
        while (my $line=<$fh>)
        {
            if ( !($line=~/^RG/) ) { next; }
            my (%vals);
            @vals{qw(type smpl chr beg end len num qual)} = split(/\s+/,$line);
            if ( $vals{type} ne 'RG' ) { next; }
            chomp($vals{qual});
            if ( $vals{len} < $$opts{min_length} ) { next; }
            if ( $vals{num} < $$opts{min_markers} ) { next; }
            if ( $vals{qual} < $$opts{min_qual} ) { next; }
            push @{$regions{$vals{chr}}{$vals{smpl}}},\%vals;
            $samples{$vals{smpl}} = 1;
            $lengths{$vals{smpl}} += $vals{end} - $vals{beg} + 1;
        }
        close($fh) or error("close failed: gunzip -c $file.txt.gz");
    }
    open(my $fh,'>',"$$opts{outdir}/merged.txt") or error("$$opts{outdir}/merged.txt: $!");
    my @samples = sort keys %samples;
    print $fh "# [1]chrom\t[2]beg\t[3]end\t[4]length (Mb)";
    my $i = 5;
    for my $smpl (@samples) { print $fh "\t[$i]$smpl"; $i++; } 
    print $fh "\n";
    while (my $min = next_region(\%regions))
    {
        my $chr = $$min{chr};
        my $beg = $$min{beg};
        my $end = $$min{end};
        printf $fh "$chr\t$beg\t$end\t%.2f",($end-$beg+1)/1e6;
        for my $smpl (@samples)
        {
            if ( !exists($regions{$chr}{$smpl}) ) { goto not_present; }
            my $reg = $regions{$chr}{$smpl}[0];
            if ( $$reg{beg} > $end ) { goto not_present; }
            if ( $$reg{end} > $end ) { $regions{$chr}{$smpl}[0]{beg} = $end + 1; }
            else { shift @{$regions{$chr}{$smpl}}; }
            if ( !@{$regions{$chr}{$smpl}} ) { delete($regions{$chr}{$smpl}); }
            $lengths{$smpl} -= $end - $beg + 1;
            print $fh "\t1";
            next;
not_present:
            print $fh "\t0";
        }
        print $fh "\n";
    }
    close($fh) or error("close failed: $$opts{outdir}/merged.txt");
    for my $smpl (@samples)
    {
        if ( $lengths{$smpl}!=0 )
        {
            print STDERR "ERROR: a bug detected, sanity check failed, expected zero length : $smpl .. $lengths{$smpl}\n"; 
        }
    }
    print STDERR "The merged regions are in $$opts{outdir}/merged.txt\n";
}


