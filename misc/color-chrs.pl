#!/usr/bin/env perl
#
#   Copyright (C) 2015, 2017 Genome Research Ltd.
#
#   Author: Petr Danecek <pd3@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

use strict;
use warnings;
use Carp;

my $opts = parse_params();
chr_dims($opts);
for my $file (@{$$opts{files}})
{
    read_dat($opts,$file);
}
merge_regs($opts);
plot($opts);

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print 
        "About: Plot output of \"bcftools +color-chrs\"\n",
        "Usage: color-chrs.pl [OPTIONS] output.dat\n",
        "Options:\n",
        "   -c, --colors <file>         File with list of \"chr hap color\".\n",
        "   -p, --prefix <name>         Prefix of output files.\n",
        "   -h, -?, --help              This help message.\n",
        "\n";
    exit -1;
}

sub parse_params
{
    my $opts =
    {
        #colors => ['red','green','blue','yellow'],
        colors => ['#ff0000','#008000','#0000ff','#ffff00'],
        height => 350,
        pad    => 10,
        dimL   => 300,     # arm length
        dimD   => 10,      # arm width
        dimE   => 7,       # arm pad
        dimB   => 5,       # arm end curve
    };
    $$opts{chr_width} = 2*$$opts{pad} + 2*$$opts{dimD} + $$opts{dimE};
    $$opts{width} = $$opts{chr_width}*23;
    while (defined(my $arg=shift(@ARGV)))
    {
        if ( -e $arg ) { push @{$$opts{files}},$arg; next }
        if ( $arg eq '-p' || $arg eq '--prefix' ) { $$opts{prefix}=shift(@ARGV); next }
        if ( $arg eq '-c' || $arg eq '--colors' ) { read_colors($opts,shift(@ARGV)); next; }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    if ( !exists($$opts{files}) ) { error("No files given?\n") }
    if ( !exists($$opts{prefix}) ) { error("Expected -p option\n") }
    return $opts;
}

sub read_colors
{
    my ($opts,$fname) = @_;
    open(my $fh,'<',$fname) or error("$fname: $!");
    while (my $line=<$fh>)
    {
        if ( $line=~/^\s*$/ ) { next; }
        $line =~ s/^\s*//;
        $line =~ s/\s*$//;
        my ($chr,$hap,$col) = split(/\s+/,$line);
        $$opts{hap_cols}{$chr}{$hap} = $col;
    }
    close($fh);
}

sub plot
{
    my ($opts) = @_;

    my $width  = $$opts{width};
    my $height = $$opts{height};
    my $pad    = $$opts{pad};

    my $svg = $$opts{svg} = [];
    push @$svg, q[<?xml version="1.0" encoding="UTF-8" standalone="yes"?>];
    push @$svg, q[<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">];
    push @$svg, qq[<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" height="100%" viewBox="0 0 $width $height" width="100%">];

    my $xpos = $pad;
    for my $chr (qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X))
    {
        draw_chr($opts,$chr,$xpos,$pad);
        $xpos += $$opts{chr_width};
    }
    $xpos = $pad + 10*$$opts{chr_width};
    for my $sample (keys %{$$opts{samples}})
    {
        $$opts{chrs}{$sample} = { len=>59e6, acen=>[24e6,26e6,28e6] };
        my $chr  = $$opts{chrs}{$sample};
        my $len1 = $$opts{dimL}*$$chr{acen}[0]/$$opts{max_chr_len};
        my $len2 = $$opts{dimL}*$$chr{acen}[2]/$$opts{max_chr_len};
        my $len3 = $$opts{dimL}*$$chr{len}/$$opts{max_chr_len};
        if ( exists($$opts{samples}{$sample}{1}) )
        {
            my $col = $$opts{samples}{$sample}{1};
            $$chr{px1} = [[0,$len1,{$col=>1.0}],[$len2,$len3,{$col=>1.0}]]; 
        }
        if ( exists($$opts{samples}{$sample}{2}) )
        {
            my $col = $$opts{samples}{$sample}{2};
            $$chr{px2} = [[0,$len1,{$col=>1.0}],[$len2,$len3,{$col=>1.0}]]; 
        }
        my $ypos = $$opts{dimL} - $len3 + $pad;
        draw_chr($opts,$sample,$xpos,$ypos);
        $xpos += 2*$$opts{chr_width};
    }

    push @$svg, q[</svg>];
    open(my $fh,'>',"$$opts{prefix}.svg");
    print $fh join("\n",@$svg);
    close($fh);
}

sub draw_chr
{
    my ($opts,$chr_name,$xpos,$ypos) = @_;

    my $svg = $$opts{svg};
    my $chr = $$opts{chrs}{$chr_name};
    my $dimL  = $$opts{dimL};
    my $dimC  = $dimL*($$chr{acen}[2]-$$chr{acen}[0])/$$opts{max_chr_len};
    my $dimL1 = $dimL*$$chr{acen}[0]/$$opts{max_chr_len};
    my $dimL2 = $dimL*($$chr{len} - $$chr{acen}[2])/$$opts{max_chr_len};
    my $dimD  = $$opts{dimD};   # arm width
    my $dimB  = $$opts{dimB};
    my $dimE  = $$opts{dimE};   # the gap width between arms
    my $dimE05 = $dimE*0.5;
    my $dimD05 = $dimD*0.5;
    my $dimB05 = $dimB*0.5;
    my $dimC05 = 0.5*$dimC;

    my $xtext = $xpos + $dimD + $dimE*0.5;
    my $ytext = $ypos;
    push @$svg, qq[<text text-anchor="middle" x="$xtext" y="$ytext">$chr_name</text>\n];

    $ypos += $$opts{pad};

    push @$svg, qq[
        <path d="M$xpos $ypos
                 l0 $dimL1
                 q$dimB $dimC05 0 $dimC
                 l0 $dimL2
                 q$dimD05 $dimB $dimD 0
                 l0 -$dimL2
                 q$dimE05 -$dimB $dimE 0 ];

    # bottom part of the right arm 
    push @$svg, qq[
                 l0 $dimL2
                 q$dimD05 $dimB $dimD 0
                 l0 -$dimL2
                 q-$dimB -$dimC05 0 -$dimC
                 l0 -$dimL1
                 q-$dimD05 -$dimB -$dimD 0
                 l0 $dimL1
                 q-$dimE05 $dimB -$dimE 0
                 l0 -$dimL1
                 q-$dimD05 -$dimB -$dimD 0
                 "
              style="stroke:#333; fill:#aaa;"/> ];

    if ( $$chr{px1} ) { draw_regs($opts,$chr_name,$$chr{px1},$xpos,$ypos); }
    if ( $$chr{px2} ) { draw_regs($opts,$chr_name,$$chr{px2},$xpos+$dimD+$dimE,$ypos); }
}

sub parse_hex_color
{
    my ($color) = @_;
    if ( !($color=~/^#/) ) { error("Could not parse: $color"); }
    $color =~ s/^#//;
    my @vals = split(//, $color);
    my ($r,$g,$b);
    if ( @vals==6 )
    {
        $r = hex($vals[0].$vals[1]);
        $g = hex($vals[2].$vals[3]);
        $b = hex($vals[4].$vals[5]);
    }
    elsif ( @vals==3 )
    {
        $r = hex($vals[0].$vals[0]);
        $g = hex($vals[1].$vals[1]);
        $b = hex($vals[2].$vals[2]);
    }
    else { error("Could not parse: $color\n"); }
    return ($r,$g,$b);
}

sub scale_color
{
    my ($color,$scale) = @_;
    my ($r1,$g1,$b1) = parse_hex_color($color);
    my ($r0,$g0,$b0) = parse_hex_color('#aaa');
    my $r = $scale*($r1-$r0) + $r0;
    my $g = $scale*($g1-$g0) + $g0;
    my $b = $scale*($b1-$b0) + $b0;
    return sprintf("#%0.2x%0.2x%0.2x",$r,$g,$b);
}

sub draw_regs
{
    my ($opts,$chr_name,$regs,$xpos,$ypos) = @_;

    my $svg = $$opts{svg};
    my $chr = $$opts{chrs}{$chr_name};
    my $dimL  = $$opts{dimL};
    my $dimC  = $dimL*($$chr{acen}[2]-$$chr{acen}[0])/$$opts{max_chr_len};
    my $dimL1 = $dimL*$$chr{acen}[0]/$$opts{max_chr_len};
    my $dimL2 = $$opts{dimL} - $dimL*$$chr{acen}[2]/$$opts{max_chr_len};
    my $dimD  = $$opts{dimD};   # arm width
    my $dimB  = $$opts{dimB};
    my $dimE  = $$opts{dimE};
    my $dimE05 = $dimE*0.5;
    my $dimD05 = $dimD*0.5;
    my $dimB05 = $dimB*0.5;
    my $dimC05 = 0.5*$dimC;

    for my $reg (@$regs)
    {
        my $sum = 0;
        my $cmax = undef;
        for my $c (keys %{$$reg[2]})
        {
            $sum += $$reg[2]{$c};
            if ( !defined $cmax or $$reg[2]{$cmax} < $$reg[2]{$c} ) { $cmax = $c; }
        }
        my $color = scale_color($cmax,$$reg[2]{$cmax}/$sum);
        my $y  = $ypos + $$reg[0];
        my $dy = $$reg[1] - $$reg[0] + 1;
        push @$svg, qq[
            <path d="M$xpos $y
                     l0 $dy
                     l$dimD 0
                     l0 -$dy
                     l-$dimD 0
                    "
                  style="stroke:$color;fill:$color;stroke-width:0;"/> ];
    }

}

sub hap2color
{
    my ($opts,$chr,$hap) = @_;
    if ( exists($$opts{hap_cols}{$chr}{$hap}) )
    {
        if ( !exists($$opts{hap_cols}{'*'}{$hap}) ) 
        {
            $$opts{hap_cols}{'*'}{$hap} = $$opts{hap_cols}{$chr}{$hap};
        }
        return $$opts{hap_cols}{$chr}{$hap};
    }
    elsif ( exists($$opts{hap_cols}{'*'}{$hap}) )
    {
        return $$opts{hap_cols}{'*'}{$hap};
    }
    if ( !exists($$opts{haps}) ) { $$opts{haps} = {}; }
    if ( !exists($$opts{haps}{$hap}) )
    {
        my $nhaps = scalar keys %{$$opts{haps}};
        my $color = $$opts{colors}[$nhaps];
        $$opts{haps}{$hap} = $color;
    }
    return $$opts{haps}{$hap};
}

sub read_dat
{
    my ($opts,$file) = @_;
    open(my $fh,'<',$file) or error("$file: $!");
    while (my $line=<$fh>)
    {
        if ( !($line=~/^SG/) ) { next; }
        my @items = split(/\s+/,$line);
        my $chr   = $items[1];
        my $start = $items[2];
        my $end   = $items[3];
        my $hap1  = $items[4];
        my $hap2  = $items[5];
        my $col1  = hap2color($opts,$chr,$hap1); 
        my $col2  = hap2color($opts,$chr,$hap2); 
        push @{$$opts{chrs}{$chr}{regs1}},[$start,$end,$col1];
        push @{$$opts{chrs}{$chr}{regs2}},[$start,$end,$col2];
        my ($smpl,$hap);
        ($smpl,$hap) = split(/:/,$hap1); $$opts{samples}{$smpl}{$hap} = $col1;
        ($smpl,$hap) = split(/:/,$hap2); $$opts{samples}{$smpl}{$hap} = $col2;
    }
    close($fh);
}

sub merge_regs
{
    my ($opts) = @_;
    for my $chr (keys %{$$opts{chrs}})
    {
        my $dat = $$opts{chrs}{$chr};
        if ( exists($$dat{regs1}) ) { $$dat{px1} = merge($opts,$dat,$$dat{regs1}); }
        if ( exists($$dat{regs2}) ) { $$dat{px2} = merge($opts,$dat,$$dat{regs2}); }
    }
}

sub merge
{
    my ($opts,$chr_dat,$regs) = @_;

    # merge adjacent regions of the same color
    for (my $i=1; $i<@$regs; $i++)
    {
        if ( $$regs[$i-1][2] eq $$regs[$i][2] ) 
        { 
            $$regs[$i-1][1] = $$regs[$i][1];
            splice(@$regs,$i,1);
            $i--;
            next; 
        }
    }
    for (my $i=0; $i<@$regs; $i++)
    {
        if ( $$regs[$i][1] < $$chr_dat{acen}[0] ) { next; }
        if ( $$regs[$i][0] > $$chr_dat{acen}[2] ) { next; }
        if ( $$regs[$i][0] >= $$chr_dat{acen}[0] && $$regs[$i][1] <= $$chr_dat{acen}[2] )
        {
            splice(@$regs,$i,1);
            $i--;
            next;
        }
        if ( $$regs[$i][1] < $$chr_dat{acen}[2] )
        {
            $$regs[$i][1] = $$chr_dat{acen}[0];
            next;
        }
        if ( $$regs[$i][0] > $$chr_dat{acen}[0] )
        {
            $$regs[$i][0] = $$chr_dat{acen}[2];
            next;
        }
        my $start = $$chr_dat{acen}[2];
        my $end   = $$regs[$i][1];
        my $color = $$regs[$i][2];
        $$regs[$i][1] = $$chr_dat{acen}[0];
        splice(@$regs,$i+1,0,[$start,$end,$color]);
    }

    # pixelize
    my @px = ([0,0,{}]);  # start,end,hash of contributions from colors
    my $dy = $$opts{max_chr_len}/$$opts{dimL};  # 1px covers $dy base pairs
    for my $reg (@$regs)
    {
        my $px_start = int($$reg[0] / $dy);
        my $px_end   = int($$reg[1] / $dy);
        my $color    = $$reg[2];
        my $contrib  = ($$reg[1] - $$reg[0])/$dy;

        if ( $px_start==$px_end )
        {
            if ( $px_end<=$px[-1][1] ) { $px[-1][2]{$color} += $contrib; }
            else { push @px, [$px_start,$px_end,{$color=>$contrib}]; }
            next;
        }

        my $rem_start = 1 - ($$reg[0]/$dy - $px_start);     # fractional remainders
        my $rem_end   = $$reg[1]/$dy - $px_end;

        if ( $px_start==$px[-1][1] )
        {
            $px[-1][2]{$color} += $rem_start;
        }

        if ( $px_end==$px[-1][1] )
        {
            $px[-1][2]{$color} += $rem_end;
            if ( $px_start+1<$px_end ) { error("really?!\n"); }
        }
        elsif ( $px_end>$px[-1][1] )
        {
            if ( $rem_start ) { $px_start++; $contrib -= $rem_start; }
            push @px, [ $px_start, $px_end, { $color=>$contrib } ];
        }
        else { error("really?!\n"); }
    }
    return \@px;
}

sub chr_dims
{
    my ($opts) = @_;

    # http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
    my $ktype = q[
        chr1	121500000	125000000	p11.1	acen
        chr1	125000000	128900000	q11	acen
        chr1	243700000	249250621	q44	gneg
        chr10	38000000	40200000	p11.1	acen
        chr10	40200000	42300000	q11.1	acen
        chr10	130600000	135534747	q26.3	gneg
        chr11	51600000	53700000	p11.11	acen
        chr11	53700000	55700000	q11	acen
        chr11	130800000	135006516	q25	gneg
        chr12	33300000	35800000	p11.1	acen
        chr12	35800000	38200000	q11	acen
        chr12	129300000	133851895	q24.33	gneg
        chr13	16300000	17900000	p11.1	acen
        chr13	17900000	19500000	q11	acen
        chr13	110300000	115169878	q34	gneg
        chr14	16100000	17600000	p11.1	acen
        chr14	17600000	19100000	q11.1	acen
        chr14	104000000	107349540	q32.33	gneg
        chr15	15800000	19000000	p11.1	acen
        chr15	19000000	20700000	q11.1	acen
        chr15	98500000	102531392	q26.3	gneg
        chr16	34600000	36600000	p11.1	acen
        chr16	36600000	38600000	q11.1	acen
        chr16	88700000	90354753	q24.3	gneg
        chr17	22200000	24000000	p11.1	acen
        chr17	24000000	25800000	q11.1	acen
        chr17	75300000	81195210	q25.3	gneg
        chr18	15400000	17200000	p11.1	acen
        chr18	17200000	19000000	q11.1	acen
        chr18	73100000	78077248	q23	gneg
        chr19	24400000	26500000	p11	acen
        chr19	26500000	28600000	q11	acen
        chr19	56300000	59128983	q13.43	gpos25
        chr2	90500000	93300000	p11.1	acen
        chr2	93300000	96800000	q11.1	acen
        chr2	237300000	243199373	q37.3	gneg
        chr20	25600000	27500000	p11.1	acen
        chr20	27500000	29400000	q11.1	acen
        chr20	58400000	63025520	q13.33	gneg
        chr21	10900000	13200000	p11.1	acen
        chr21	13200000	14300000	q11.1	acen
        chr21	42600000	48129895	q22.3	gneg
        chr22	12200000	14700000	p11.1	acen
        chr22	14700000	17900000	q11.1	acen
        chr22	49400000	51304566	q13.33	gneg
        chr3	87900000	91000000	p11.1	acen
        chr3	91000000	93900000	q11.1	acen
        chr3	192300000	198022430	q29	gneg
        chr4	48200000	50400000	p11	acen
        chr4	50400000	52700000	q11	acen
        chr4	187100000	191154276	q35.2	gpos25
        chr5	46100000	48400000	p11	acen
        chr5	48400000	50700000	q11.1	acen
        chr5	176600000	180915260	q35.3	gneg
        chr6	58700000	61000000	p11.1	acen
        chr6	61000000	63300000	q11.1	acen
        chr6	164500000	171115067	q27	gneg
        chr7	58000000	59900000	p11.1	acen
        chr7	59900000	61700000	q11.1	acen
        chr7	155100000	159138663	q36.3	gneg
        chr8	43100000	45600000	p11.1	acen
        chr8	45600000	48100000	q11.1	acen
        chr8	139900000	146364022	q24.3	gneg
        chr9	47300000	49000000	p11.1	acen
        chr9	49000000	50700000	q11	acen
        chr9	137400000	141213431	q34.3	gneg
        chrX	58100000	60600000	p11.1	acen
        chrX	60600000	63000000	q11.1	acen
        chrX	147100000	155270560	q28	gneg
    ];
    my @lines = split(/\n/,$ktype);
    my %chrs  = ();
    my $max_len = 0;
    for my $line (@lines)
    {
        if ( $line =~ /^\s*$/ ) { next; }
        $line =~ s/^\s*//;
        $line =~ s/\s*$//;
        my @items = split(/\s+/,$line);
        my $chr   = $items[0];
        my $start = $items[1];
        my $end   = $items[2];
        $chr =~ s/^chr//;
        if ( $items[-1] eq 'acen' )
        {
            push @{$chrs{$chr}{acen}},$start;
            push @{$chrs{$chr}{acen}},$end;
        }
        if ( !defined $chrs{$chr}{len} or $chrs{$chr}{len} < $end )
        { 
            $chrs{$chr}{len} = $end; 
        }
        if ( $max_len < $end ) { $max_len = $end; }
    }
    for my $chr (keys %chrs)
    {
        splice(@{$chrs{$chr}{acen}},1,1);
    }
    $$opts{chrs} = \%chrs;
    $$opts{max_chr_len} = $max_len;
}




