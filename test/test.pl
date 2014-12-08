#!/usr/bin/env perl
#
#   Copyright (C) 2012-2014 Genome Research Ltd.
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
#use IPC::Open2;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;

my $opts = parse_params();

test_usage($opts,cmd=>'bcftools');
test_tabix($opts,in=>'merge.a',reg=>'2:3199812-3199812',out=>'tabix.2.3199812.out');
test_tabix($opts,in=>'merge.a',reg=>'1:3000151-3000151',out=>'tabix.1.3000151.out');
test_tabix($opts,in=>'large_chrom_tbi_limit',reg=>'chr11:1-536870912',out=>'large_chrom_tbi_limit.20.1.536870912.out'); # 536870912 (1<<29) is the current limit for tbi. cannot retrieve regions larger than that
test_index($opts,in=>'large_chrom_csi_limit',reg=>'chr20:1-2147483647',out=>'large_chrom_csi_limit.20.1.2147483647.out'); # 2147483647 (1<<31-1) is the current chrom limit for csi. bcf conversion and indexing fail above this
test_index($opts,in=>'large_chrom_csi_limit',reg=>'chr20',out=>'large_chrom.20.1.2147483647.out'); # this fails until bug resolved
test_vcf_idxstats($opts,in=>'idx',args=>'-s',out=>'idx.out');
test_vcf_idxstats($opts,in=>'idx',args=>'-n',out=>'idx_count.out');
test_vcf_idxstats($opts,in=>'empty',args=>'-s',out=>'empty.idx.out');
test_vcf_idxstats($opts,in=>'empty',args=>'-n',out=>'empty.idx_count.out');
test_vcf_check($opts,in=>'check',out=>'check.chk');
test_vcf_stats($opts,in=>['stats.a','stats.b'],out=>'stats.chk',args=>'-s -');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.out',args=>'-n =2');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.flt.out',args=>'-n =2 -i"STRLEN(REF)==2"');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.both.out',args=>'-n =2 -c both');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.any.out',args=>'-n =2 -c any');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.C.out',args=>'-C -c any');
test_vcf_isec2($opts,vcf_in=>['isec.a'],tab_in=>'isec',out=>'isec.tab.out',args=>'');
test_vcf_merge($opts,in=>['merge.a','merge.b','merge.c'],out=>'merge.abc.out',args=>'--force-samples');
test_vcf_merge($opts,in=>['merge.2.a','merge.2.b'],out=>'merge.2.none.out',args=>'--force-samples -m none');
test_vcf_merge($opts,in=>['merge.2.a','merge.2.b'],out=>'merge.2.both.out',args=>'--force-samples -m both');
test_vcf_merge($opts,in=>['merge.2.a','merge.2.b'],out=>'merge.2.all.out',args=>'--force-samples -m all');
test_vcf_merge($opts,in=>['merge.3.a','merge.3.b'],out=>'merge.3.out',args=>'--force-samples -i TR:sum,TA:sum,TG:sum');
test_vcf_merge($opts,in=>['merge.4.a','merge.4.b'],out=>'merge.4.out',args=>'--force-samples -m id');
test_vcf_query($opts,in=>'query',out=>'query.out',args=>q[-f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%DP4\\t%AN[\\t%GT\\t%TGT]\\n']);
test_vcf_query($opts,in=>'view.filter',out=>'query.2.out',args=>q[-f'%XRI\\n' -i'XRI[*]>1111']);
test_vcf_query($opts,in=>'view.filter',out=>'query.3.out',args=>q[-f'%XRF\\n' -i'XRF[*]=2e6']);
test_vcf_query($opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[5]="PQR"']);
test_vcf_query($opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[*]="GHI"']);
test_vcf_query($opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[2]~"H"']);
test_vcf_query($opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[3]!~"H"']);
test_vcf_query($opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[*]~"H"']);
test_vcf_query($opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[*]~"HI,JK"']);
test_vcf_query($opts,in=>'query',out=>'query.5.out',args=>q[-f'%POS %REF %ALT\\n' -i'REF~"C" && ALT~"CT"']);
test_vcf_query($opts,in=>'query',out=>'query.6.out',args=>q[-f'%POS %REF %ALT\\n' -i'N_ALT=2']);
test_vcf_query($opts,in=>'query',out=>'query.7.out',args=>q[-f'%POS %AN\\n' -i'AN!=2*N_SAMPLES']);
test_vcf_query($opts,in=>'query',out=>'query.8.out',args=>q[-f'%POS[ %GL]\\n' -i'min(abs(GL[0]))=10']);
test_vcf_query($opts,in=>'view.filter',out=>'query.9.out',args=>q[-f'%POS %CIGAR\\n' -i'strlen(CIGAR[*])=4']);
test_vcf_query($opts,in=>'query',out=>'query.10.out',args=>q[-f'%POS[ %GT]\\n' -i'AC[0]=3']);
test_vcf_query($opts,in=>'query',out=>'query.10.out',args=>q[-f'%POS[ %GT]\\n' -i'AF[0]=3/4']);
test_vcf_query($opts,in=>'query',out=>'query.11.out',args=>q[-f'%POS[ %GT]\\n' -i'MAC[0]=1']);
test_vcf_query($opts,in=>'query',out=>'query.11.out',args=>q[-f'%POS[ %GT]\\n' -i'MAF[0]=1/4']);
test_vcf_query($opts,in=>'view.vectors',out=>'query.12.out',args=>q[-f'I8=%I8 I16=%I16 I32=%I32 IF=%IF IA8=%IA8 IA16=%IA16 IA32=%IA32 IAF=%IAF IA8=%IA8{1} IA16=%IA16{1} IA32=%IA32{1} IAF=%IAF{1} [ %F8:%F16:%F32:%FF]\\n']);
test_vcf_query($opts,in=>'query.filter',out=>'query.13.out',args=>q[-f'%POS[ %GT]\\n' -i'GT ="1"']);
test_vcf_query($opts,in=>'query.filter',out=>'query.14.out',args=>q[-f'%POS[ %GT]\\n' -i'GT!="1"']);
test_vcf_query($opts,in=>'query.filter',out=>'query.15.out',args=>q[-f'%POS[ %GT]\\n' -e'GT ="1"']);
test_vcf_query($opts,in=>'query.filter',out=>'query.16.out',args=>q[-f'%POS[ %GT]\\n' -e'GT!="1"']);
test_vcf_query($opts,in=>'query.2',out=>'query.17.out',args=>q[-f'%XX_A %XX.A %XX.A0 %xx.a0\\n']);
test_vcf_norm($opts,in=>'norm',out=>'norm.out',fai=>'norm');
test_vcf_norm($opts,in=>'norm.split',out=>'norm.split.out',args=>'-m-');
test_vcf_norm($opts,in=>'norm.merge',out=>'norm.merge.out',args=>'-m+');
test_vcf_norm($opts,in=>'norm.merge',out=>'norm.merge.strict.out',args=>'-m+ -s');
test_vcf_view($opts,in=>'view',out=>'view.1.out',args=>'-aUc1 -C1 -s NA00002 -v snps',reg=>'');
test_vcf_view($opts,in=>'view',out=>'view.2.out',args=>'-f PASS -Xks NA00003',reg=>'-r20,Y');
test_vcf_view($opts,in=>'view',out=>'view.3.out',args=>'-xs NA00003',reg=>'');
test_vcf_view($opts,in=>'view',out=>'view.4.out',args=>q[-i 'QUAL==999 && (FS<20 || FS>=41.02) && ICF>-0.1 && HWE*2>1.2'],reg=>'');
test_vcf_view($opts,in=>'view',out=>'view.5.out',args=>q[-p],reg=>'');
test_vcf_view($opts,in=>'view',out=>'view.6.out',args=>q[-P],reg=>'');
test_vcf_view($opts,in=>'view',out=>'view.7.out',args=>q[-hm2 -M2 -q0.3 -Q0.7],reg=>'');
test_vcf_view($opts,in=>'view',out=>'view.8.out',args=>q[-Hu],reg=>'');
test_vcf_view($opts,in=>'view',out=>'view.9.out',args=>q[-GVsnps],reg=>'');
test_vcf_view($opts,in=>'view',out=>'view.10.out',args=>q[-ne 'INDEL=1 || PV4[0]<0.006'],reg=>'');
test_vcf_view($opts,in=>'view',out=>'view.exclude.out',args=>'-s ^NA00003',reg=>'');
test_vcf_view($opts,in=>'view.omitgenotypes',out=>'view.omitgenotypes.out',args=>'',reg=>'');
test_vcf_view($opts,in=>'view.vectors',out=>'view.vectors.A.out',args=>'-asA',reg=>'');
test_vcf_view($opts,in=>'view.vectors',out=>'view.vectors.B.out',args=>'-asB',reg=>'');
test_vcf_view($opts,in=>'view.filter',out=>'view.filter.1.out',args=>q[-H -i'FMT/FGS[0]="AAAAAA"'],reg=>'');    # test expressions
test_vcf_view($opts,in=>'view.filter',out=>'view.filter.2.out',args=>q[-H -i'FMT/FGS[2]="C"'],reg=>'');
test_vcf_view($opts,in=>'view.filter',out=>'view.filter.3.out',args=>q[-H -i'FMT/FGS[4]="EE"'],reg=>'');
test_vcf_view($opts,in=>'view.filter',out=>'view.filter.4.out',args=>q[-H -i'FMT/FRS[1]="BB"'],reg=>'');
test_vcf_view($opts,in=>'view.filter',out=>'view.filter.5.out',args=>q[-H -i'TXT0="text"'],reg=>'');
test_vcf_view($opts,in=>'view.chrs',out=>'view.chrs.out',args=>'',reg=>'',tgts=>'view.chrs.tab');
test_vcf_filter($opts,in=>'view.filter',out=>'view.filter.6.out',args=>q[-S. -e'TXT0="text"'],reg=>'');
test_vcf_filter($opts,in=>'view.filter',out=>'view.filter.7.out',args=>q[-S. -e'FMT/FRS[1]="BB"'],reg=>'');
test_vcf_filter($opts,in=>'view.filter',out=>'view.filter.8.out',args=>q[-S. -e'FMT/FGS[0]="AAAAAA"'],reg=>'');
test_vcf_filter($opts,in=>'view.filter',out=>'view.filter.9.out',args=>q[-S. -e'FMT/FGS[1]="BBB"'],reg=>'');
test_vcf_filter($opts,in=>'view.filter',out=>'view.filter.10.out',args=>q[-S. -e'FMT/FGS[4]="EE"'],reg=>'');
test_vcf_filter($opts,in=>'view.filter',out=>'view.filter.11.out',args=>q[-S. -e'FMT/STR="XX"'],reg=>'');
test_vcf_view($opts,in=>'view.minmaxac',out=>'view.minmaxac.1.out',args=>q[-H -C5:nonmajor],reg=>'');
test_vcf_view($opts,in=>'view.minmaxac',out=>'view.minmaxac.2.out',args=>q[-H -c6:nonmajor],reg=>'');
test_vcf_view($opts,in=>'view.minmaxac',out=>'view.minmaxac.1.out',args=>q[-H -q0.3:major],reg=>'');
test_vcf_call($opts,in=>'mpileup',out=>'mpileup.1.out',args=>'-mv');
test_vcf_call($opts,in=>'mpileup',out=>'mpileup.2.out',args=>'-mvg0');
test_vcf_call_cAls($opts,in=>'mpileup',out=>'mpileup.cAls.out',tab=>'mpileup');
test_vcf_filter($opts,in=>'filter.1',out=>'filter.1.out',args=>'-mx -g2 -G2');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.2.out',args=>q[-e'QUAL==59.2 || (INDEL=0 & (FMT/GQ=25 | FMT/DP=10))' -sModified -S.]);
test_vcf_filter($opts,in=>'filter.3',out=>'filter.3.out',args=>q[-e'DP=19'],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.4.out',args=>q[-e'DP=19' -s XX],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.5.out',args=>q[-e'DP=19' -s XX -m+],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.6.out',args=>q[-e'DP=19' -s XX -mx],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.7.out',args=>q[-e'DP=19' -s XX -m+x],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.3.out',args=>q[-e'FMT/GT="0/2"'],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.4.out',args=>q[-e'FMT/GT="0/2"' -s XX],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.5.out',args=>q[-e'FMT/GT="0/2"' -s XX -m+],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.6.out',args=>q[-e'FMT/GT="0/2"' -s XX -mx],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.7.out',args=>q[-e'FMT/GT="0/2"' -s XX -m+x],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.8.out',args=>q[-i'FMT/GT="0/0" && AC[*]=2'],fmt=>'%POS\\t%AC[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.8.out',args=>q[-i'AC[*]=2 && FMT/GT="0/0"'],fmt=>'%POS\\t%AC[\\t%GT]\\n');
test_vcf_regions($opts,in=>'regions');
test_vcf_annotate($opts,in=>'annotate',tab=>'annotate',out=>'annotate.out',args=>'-c CHROM,POS,REF,ALT,ID,QUAL,INFO/T_INT,INFO/T_FLOAT,INDEL');
test_vcf_annotate($opts,in=>'annotate',tab=>'annotate2',out=>'annotate2.out',args=>'-c CHROM,FROM,TO,T_STR');
test_vcf_annotate($opts,in=>'annotate',vcf=>'annots',out=>'annotate3.out',args=>'-c STR,ID,QUAL,FILTER');
test_vcf_annotate($opts,in=>'annotate2',vcf=>'annots2',out=>'annotate4.out',args=>'-c ID,QUAL,FILTER,INFO,FMT');
test_vcf_annotate($opts,in=>'annotate2',vcf=>'annots2',out=>'annotate5.out',args=>'-c ID,QUAL,+FILTER,+INFO,FMT/GT -s A');
test_vcf_annotate($opts,in=>'annotate3',out=>'annotate6.out',args=>'-x ID,QUAL,^FILTER/fltA,FILTER/fltB,^INFO/AA,INFO/BB,^FMT/GT,FMT/PL');
test_vcf_annotate($opts,in=>'annotate3',out=>'annotate7.out',args=>'-x FORMAT');
test_vcf_plugin($opts,in=>'plugin1',out=>'missing2ref.out',cmd=>'+missing2ref');
test_vcf_plugin($opts,in=>'plugin1',out=>'fill-AN-AC.out',cmd=>'+fill-AN-AC');
test_vcf_plugin($opts,in=>'plugin1',out=>'dosage.out',cmd=>'+dosage');
test_vcf_plugin($opts,in=>'fixploidy',out=>'fixploidy.out',cmd=>'+fixploidy',args=>'-- -s {PATH}/fixploidy.samples -p {PATH}/fixploidy.ploidy');
test_vcf_plugin($opts,in=>'vcf2sex',out=>'vcf2sex.out',cmd=>'+vcf2sex',args=>'-- -n 5');
test_vcf_plugin($opts,in=>'vcf2sex',out=>'vcf2sex.out',cmd=>'+vcf2sex',args=>'-- -gn 5');
test_vcf_concat($opts,in=>['concat.1.a','concat.1.b'],out=>'concat.1.vcf.out',do_bcf=>0,args=>'');
test_vcf_concat($opts,in=>['concat.1.a','concat.1.b'],out=>'concat.1.bcf.out',do_bcf=>1,args=>'');
test_vcf_concat($opts,in=>['concat.2.a','concat.2.b'],out=>'concat.2.vcf.out',do_bcf=>0,args=>'-a');
test_vcf_concat($opts,in=>['concat.2.a','concat.2.b'],out=>'concat.2.bcf.out',do_bcf=>1,args=>'-a');
test_vcf_concat($opts,in=>['concat.2.a','concat.2.b'],out=>'concat.4.vcf.out',do_bcf=>0,args=>'-aD');
test_vcf_concat($opts,in=>['concat.2.a','concat.2.b'],out=>'concat.4.bcf.out',do_bcf=>1,args=>'-aD');
test_vcf_concat($opts,in=>['concat.3.a','concat.3.b','concat.3.0','concat.3.c','concat.3.d','concat.3.e','concat.3.f'],out=>'concat.3.vcf.out',do_bcf=>0,args=>'-l');
test_vcf_concat($opts,in=>['concat.3.a','concat.3.b','concat.3.0','concat.3.c','concat.3.d','concat.3.e','concat.3.f'],out=>'concat.3.bcf.out',do_bcf=>1,args=>'-l');
test_vcf_reheader($opts,in=>'reheader',out=>'reheader.1.out',header=>'reheader.hdr');
test_vcf_reheader($opts,in=>'reheader',out=>'reheader.2.out',samples=>'reheader.samples');
test_vcf_reheader($opts,in=>'reheader',out=>'reheader.2.out',samples=>'reheader.samples2');
test_rename_chrs($opts,in=>'annotate');
test_vcf_convert($opts,in=>'convert',out=>'convert.gs.gt.gen',args=>'-g -,.');
test_vcf_convert($opts,in=>'convert',out=>'convert.gs.gt.samples',args=>'-g .,-');
test_vcf_convert($opts,in=>'convert',out=>'convert.gs.pl.gen',args=>'-g -,. --tag PL');
test_vcf_convert($opts,in=>'convert',out=>'convert.gs.pl.samples',args=>'-g .,- --tag PL');
test_vcf_convert($opts,in=>'convert',out=>'convert.hls.haps',args=>'-h -,.,.');
test_vcf_convert($opts,in=>'convert',out=>'convert.hls.legend',args=>'-h .,-,.');
test_vcf_convert($opts,in=>'convert',out=>'convert.hls.samples',args=>'-h .,.,-');
test_vcf_convert_tsv2vcf($opts,in=>'convert.23andme',out=>'convert.23andme.vcf',args=>'-c ID,CHROM,POS,AA -s SAMPLE1',fai=>'23andme');
test_vcf_consensus($opts,in=>'consensus',out=>'consensus.1.out',fa=>'consensus.fa',mask=>'consensus.tab',args=>'');
test_vcf_consensus_chain($opts,in=>'consensus',out=>'consensus.1.chain',chain=>'consensus.1.chain',fa=>'consensus.fa',mask=>'consensus.tab',args=>'');
test_vcf_consensus($opts,in=>'consensus',out=>'consensus.2.out',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-H 1');
test_vcf_consensus_chain($opts,in=>'consensus',out=>'consensus.2.chain',chain=>'consensus.2.chain',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-H 1');
test_vcf_consensus($opts,in=>'consensus',out=>'consensus.3.out',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-i');
test_vcf_consensus_chain($opts,in=>'consensus',out=>'consensus.3.chain',chain=>'consensus.3.chain',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-i');
test_vcf_consensus($opts,in=>'consensus',out=>'consensus.4.out',fa=>'consensus.fa',args=>'-H 1');
test_vcf_consensus_chain($opts,in=>'consensus',out=>'consensus.4.chain',chain=>'consensus.4.chain',fa=>'consensus.fa',args=>'-H 1');

print "\nNumber of tests:\n";
printf "    total   .. %d\n", $$opts{nok}+$$opts{nfailed};
printf "    passed  .. %d\n", $$opts{nok};
printf "    failed  .. %d\n", $$opts{nfailed};
print "\n";

exit ($$opts{nfailed} != 0);

#--------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print
        "About: htslib consistency test script\n",
        "Usage: test.pl [OPTIONS]\n",
        "Options:\n",
        "   -p, --plugins                   Test also plugins, requires libhts.so.\n",
        "   -r, --redo-outputs              Recreate expected output files.\n",
        "   -t, --temp-dir <path>           When given, temporary files will not be removed.\n",
        "   -h, -?, --help                  This help message.\n",
        "\n";
    exit -1;
}
sub parse_params
{
    my $opts = { bgzip=>"bgzip", keep_files=>0, nok=>0, nfailed=>0, tabix=>"tabix", plugins=>0 };
    my $help;
    Getopt::Long::Configure('bundling');
    my $ret = GetOptions (
            'e|exec=s' => sub { my ($tool, $path) = split /=/, $_[1]; $$opts{$tool} = $path if $path },
            't|temp-dir:s' => \$$opts{keep_files},
            'p|plugins' => \$$opts{test_plugins},
            'r|redo-outputs' => \$$opts{redo_outputs},
            'h|?|help' => \$help
            );
    if ( !$ret or $help ) { error(); }
    $$opts{tmp} = $$opts{keep_files} ? $$opts{keep_files} : tempdir(CLEANUP=>1);
    if ( $$opts{keep_files} ) { cmd("mkdir -p $$opts{keep_files}"); }
    $$opts{path} = $FindBin::RealBin;
    $$opts{bin}  = $FindBin::RealBin;
    $$opts{bin}  =~ s{/test/?$}{};
    return $opts;
}
sub _cmd
{
    my ($cmd) = @_;
    my $kid_io;
    my @out;
    my $pid = open($kid_io, "-|");
    if ( !defined $pid ) { error("Cannot fork: $!"); }
    if ($pid)
    {
        # parent
        @out = <$kid_io>;
        close($kid_io);
    }
    else
    {
        # child
        exec('/bin/bash', '-o','pipefail','-c', $cmd) or error("Cannot execute the command [/bin/sh -o pipefail -c $cmd]: $!");
    }
    return ($? >> 8, join('',@out));
}
sub cmd
{
    my ($cmd) = @_;
    my ($ret,$out) = _cmd($cmd);
    if ( $ret ) { error("The command failed: $cmd\n", $out); }
    return $out;
}
sub test_cmd
{
    my ($opts,%args) = @_;
    if ( !exists($args{out}) )
    {
        if ( !exists($args{in}) ) { error("FIXME: expected out or in key\n"); }
        $args{out} = "$args{in}.out";
    }
    my ($package, $filename, $line, $test)=caller(1);
    $test =~ s/^.+:://;

    print "$test:\n";
    print "\t$args{cmd}\n";

    my ($ret,$out) = _cmd("$args{cmd}");
    if ( $ret ) { failed($opts,$test,"Non-zero status $ret"); return; }
    if ( $$opts{redo_outputs} && -e "$$opts{path}/$args{out}" )
    {
        rename("$$opts{path}/$args{out}","$$opts{path}/$args{out}.old");
        open(my $fh,'>',"$$opts{path}/$args{out}") or error("$$opts{path}/$args{out}: $!");
        print $fh $out;
        close($fh);
        my ($ret,$out) = _cmd("diff -q $$opts{path}/$args{out} $$opts{path}/$args{out}.old");
        if ( !$ret && $out eq '' ) { unlink("$$opts{path}/$args{out}.old"); }
        else
        {
            print "\tthe expected output changed, saving:\n";
            print "\t  old .. $$opts{path}/$args{out}.old\n";
            print "\t  new .. $$opts{path}/$args{out}\n";
        }
    }
    my $exp = '';
    if ( open(my $fh,'<',"$$opts{path}/$args{out}") )
    {
        my @exp = <$fh>;
        $exp = join('',@exp);
        close($fh);
    }
    elsif ( !$$opts{redo_outputs} ) { failed($opts,$test,"$$opts{path}/$args{out}: $!"); return; }

    if ( $exp ne $out )
    {
        open(my $fh,'>',"$$opts{path}/$args{out}.new") or error("$$opts{path}/$args{out}.new");
        print $fh $out;
        close($fh);
        if ( !-e "$$opts{path}/$args{out}" )
        {
            rename("$$opts{path}/$args{out}.new","$$opts{path}/$args{out}") or error("rename $$opts{path}/$args{out}.new $$opts{path}/$args{out}: $!");
            print "\tthe file with expected output does not exist, creating new one:\n";
            print "\t\t$$opts{path}/$args{out}\n";
        }
        else
        {
            failed($opts,$test,"The outputs differ:\n\t\t$$opts{path}/$args{out}\n\t\t$$opts{path}/$args{out}.new");
        }
        return;
    }
    passed($opts,$test);
}
sub failed
{
    my ($opts,$test,$reason) = @_;
    $$opts{nfailed}++;
    if ( defined $reason ) { print "\n\t$reason"; }
    print "\n.. failed ...\n\n";
}
sub passed
{
    my ($opts,$test) = @_;
    $$opts{nok}++;
    print ".. ok\n\n";
}
sub is_file_newer
{
    my ($afile,$bfile) = @_;
    my (@astat) = stat($afile) or return 0;
    my (@bstat) = stat($bfile) or return 0;
    if ( $astat[9]>$bstat[9] ) { return 1 }
    return 0;
}
sub bgzip_tabix
{
    my ($opts,%args) = @_;
    my $file = "$args{file}.$args{suffix}";
    if ( $$opts{redo_outputs} or !-e "$$opts{tmp}/$file.gz" or is_file_newer("$$opts{path}/$file","$$opts{tmp}/$file.gz") )
    {
        cmd("cat $$opts{path}/$file | $$opts{bgzip} -c > $$opts{tmp}/$file.gz");
    }
    if ( $$opts{redo_outputs} or !-e "$$opts{tmp}/$file.gz.tbi" or is_file_newer("$$opts{tmp}/$file.gz","$$opts{tmp}/$file.gz.tbi") )
    {
        cmd("$$opts{tabix} -f $args{args} $$opts{tmp}/$file.gz");
    }
}
sub bgzip_tabix_vcf
{
    my ($opts,$file) = @_;
    bgzip_tabix($opts,file=>$file,suffix=>'vcf',args=>'-p vcf');
}


# The tests --------------------------

sub test_tabix
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{tabix} $$opts{tmp}/$args{in}.vcf.gz $args{reg}");

    cmd("$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz > $$opts{tmp}/$args{in}.bcf");
    cmd("$$opts{bin}/bcftools index -f $$opts{tmp}/$args{in}.bcf");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -H $$opts{tmp}/$args{in}.bcf $args{reg}");
}
sub test_index
{
    my ($opts,%args) = @_;
    cmd("$$opts{bin}/bcftools view -Oz $$opts{path}/$args{in}.vcf > $$opts{tmp}/$args{in}.vcf.gz");
    cmd("$$opts{bin}/bcftools index -f $$opts{tmp}/$args{in}.vcf.gz");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -H $$opts{tmp}/$args{in}.vcf.gz $args{reg}");

    cmd("$$opts{bin}/bcftools view -Ob $$opts{path}/$args{in}.vcf > $$opts{tmp}/$args{in}.bcf");
    cmd("$$opts{bin}/bcftools index -f $$opts{tmp}/$args{in}.bcf");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -H $$opts{tmp}/$args{in}.bcf $args{reg}");
}

sub test_vcf_idxstats
{
    my ($opts,%args) = @_;
    cmd("$$opts{bin}/bcftools view -Oz $$opts{path}/$args{in}.vcf > $$opts{tmp}/$args{in}.vcf.gz");
    cmd("$$opts{bin}/bcftools index --tbi -f $$opts{tmp}/$args{in}.vcf.gz");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools index $args{args} $$opts{tmp}/$args{in}.vcf.gz");
    cmd("$$opts{bin}/bcftools index --csi -f $$opts{tmp}/$args{in}.vcf.gz");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools index $args{args} $$opts{tmp}/$args{in}.vcf.gz");

    cmd("$$opts{bin}/bcftools view -Ob $$opts{path}/$args{in}.vcf > $$opts{tmp}/$args{in}.bcf");
    cmd("$$opts{bin}/bcftools index -f $$opts{tmp}/$args{in}.bcf");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools index $args{args} $$opts{tmp}/$args{in}.bcf");
}

sub test_vcf_check
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools stats -s - $$opts{tmp}/$args{in}.vcf.gz | grep -v '^# The command' | grep -v '^# This' | grep -v '^ID\t'");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools stats -s - | grep -v '^# The command' | grep -v '^# This' | grep -v '^ID\t'");
}
sub test_vcf_stats
{
    my ($opts,%args) = @_;
    my $files = '';
    for my $file (@{$args{in}})
    {
        bgzip_tabix_vcf($opts,$file);
        $files .= " $$opts{tmp}/$file.vcf.gz";
    }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools stats $args{args} $files | grep -v '^#' | grep -v '^ID\t'");
}
sub test_vcf_merge
{
    my ($opts,%args) = @_;
    my @files;
    for my $file (@{$args{in}})
    {
        bgzip_tabix_vcf($opts,$file);
        push @files, "$$opts{tmp}/$file.vcf.gz";
    }
    my $args  = exists($args{args}) ? $args{args} : '';
    my $files = join(' ',@files);
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools merge $args $files | grep -v ^##bcftools_");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools merge -Ob $args $files | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
}
sub test_vcf_isec
{
    my ($opts,%args) = @_;
    my @files;
    for my $file (@{$args{in}})
    {
        bgzip_tabix_vcf($opts,$file);
        push @files, "$$opts{tmp}/$file.vcf.gz";
    }
    my $files = join(' ',@files);
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools isec $args{args} $files");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools isec -Ob $args{args} $files");
}
sub test_vcf_isec2
{
    my ($opts,%args) = @_;
    my @files;
    for my $file (@{$args{vcf_in}})
    {
        bgzip_tabix_vcf($opts,$file);
        push @files, "$$opts{tmp}/$file.vcf.gz";
    }
    my $files = join(' ',@files);
    bgzip_tabix($opts,file=>$args{tab_in},suffix=>'tab',args=>'-s 1 -b 2 -e 3');
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools isec $args{args} -T $$opts{tmp}/$args{tab_in}.tab.gz $files 2>/dev/null | grep -v ^##bcftools_");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools isec -Ob $args{args} -T $$opts{tmp}/$args{tab_in}.tab.gz $files 2>/dev/null | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
}
sub test_vcf_query
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools query $args{args} $$opts{tmp}/$args{in}.vcf.gz");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools query $args{args}");
}
sub test_vcf_convert
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools convert $args{args} $$opts{tmp}/$args{in}.vcf.gz");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools convert $args{args}");
}
sub test_vcf_convert_tsv2vcf
{
    my ($opts,%args) = @_;
    my $params = '';
    if ( exists($args{args}) ) { $params .= " $args{args}"; }
    if ( exists($args{fai} ) ) { $params .= " -f $$opts{path}/$args{fai}.fa"; }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools convert $params --tsv2vcf $$opts{path}/$args{in} | grep -v ^##bcftools_");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools convert -Ou $params --tsv2vcf $$opts{path}/$args{in} | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
}
sub test_vcf_norm
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    my $params = '';
    if ( exists($args{args}) ) { $params .= " $args{args}"; }
    if ( exists($args{fai} ) ) { $params .= " -f $$opts{path}/$args{fai}.fa"; }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools norm $params $$opts{tmp}/$args{in}.vcf.gz | grep -v ^##bcftools_");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools norm -Ob $params $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
}
sub test_vcf_view
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});

    if ( !exists($args{args}) ) { $args{args} = ''; }
    if ( exists($args{tgts}) ) { $args{args} .= "-T $$opts{path}/$args{tgts}"; }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view $args{args} $$opts{tmp}/$args{in}.vcf.gz $args{reg} | grep -v ^##bcftools_");
    unless ($args{args} =~ /-H/) {
        test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -Ob $args{args} $$opts{tmp}/$args{in}.vcf.gz $args{reg} | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
    }
}
sub test_vcf_call
{
    my ($opts,%args) = @_;
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools call $args{args} $$opts{path}/$args{in}.vcf | grep -v ^##bcftools_");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools call -Ob $args{args} $$opts{path}/$args{in}.vcf | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
}
sub test_vcf_call_cAls
{
    my ($opts,%args) = @_;
    bgzip_tabix($opts,file=>$args{tab},suffix=>'tab',args=>'-s1 -b2 -e2');
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools call -mA -C alleles -T $$opts{tmp}/$args{tab}.tab.gz $$opts{path}/$args{in}.vcf | grep -v ^##bcftools_");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools call -Ob -mA -C alleles -T $$opts{tmp}/$args{tab}.tab.gz $$opts{path}/$args{in}.vcf | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
}
sub test_vcf_filter
{
    my ($opts,%args) = @_;
    my $pipe = 'grep -v ^##bcftools_';
    if ( exists($args{fmt}) )
    {
        $pipe = "$$opts{bin}/bcftools query -f '$args{fmt}'";
    }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools filter $args{args} $$opts{path}/$args{in}.vcf | $pipe");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools filter -Ob $args{args} $$opts{path}/$args{in}.vcf | $$opts{bin}/bcftools view | $pipe");
}
sub test_vcf_regions
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});

    # regions vs targets, holding tab in memory
    my $query = q[%CHROM %POS %REF,%ALT\n];
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -T $$opts{path}/$args{in}.tab $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools query -f'$query' -T $$opts{path}/$args{in}.tab],out=>'regions.out');
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -R $$opts{path}/$args{in}.tab $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');

    # regions vs targets, reading tabix-ed tab
    cmd(qq[cat $$opts{path}/$args{in}.tab | $$opts{bgzip} -c > $$opts{tmp}/$args{in}.tab.gz]);
    cmd(qq[$$opts{tabix} -f -s1 -b2 -e3 $$opts{tmp}/$args{in}.tab.gz]);
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -T $$opts{tmp}/$args{in}.tab.gz $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools query -f'$query' -T $$opts{tmp}/$args{in}.tab.gz],out=>'regions.out');
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -R $$opts{tmp}/$args{in}.tab.gz $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');

    # regions vs targets, holding bed in memory
    cmd(qq[cat $$opts{path}/$args{in}.tab | awk '{OFS="\\t"}{print \$1,\$2-1,\$3}' > $$opts{tmp}/$args{in}.bed]);
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -T $$opts{tmp}/$args{in}.bed $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools query -f'$query' -T $$opts{tmp}/$args{in}.bed],out=>'regions.out');
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -R $$opts{tmp}/$args{in}.bed $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');

    # regions vs targets, reading tabix-ed bed
    cmd(qq[cat $$opts{tmp}/$args{in}.bed | $$opts{bgzip} -c > $$opts{tmp}/$args{in}.bed.gz]);
    cmd(qq[$$opts{tabix} -f -p bed $$opts{tmp}/$args{in}.bed.gz]);
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -T $$opts{tmp}/$args{in}.bed.gz $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools query -f'$query' -T $$opts{tmp}/$args{in}.bed.gz],out=>'regions.out');
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -R $$opts{tmp}/$args{in}.bed.gz $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');
}
sub test_usage
{
    my ($opts,%args) = @_;

    my $test = "test_usage";
    print "$test:\n";
    print "\t$args{cmd}\n";

    my $tty_input;
    if (-t) {
        $args{redirection} = "";  # no redirection necessary
    }
    elsif (eval { require IO::Pty }) {
        $tty_input = new IO::Pty;
        # ensure stdin is a terminal, so that subcommands display their usage
        $args{redirection} = "<'" . $tty_input->ttyname . "'";
    }
    else {
        warn "$0: module IO::Pty not found; skipping usage tests\n";
        return;
    }

    my $command = $args{cmd};
    my $commandpath = $$opts{bin}."/".$command;
    my ($ret,$out) = _cmd("$commandpath $args{redirection} 2>&1");
    if ( $out =~ m/\/bin\/bash.*no.*such/i ) { failed($opts,$test,"could not run $commandpath: $out"); return; }

    my @sections = ($out =~ m/(^[A-Za-z]+.*?)(?:(?=^[A-Za-z]+:)|\z)/msg);

    my $have_usage = 0;
    my $have_version = 0;
    my $have_subcommands = 0;
    my $usage = "";
    my @subcommands = ();
    foreach my $section (@sections) {
        if ( $section =~ m/^usage/i ) {
            $have_usage = 1;
            $section =~ s/^[[:word:]]+[[:punct:]]?[[:space:]]*//;
            $usage = $section;
        } elsif ( $section =~ m/^version/i ) {
            $have_version = 1;
        } elsif ( $section =~ m/^command/i ) {
            $have_subcommands = 1;
            foreach my $line (split /\n/, $section) {
                push @subcommands, $1 if $line =~ /^\s{2,}(\w+)\s{2,}/;
            }
        }
    }

    if ( !$have_usage ) { failed($opts,$test,"did not have Usage:"); return; }
    if ( !$have_version ) { failed($opts,$test,"did not have Version:"); return; }
    if ( !$have_subcommands ) { failed($opts,$test,"did not have Commands:"); return; }

    if ( !($usage =~ m/$command/) ) { failed($opts,$test,"usage did not mention $command"); return; }

    if ( scalar(@subcommands) < 1 ) { failed($opts,$test,"could not parse subcommands"); return; }

    passed($opts,$test);

    # now test subcommand usage as well
    foreach my $subcommand (@subcommands) {
        test_usage_subcommand($opts,%args,subcmd=>$subcommand);
    }
}
sub test_usage_subcommand
{
    my ($opts,%args) = @_;

    my $test = "test_usage_subcommand";
    print "$test:\n";
    print "\t$args{cmd} $args{subcmd}\n";

    my $command = $args{cmd};
    my $subcommand = $args{subcmd};
    my $commandpath = $$opts{bin}."/".$command;
    my ($ret,$out) = _cmd("$commandpath $subcommand $args{redirection} 2>&1");
    if ( $out =~ m/\/bin\/bash.*no.*such/i ) { failed($opts,$test,"could not run $commandpath $subcommand: $out"); return; }

    my @sections = ($out =~ m/(^[A-Za-z]+.*?)(?:(?=^[A-Za-z]+:)|\z)/msg);

    my $have_usage = 0;
    my $usage = "";
    foreach my $section (@sections) {
        if ( $section =~ m/^usage/i ) {
            $have_usage = 1;
            $section =~ s/^[[:word:]]+[[:punct:]]?[[:space:]]*//;
            $usage = $section;
        }
    }

    if ( !$have_usage ) { failed($opts,$test,"did not have Usage:"); return; }

    if ( !($usage =~ m/$command[[:space:]]+$subcommand/) ) { failed($opts,$test,"usage did not mention $command $subcommand"); return; }

    passed($opts,$test);
}
sub test_vcf_annotate
{
    my ($opts,%args) = @_;
    my ($annot_fname,$in_fname,$hdr);
    if ( exists($args{tab}) )
    {
        bgzip_tabix($opts,file=>$args{tab},suffix=>'tab',args=>'-s1 -b2 -e2');
        $annot_fname = "-a $$opts{tmp}/$args{tab}.tab.gz";
        $in_fname = "$$opts{path}/$args{in}.vcf";
        $hdr = "-h $$opts{path}/$args{in}.hdr";
    }
    elsif ( exists($args{vcf}) )
    {
        bgzip_tabix_vcf($opts,"$args{in}");
        bgzip_tabix_vcf($opts,$args{vcf});
        $annot_fname = "-a $$opts{tmp}/$args{vcf}.vcf.gz";
        $in_fname = "$$opts{tmp}/$args{in}.vcf.gz";
        $hdr = '';
    }
    else
    {
        $in_fname = "$$opts{path}/$args{in}.vcf";
        $annot_fname = '';
        $hdr = '';
    }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools annotate $annot_fname $hdr $args{args} $in_fname | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools annotate -Ob $annot_fname $hdr $args{args} $in_fname | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
}
sub test_vcf_plugin
{
    my ($opts,%args) = @_;
    if ( !$$opts{test_plugins} ) { return; }
    $ENV{BCFTOOLS_PLUGINS} = "$$opts{bin}/plugins";
    if ( !exists($args{args}) ) { $args{args} = ''; }
    $args{args} =~ s/{PATH}/$$opts{path}/g;
    bgzip_tabix_vcf($opts,"$args{in}");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools $args{cmd} $$opts{tmp}/$args{in}.vcf.gz $args{args} | grep -v ^##bcftools_");
}
sub test_vcf_concat
{
    my ($opts,%args) = @_;
    my $files;
    for my $file (@{$args{in}})
    {
        if ( $args{do_bcf} )
        {
            cmd("$$opts{bin}/bcftools view -Ob $$opts{tmp}/$file.vcf.gz > $$opts{tmp}/$file.bcf");
            cmd("$$opts{bin}/bcftools index -f $$opts{tmp}/$file.bcf");
            $files .= " $$opts{tmp}/$file.bcf";
        }
        else
        {
            bgzip_tabix_vcf($opts,$file);
            $files .= " $$opts{tmp}/$file.vcf.gz";
        }
    }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools concat $args{args} $files | grep -v ^##bcftools_");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools concat -Ob $args{args} $files | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
}
sub test_vcf_reheader
{
    my ($opts,%args) = @_;
    cmd("$$opts{bin}/bcftools view -Ob $$opts{path}/$args{in}.vcf > $$opts{tmp}/$args{in}.bcf");
    cmd("$$opts{bin}/bcftools view -Oz $$opts{path}/$args{in}.vcf > $$opts{tmp}/$args{in}.vcf.gz");

    my $arg = exists($args{header}) ? "-h $$opts{path}/$args{header}" : "-s $$opts{path}/$args{samples}";
    for my $file ("$$opts{path}/$args{in}.vcf","$$opts{tmp}/$args{in}.bcf","$$opts{tmp}/$args{in}.vcf.gz")
    {
        # bcf header lines can come in different order
        my %bcf_args = ();
        if ( $file=~/\.bcf$/ && -e "$$opts{path}/$args{out}.bcf" ) { %bcf_args = ( out=>"$args{out}.bcf" ); }
        test_cmd($opts,%args,%bcf_args,cmd=>"$$opts{bin}/bcftools reheader $arg $file | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
        test_cmd($opts,%args,%bcf_args,cmd=>"cat $file | $$opts{bin}/bcftools reheader $arg | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
    }
}
sub test_rename_chrs
{
    my ($opts,%args) = @_;
    cmd("$$opts{bin}/bcftools view -Ob $$opts{path}/$args{in}.vcf > $$opts{tmp}/$args{in}.bcf");
    cmd("$$opts{bin}/bcftools query -f'chr%CHROM\\t%POS\\n' $$opts{path}/$args{in}.vcf > $$opts{path}/rename.out.tmp");
    cmd("$$opts{bin}/bcftools query -f'%CHROM\\tchr%CHROM\\n' $$opts{path}/$args{in}.vcf | uniq > $$opts{tmp}/rename.map");
    my $prevfailed = $$opts{nfailed};
    for my $file ("$$opts{tmp}/$args{in}.bcf","$$opts{path}/$args{in}.vcf")
    {
        test_cmd($opts,%args,out=>"rename.out.tmp",cmd=>"$$opts{bin}/bcftools annotate --rename-chrs $$opts{tmp}/rename.map -Ov $file | $$opts{bin}/bcftools query -f'%CHROM\\t%POS\\n'");
        test_cmd($opts,%args,out=>"rename.out.tmp",cmd=>"$$opts{bin}/bcftools annotate --rename-chrs $$opts{tmp}/rename.map -Ob $file | $$opts{bin}/bcftools query -f'%CHROM\\t%POS\\n'");
    }
    unlink "$$opts{path}/rename.out.tmp" if $$opts{nfailed} == $prevfailed;
}
sub test_vcf_consensus
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    my $mask = $args{mask} ? "-m $$opts{path}/$args{mask}" : '';
    my $chain = $args{chain} ? "-c $$opts{path}/$args{chain}" : '';
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools consensus $$opts{tmp}/$args{in}.vcf.gz -f $$opts{path}/$args{fa} $args{args} $mask $chain");
}
sub test_vcf_consensus_chain
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    my $mask = $args{mask} ? "-m $$opts{path}/$args{mask}" : '';
    my $chain = $args{chain} ? "-c $$opts{path}/$args{chain}.new" : '';
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools consensus $$opts{tmp}/$args{in}.vcf.gz -f $$opts{path}/$args{fa} $args{args} $mask $chain > /dev/null; cat $$opts{path}/$args{chain}.new");
}

