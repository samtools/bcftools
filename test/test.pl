#!/usr/bin/env perl
#
#   Copyright (C) 2012-2019 Genome Research Ltd.
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
use Cwd qw/ abs_path /;

my $opts = parse_params();
test_usage($opts,cmd=>'bcftools');
test_tabix($opts,in=>'merge.a',reg=>'2:3199812-3199812',out=>'tabix.2.3199812.out');
test_tabix($opts,in=>'merge.a',reg=>'1:3000151-3000151',out=>'tabix.1.3000151.out');
test_index($opts,in=>'large_chrom_csi_limit',reg=>'chr20:1-2147483647',out=>'large_chrom_csi_limit.20.1.2147483647.out'); # 2147483647 (1<<31-1) is the current chrom limit for csi. bcf conversion and indexing fail above this
test_index($opts,in=>'large_chrom_csi_limit',reg=>'chr20',out=>'large_chrom.20.1.2147483647.out'); # this fails until bug resolved
test_vcf_idxstats($opts,in=>'idx',args=>'-s',out=>'idx.out');
test_vcf_idxstats($opts,in=>'idx',args=>'-n',out=>'idx_count.out');
test_vcf_idxstats($opts,in=>'empty',args=>'-s',out=>'empty.idx.out');
test_vcf_idxstats($opts,in=>'empty',args=>'-n',out=>'empty.idx_count.out');
test_vcf_check($opts,in=>'check',out=>'check.chk');
test_vcf_check_merge($opts,in=>'check',out=>'check_merge.chk');
test_vcf_stats($opts,in=>['stats.a','stats.b'],out=>'stats.chk',args=>'-s -');
test_vcf_stats($opts,in=>['stats.a','stats.b'],out=>'stats.B.chk',args=>'-s B');
test_vcf_stats($opts,in=>['stats.counts'],out=>'stats.counts.chk',args=>'-s -');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.out',args=>'-n =2');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.flt.out',args=>'-n =2 -i"STRLEN(REF)==2"');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.both.out',args=>'-n =2 -c both');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.any.out',args=>'-n =2 -c any');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.C.out',args=>'-C -c any');
test_vcf_isec2($opts,vcf_in=>['isec.a'],tab_in=>'isec',out=>'isec.tab.out',args=>'');
test_vcf_isec($opts,in=>['isec-miss.1.1','isec-miss.1.2','isec-miss.1.3'],out=>'isec-miss.1.1.out',args=>'-n +1 -r 20:100,20:140,12:55,20:140,20:100');
test_vcf_isec($opts,in=>['isec-miss.1.1','isec-miss.1.2','isec-miss.1.3'],out=>'isec-miss.1.1.out',args=>'-R {PATH}/isec-miss.1.regs.txt -n +1');
test_vcf_isec($opts,in=>['isec-miss.2.1','isec-miss.2.2','isec-miss.2.3'],out=>'isec-miss.2.1.out',args=>'-n +1 -r 20:100,20:140,12:55,20:140,20:100');
test_vcf_isec($opts,in=>['isec-miss.2.1','isec-miss.2.2','isec-miss.2.3'],out=>'isec-miss.2.1.out',args=>'-R {PATH}/isec-miss.1.regs.txt -n +1');
test_vcf_merge($opts,in=>['merge.a','merge.b','merge.c'],out=>'merge.abc.out',args=>'--force-samples');
test_vcf_merge($opts,in=>['merge.a','merge.b','merge.c'],out=>'merge.abc.2.out',args=>'--force-samples -Fx');
test_vcf_merge($opts,in=>['merge.a','merge.b','merge.c'],out=>'merge.abc.3.out',args=>'--force-samples -0');
test_vcf_merge($opts,in=>['merge.2.a','merge.2.b'],out=>'merge.2.none.out',args=>'--force-samples -m none');
test_vcf_merge($opts,in=>['merge.2.a','merge.2.b'],out=>'merge.2.both.out',args=>'--force-samples -m both');
test_vcf_merge($opts,in=>['merge.2.a','merge.2.b'],out=>'merge.2.all.out',args=>'--force-samples -m all');
test_vcf_merge($opts,in=>['merge.3.a','merge.3.b'],out=>'merge.3.out',args=>'--force-samples -i TR:sum,TA:sum,TG:sum');
test_vcf_merge($opts,in=>['merge.4.a','merge.4.b'],out=>'merge.4.out',args=>'--force-samples -m id');
test_vcf_merge($opts,in=>['gvcf.merge.1','gvcf.merge.2','gvcf.merge.3'],out=>'gvcf.merge.1.out',args=>'--gvcf -');
test_vcf_merge($opts,in=>['merge.gvcf.2.a','merge.gvcf.2.b','merge.gvcf.2.c'],out=>'merge.gvcf.2.out',args=>'--gvcf -');
test_vcf_merge($opts,in=>['merge.gvcf.3.a','merge.gvcf.3.b'],out=>'merge.gvcf.3.out',args=>'--gvcf - -i SRC:join');
test_vcf_merge($opts,in=>['merge.gvcf.4.a','merge.gvcf.4.b'],out=>'merge.gvcf.4.out',args=>'--gvcf -');
test_vcf_merge($opts,in=>['merge.5.a','merge.5.b'],out=>'merge.5.out');
test_vcf_merge($opts,in=>['merge.6.a','merge.6.b'],out=>'merge.6.out');
test_vcf_merge($opts,in=>['merge.gvcf.7.a','merge.gvcf.7.b'],out=>'merge.gvcf.7.out',args=>'--gvcf -');
test_vcf_merge($opts,in=>['merge.gvcf.8.a','merge.gvcf.8.b'],out=>'merge.gvcf.8.out',args=>'--gvcf -');
test_vcf_merge($opts,in=>['merge.7.a','merge.7.b'],out=>'merge.9.out',args=>'--force-samples');
test_vcf_query($opts,in=>'query',out=>'query.out',args=>q[-f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%DP4\\t%AN[\\t%GT\\t%TGT]\\n']);
test_vcf_query($opts,in=>'query.variantkey',out=>'query.variantkey.hex.out',args=>q[-f '%RSX\\t%VKX\\n']);
test_vcf_query($opts,in=>'view.filter',out=>'query.2.out',args=>q[-f'%XRI\\n' -i'XRI[*]>1111']);
test_vcf_query($opts,in=>'view.filter',out=>'query.3.out',args=>q[-f'%XRF\\n' -i'XRF[*]=2e6']);
test_vcf_query($opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[5]="PQR"']);
test_vcf_query($opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[*]="GHI"']);
test_vcf_query($opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[2]~"H"']);
test_vcf_query($opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[3]!~"H" && XGS!="."']);
test_vcf_query($opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[*]~"H"']);
test_vcf_query($opts,in=>'query',out=>'query.5.out',args=>q[-f'%POS %REF %ALT\\n' -i'REF~"C" && ALT[*]~"CT"']);
test_vcf_query($opts,in=>'query',out=>'query.6.out',args=>q[-f'%POS %REF %ALT\\n' -i'N_ALT=2']);
test_vcf_query($opts,in=>'query',out=>'query.7.out',args=>q[-f'%POS %AN\\n' -i'AN!=2*N_SAMPLES']);
test_vcf_query($opts,in=>'query',out=>'query.8.out',args=>q[-f'%POS[ %GL]\\n' -i'min(abs(GL[*:0]))=10']);
test_vcf_query($opts,in=>'view.filter',out=>'query.9.out',args=>q[-f'%POS %CIGAR\\n' -i'strlen(CIGAR[*])=4']);
test_vcf_query($opts,in=>'query',out=>'query.10.out',args=>q[-f'%POS[ %GT]\\n' -i'AC[0]=3']);
test_vcf_query($opts,in=>'query',out=>'query.10.out',args=>q[-f'%POS[ %GT]\\n' -i'AF[0]=3/4']);
test_vcf_query($opts,in=>'query',out=>'query.11.out',args=>q[-f'%POS[ %GT]\\n' -i'MAC[0]=1']);
test_vcf_query($opts,in=>'query',out=>'query.11.out',args=>q[-f'%POS[ %GT]\\n' -i'MAF[0]=1/4']);
test_vcf_query($opts,in=>'view.vectors',out=>'query.12.out',args=>q[-f'I8=%I8 I16=%I16 I32=%I32 IF=%IF IA8=%IA8 IA16=%IA16 IA32=%IA32 IAF=%IAF IA8=%IA8{1} IA16=%IA16{1} IA32=%IA32{1} IAF=%IAF{1} [ %F8:%F16:%F32:%FF]\\n']);
test_vcf_query($opts,in=>'query.filter',out=>'query.13.out',args=>q[-f'[%POS %SAMPLE %GT\\n]' -i'GT ="1"']);
test_vcf_query($opts,in=>'query.filter',out=>'query.14.out',args=>q[-f'[%POS %SAMPLE %GT\\n]' -i'GT!="1"']);
test_vcf_query($opts,in=>'query.filter',out=>'query.15.out',args=>q[-f'[%POS %SAMPLE %GT\\n]' -e'GT ="1"']);
test_vcf_query($opts,in=>'query.filter',out=>'query.16.out',args=>q[-f'[%POS %SAMPLE %GT\\n]' -e'GT!="1"']);
test_vcf_query($opts,in=>'query.2',out=>'query.17.out',args=>q[-f'%XX_A %XX.A %XX.A0 %xx.a0\\n']);
test_vcf_query($opts,in=>'missing',out=>'query.18.out',args=>q[-i'IINT="."'  -f'%POS %IINT\\n']);
test_vcf_query($opts,in=>'missing',out=>'query.19.out',args=>q[-i'IINT!="."' -f'%POS %IINT\\n']);
test_vcf_query($opts,in=>'missing',out=>'query.18.out',args=>q[-e'IINT!="."' -f'%POS %IINT\\n']);
test_vcf_query($opts,in=>'missing',out=>'query.19.out',args=>q[-e'IINT="."'  -f'%POS %IINT\\n']);
test_vcf_query($opts,in=>'missing',out=>'query.20.out',args=>q[-i'IFLT="."'  -f'%POS %IFLT\\n']);
test_vcf_query($opts,in=>'missing',out=>'query.21.out',args=>q[-i'IFLT!="."' -f'%POS %IFLT\\n']);
test_vcf_query($opts,in=>'missing',out=>'query.20.out',args=>q[-e'IFLT!="."' -f'%POS %IFLT\\n']);
test_vcf_query($opts,in=>'missing',out=>'query.21.out',args=>q[-e'IFLT="."'  -f'%POS %IFLT\\n']);
test_vcf_query($opts,in=>'missing',out=>'query.22.out',args=>q[-i'ISTR="."'  -f'%POS %ISTR\\n']);
test_vcf_query($opts,in=>'missing',out=>'query.23.out',args=>q[-i'ISTR!="."' -f'%POS %ISTR\\n']);
test_vcf_query($opts,in=>'missing',out=>'query.23.out',args=>q[-e'ISTR="."'  -f'%POS %ISTR\\n']);
test_vcf_query($opts,in=>'missing',out=>'query.22.out',args=>q[-e'ISTR!="."' -f'%POS %ISTR\\n']);
test_vcf_query($opts,in=>'missing',out=>'query.24.out',args=>q[-i'FILTER="q11"' -f'%POS %ISTR\\n']);
test_vcf_query($opts,in=>'query',out=>'query.25.out',args=>q[-f'%LINE']);
test_vcf_query($opts,in=>'query.filter-type',out=>'query.26.out',args=>q[-f'%POS\\t%REF\\t%ALT\\n' -i'type="snp"']);
test_vcf_query($opts,in=>'query.filter-type',out=>'query.27.out',args=>q[-f'%POS\\t%REF\\t%ALT\\n' -i'type~"snp"']);
test_vcf_query($opts,in=>'query.filter-type',out=>'query.28.out',args=>q[-f'%POS\\t%REF\\t%ALT\\n' -i'type!="snp"']);
test_vcf_query($opts,in=>'query.filter-type',out=>'query.29.out',args=>q[-f'%POS\\t%REF\\t%ALT\\n' -i'type!~"snp"']);
test_vcf_query($opts,in=>'query.filter-type',out=>'query.67.out',args=>q[-f'%POS\\t%REF\\t%ALT\\n' -i'INFO/TYPE="xxx"']);
test_vcf_query($opts,in=>'filter-missing-floats',out=>'query.30.out',args=>q[-f'%POS\\t%A_AF\\t%B_AF\\t%C_AF\\n' -i'A_AF>=0.0001 || B_AF >= 0.0001 || C_AF >= 0.0001']);
test_vcf_query($opts,in=>'filter-missing-floats',out=>'query.31.out',args=>q[-f'%POS\\t%A_AF\\t%B_AF\\t%C_AF\\n' -e'A_AF>=0.0001 || B_AF >= 0.0001 || C_AF >= 0.0001']);
test_vcf_query($opts,in=>'missing',out=>'query.32.out',args=>q[-i'FMT/FINT!="."' -f'[\t%FINT]\\n']);
test_vcf_query($opts,in=>'query.filter.2',out=>'query.33.out',args=>q[-f'[%GT]\\n' -i'GT~"0/[1-9]" || GT~"[1-9]/0"']);
test_vcf_query($opts,in=>'view.filter',out=>'query.34.out',args=>q[-f'[%POS %SAMPLE %FGS\\n]\\n' -i'FGS[*:1-2,4]="EE"']);
test_vcf_query($opts,in=>'view.filter',out=>'query.34.out',args=>q[-f'[%POS %SAMPLE %FGS\\n]\\n' -i'"EE"=FGS[*:1-2,4]']);
test_vcf_query($opts,in=>'view.filter',out=>'query.35.out',args=>q[-f'[%POS %SAMPLE %FGI\\n]\\n' -i'FGI[*:1-2,5]=6']);
test_vcf_query($opts,in=>'view.filter',out=>'query.35.out',args=>q[-f'[%POS %SAMPLE %FGI\\n]\\n' -i'6=FGI[*:1-2,5]']);
test_vcf_query($opts,in=>'view.filter',out=>'query.36.out',args=>q[-f'[%POS %SAMPLE %FGF\\n]\\n' -i'FGF[*:1-3,4]=5e-5']);
test_vcf_query($opts,in=>'view.filter',out=>'query.37.out',args=>q[-f'[%POS %SAMPLE %FGS\\n]\\n' -i'FGS[*:4]!="EE"']);
test_vcf_query($opts,in=>'view.filter',out=>'query.38.out',args=>q[-f'[%POS %SAMPLE %FGS\\n]\\n' -i'FGS[*:4]="."']);
test_vcf_query($opts,in=>'view.filter',out=>'query.39.out',args=>q[-f'[%POS %SAMPLE %FGS\\n]\\n' -i'FGS[*:4]!="."']);
test_vcf_query($opts,in=>'view.filter',out=>'query.40.out',args=>q[-f'[%POS %SAMPLE %FGI\\n]\\n' -i'FGI[*:1-2,5]="."']);
test_vcf_query($opts,in=>'view.filter',out=>'query.41.out',args=>q[-f'[%POS %SAMPLE %FGI\\n]\\n' -i'FGI[*:1-2,5]!="."']);
test_vcf_query($opts,in=>'filter.5',out=>'query.42.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP=1 | FMT/DP=2']);
test_vcf_query($opts,in=>'filter.5',out=>'query.43.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP=1 | FMT/DP="."']);
test_vcf_query($opts,in=>'view.filter',out=>'query.44.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP=19']);
test_vcf_query($opts,in=>'view.filter',out=>'query.45.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP=19 |  FMT/DP="."']);
test_vcf_query($opts,in=>'view.filter',out=>'query.46.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP=19 || FMT/DP="."']);
test_vcf_query($opts,in=>'view.filter',out=>'query.47.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP[0]=19']);
test_vcf_query($opts,in=>'view.filter',out=>'query.48.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP[0]!=19']);
test_vcf_query($opts,in=>'view.filter',out=>'query.49.out',args=>q[-f'[%POS %SAMPLE %DP %GQ\\n]\\n' -i'FMT/DP=19 &  FMT/GQ=589']);
test_vcf_query($opts,in=>'view.filter',out=>'query.50.out',args=>q[-f'[%POS %SAMPLE %DP %GQ\\n]\\n' -i'FMT/DP=19 && FMT/GQ=1']);
test_vcf_query($opts,in=>'query.filter.3',out=>'query.51.out',args=>q[-f'[\\t%GT\\n]\\n' -i'GT~"1" && GT~"2"']);
test_vcf_query($opts,in=>'query.filter.3',out=>'query.52.out',args=>q[-f'[\\t%GT\\n]\\n' -i'GT~"1" &  GT~"2"']);
test_vcf_query($opts,in=>'query.filter.3',out=>'query.53.out',args=>q[-f'%POS[\\t%GT]\\n' -i'COUNT(GT="het")=1']);
test_vcf_query($opts,in=>'filter.5',out=>'query.54.out',args=>q[-f'[%POS  %SAMPLE  %AD\\n]\\n' -i'AD[:0]+AD[:1] > 12']);
test_vcf_query($opts,in=>'query.filter.4',out=>'query.55.out',args=>q[-f'%POS\\t%REF\\t%ALT[\\t%GT]\\n' -e'TYPE!="snp" || ALT="*"']);
test_vcf_query($opts,in=>'view',out=>'query.56.out',args=>q[-f'%ID\\n' -i 'ID=@].$$opts{path}.q[/query.56.out']);
test_vcf_query($opts,in=>'query.filter.5',out=>'query.57.out',args=>q[-f'[%POS\\t%SAMPLE\\t%GT\\t%AD\\n]' -i'GT="het" & binom(FMT/AD)>0.01']);
test_vcf_query($opts,in=>'query.filter.5',out=>'query.58.out',args=>q[-f'[%POS\\t%SAMPLE\\t%GT\\t%AD\\n]' -i'GT="het" & binom(FMT/AD[:0],FMT/AD[:1])>0.01']);
test_vcf_query($opts,in=>'query.filter.5',out=>'query.59.out',args=>q[-f'%POS\\t%AD\\n' -i'binom(INFO/AD[0],INFO/AD[1])>0.01']);
test_vcf_query($opts,in=>'query',out=>'query.60.out',args=>q[-f'%CHROM %POS\\n' -i'CHROM="4"']);
test_vcf_query($opts,in=>'query.negative',out=>'query.61.out',args=>q[-f'%POS\\t%TAG1\\n' -i'(TAG1>=-129 && TAG1<=-120) || (TAG1>=-32769 && TAG1<=-32760)']);
test_vcf_query($opts,in=>'query.negative',out=>'query.61.out',args=>q[-f'%POS\\t%TAGV1\\n' -i'(TAGV1>=-129 && TAGV1<=-120) || (TAGV1>=-32769 && TAGV1<=-32760)']);
test_vcf_query($opts,in=>'query.negative',out=>'query.62.out',args=>q[-f'%POS\\t%TAG2\\n' -i'(TAG2>=-129 && TAG2<=-120) || (TAG2>=-32769 && TAG2<=-32760)']);
test_vcf_query($opts,in=>'query.negative',out=>'query.62.out',args=>q[-f'%POS\\t%TAGV2\\n' -i'(TAGV2>=-129 && TAGV2<=-120) || (TAGV2>=-32769 && TAGV2<=-32760)']);
test_vcf_query($opts,in=>'query',out=>'query.63.out',args=>q[-f'[%POS\\t%SAMPLE\\t%GQ\\n]' -i'N_PASS(GQ<20)==1']);
test_vcf_query($opts,in=>'query',out=>'query.64.out',args=>q[-f'%CHROM\\t%POS\\t%INFO\\t%FORMAT\\n' -s D,C]);
test_vcf_query($opts,in=>'query.pbinom.1',out=>'query.65.out',args=>q[-f'[%POS %SAMPLE %GT %AD %PBINOM(AD)\\n]' -i'phred(binom(FMT/AD))>=0']);
test_vcf_query($opts,in=>'query.filter.6',out=>'query.66.out',args=>q[-f'%POS\\n' -i'POS==16777217 || POS==33554432 || POS=118673904']);
test_vcf_query($opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'sum(II)==6']);
test_vcf_query($opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'sum(FORMAT/FI)==7']);
test_vcf_query($opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'avg(II)==2']);
test_vcf_query($opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'avg(FORMAT/FI)==1.75']);
test_vcf_query($opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'mean(II)==2']);
test_vcf_query($opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'mean(FORMAT/FI)==1.75']);
test_vcf_query($opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'median(II)==2']);
test_vcf_query($opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'median(FORMAT/FI)==1.5']);
test_vcf_query($opts,in=>'query.filter.8',out=>'query.69.out',args=>q[-f'%POS\\t%REF\\t%ALT\\t%ILEN\\n' -i'%ILEN==1']);
test_vcf_query($opts,in=>'query.filter.8',out=>'query.70.out',args=>q[-f'%POS\\t%REF\\t%ALT\\t%ILEN\\n' -i'ILEN==1']);
test_vcf_query($opts,in=>'query.filter.9',out=>'query.71.out',args=>q[-f'[%POS  %SAMPLE  %AD\\n]' -i'FMT/AD[:0] < FMT/AD[:1]']);
test_vcf_query($opts,in=>'query.filter.9',out=>'query.72.out',args=>q[-f'[%POS  %SAMPLE  %AD\\n]' -i'FMT/AD[:0] > FMT/AD[:1]']);
test_vcf_query($opts,in=>'query.filter.10',out=>'query.73.out',args=>q[-f'%POS  %NUM_TAG\\n' -i'COUNT(INFO/NUM_TAG)=2']);
test_vcf_query($opts,in=>'query.filter.10',out=>'query.74.out',args=>q[-f'%POS  %STR_TAG\\n' -i'COUNT(INFO/STR_TAG)=2']);
test_vcf_norm($opts,in=>'norm',out=>'norm.out',fai=>'norm',args=>'-cx');
test_vcf_norm($opts,in=>'norm.split',out=>'norm.split.out',args=>'-m-');
test_vcf_norm($opts,in=>'norm.split.2',out=>'norm.split.2.out',args=>'-m-');
test_vcf_norm($opts,in=>'norm.split.3',out=>'norm.split.3.out',args=>'-m- --force');
test_vcf_norm($opts,in=>'norm.split',fai=>'norm',out=>'norm.split.and.norm.out',args=>'-m-');
test_vcf_norm($opts,in=>'norm.merge',out=>'norm.merge.out',args=>'-m+');
test_vcf_norm($opts,in=>'norm.merge.2',out=>'norm.merge.2.out',args=>'-m+');
test_vcf_norm($opts,in=>'norm.merge.3',out=>'norm.merge.3.out',args=>'-m+');
test_vcf_norm($opts,in=>'norm.merge',out=>'norm.merge.strict.out',args=>'-m+ -s');
test_vcf_norm($opts,in=>'norm.setref',out=>'norm.setref.out',args=>'-Nc s',fai=>'norm');
test_vcf_norm($opts,in=>'norm.telomere',out=>'norm.telomere.out',fai=>'norm');
test_vcf_norm($opts,in=>'norm.rmdup',out=>'norm.rmdup.1.out',args=>'-d snps');
test_vcf_norm($opts,in=>'norm.rmdup',out=>'norm.rmdup.2.out',args=>'-d indels');
test_vcf_norm($opts,in=>'norm.rmdup',out=>'norm.rmdup.3.out',args=>'-d both');
test_vcf_norm($opts,in=>'norm.rmdup',out=>'norm.rmdup.4.out',args=>'-d all');
test_vcf_norm($opts,in=>'norm.rmdup',out=>'norm.rmdup.5.out',args=>'-d none');
test_vcf_norm($opts,in=>'norm.rmdup',out=>'norm.rmdup.5.out',args=>'-d exact');
test_vcf_norm($opts,in=>'norm.rmdup.2',out=>'norm.rmdup.2.1.out',args=>'-d none');
test_vcf_norm($opts,in=>'norm.rmdup.2',out=>'norm.rmdup.2.1.out',args=>'-d exact');
test_vcf_norm($opts,in=>'norm.rmdup.2',out=>'norm.rmdup.2.1.out',args=>'-d indels');
test_vcf_norm($opts,in=>'norm.rmdup.2',out=>'norm.rmdup.2.2.out',args=>'-d any');
test_vcf_norm($opts,in=>'norm.rmdup.2',out=>'norm.rmdup.2.2.out',args=>'-d both');
test_vcf_norm($opts,in=>'norm.rmdup.2',out=>'norm.rmdup.2.2.out',args=>'-d snps');
test_vcf_norm($opts,in=>'norm.2',fai=>'norm.2',out=>'norm.2.out',args=>'');
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
test_vcf_view($opts,in=>'view.omitgenotypes',out=>'view.dropgenotypes.out',args=>'-G',reg=>'');
test_vcf_view($opts,in=>'view.omitgenotypes',out=>'view.dropgenotypes.noheader.out',args=>'-HG',reg=>'');
test_vcf_view($opts,in=>'many.alleles',out=>'many.alleles.trim.out',args=>'-a',reg=>'');
test_vcf_view($opts,in=>'view.vectors',out=>'view.vectors.A.out',args=>'-asA',reg=>'');
test_vcf_view($opts,in=>'view.vectors',out=>'view.vectors.B.out',args=>'-asB',reg=>'');
test_vcf_view($opts,in=>'view.vectors.2',out=>'view.vectors.C.out',args=>'-asA',reg=>'');
test_vcf_view($opts,in=>'view.filter',out=>'view.filter.1.out',args=>q[-H -i'FMT/FGS[*:0]="AAAAAA"'],reg=>'');    # test expressions
test_vcf_view($opts,in=>'view.filter',out=>'view.filter.2.out',args=>q[-H -i'FMT/FGS[*:2]="C"'],reg=>'');
test_vcf_view($opts,in=>'view.filter',out=>'view.filter.3.out',args=>q[-H -i'FMT/FGS[*:4]="EE"'],reg=>'');
test_vcf_view($opts,in=>'view.filter',out=>'view.filter.4.out',args=>q[-H -i'FMT/FRS[*:1]="BB"'],reg=>'');
test_vcf_view($opts,in=>'view.filter',out=>'view.filter.5.out',args=>q[-H -i'TXT0="text"'],reg=>'');
test_vcf_view($opts,in=>'view.chrs',out=>'view.chrs.out',args=>'',reg=>'',tgts=>'view.chrs.tab');
test_vcf_view($opts,in=>'filter.2',out=>'filter.11.out',args=>q[-i 'POS>=3062917'],reg=>'1:3062917-3157410');
test_vcf_view($opts,in=>'idx.1',out=>'idx.1.out',args=>q[-H -r 1:10,1:12,1:10]);
test_vcf_view($opts,in=>'idx.2',out=>'idx.2.out',args=>q[-H -r 1:1172777-1172804,1:1172806-1172808]);
test_vcf_view($opts,in=>'idx.2',out=>'idx.2.out',args=>q[-H -R {PATH}/idx.2.bed]);
test_vcf_view($opts,in=>'idx.3',out=>'idx.3.out',args=>q[-H -R {PATH}/idx.3.bed]);
test_vcf_view($opts,in=>'idx.4',out=>'idx.4.out',args=>q[-H -R {PATH}/idx.4.bed]);
test_vcf_64bit($opts,in=>'view64bit.1',out=>'view64bit.1.out',do_bcf=>1);
test_vcf_64bit($opts,in=>'view64bit.2',out=>'view64bit.2.out',do_bcf=>1);
test_vcf_64bit($opts,in=>'view64bit.3',out=>'view64bit.3.out');     # large coordinates don't work with BCF
test_vcf_64bit($opts,in=>'view64bit.4',out=>'view64bit.4.out',do_bcf=>1);
test_vcf_64bit($opts,in=>'view64bit.5',out=>'view64bit.5.out',do_bcf=>1);
test_vcf_filter($opts,in=>'view.filter',out=>'view.filter.6.out',args=>q[-S. -e'TXT0="text"'],reg=>'');
test_vcf_filter($opts,in=>'view.filter',out=>'view.filter.7.out',args=>q[-S. -e'FMT/FRS[*:1]="BB"'],reg=>'');
test_vcf_filter($opts,in=>'view.filter',out=>'view.filter.8.out',args=>q[-S. -e'FMT/FGS[*:0]="AAAAAA"'],reg=>'');
test_vcf_filter($opts,in=>'view.filter',out=>'view.filter.9.out',args=>q[-S. -e'FMT/FGS[*:1]="BBB"'],reg=>'');
test_vcf_filter($opts,in=>'view.filter',out=>'view.filter.10.out',args=>q[-S. -e'FMT/FGS[*:4]="EE"'],reg=>'');
test_vcf_filter($opts,in=>'view.filter',out=>'view.filter.11.out',args=>q[-S. -e'FMT/STR="XX"'],reg=>'');
test_vcf_filter($opts,in=>'view.filter.2',out=>'view.filter.12.out',args=>q[-S. -e'FMT/FILTER="aaa"'],reg=>'');
test_vcf_view($opts,in=>'view.minmaxac',out=>'view.minmaxac.1.out',args=>q[-H -C5:nonmajor],reg=>'');
test_vcf_view($opts,in=>'view.minmaxac',out=>'view.minmaxac.2.out',args=>q[-H -c6:nonmajor],reg=>'');
test_vcf_view($opts,in=>'view.minmaxac',out=>'view.minmaxac.1.out',args=>q[-H -q0.3:major],reg=>'');
test_vcf_view($opts,in=>'view.filter.annovar',out=>'view.filter.annovar.1.out',args=>q[-H -i 'Gene.refGene=="RAD21L1"'],reg=>'');
test_vcf_view($opts,in=>'view.filter.annovar',out=>'view.filter.annovar.2.out',args=>q[-H -i 'Gene.refGene~"NOD"'],reg=>'');
test_vcf_view($opts,in=>'view.filter.annovar',out=>'view.filter.annovar.3.out',args=>q[-H -i 'LJB2_MutationTaster=="0.291000"'],reg=>'');
test_vcf_call($opts,in=>'mpileup',out=>'mpileup.1.out',args=>'-mv');
test_vcf_call($opts,in=>'mpileup',out=>'mpileup.2.out',args=>'-mg0');
test_vcf_call($opts,in=>'mpileup',out=>'mpileup.3.out',args=>'-mv -S {PATH}/mpileup.3.samples');
test_vcf_call($opts,in=>'mpileup',out=>'mpileup.4.out',args=>'-mv -S {PATH}/mpileup.4.samples');
test_vcf_call($opts,in=>'mpileup',out=>'mpileup.5.out',args=>'-mv -S {PATH}/mpileup.5.samples');
test_vcf_call($opts,in=>'mpileup.X',out=>'mpileup.X.out',args=>'-mv --ploidy-file {PATH}/mpileup.ploidy -S {PATH}/mpileup.samples');
test_vcf_call($opts,in=>'mpileup.X',out=>'mpileup.X.out',args=>'-mv --ploidy-file {PATH}/mpileup.ploidy -S {PATH}/mpileup.ped');
test_vcf_call($opts,in=>'mpileup.X',out=>'mpileup.X.2.out',args=>'-mv --ploidy-file {PATH}/mpileup.ploidy -S {PATH}/mpileup.2.samples');
test_vcf_call($opts,in=>'mpileup.NA19213.NA19129',out=>'mpileup.hwe.1.out',args=>'-mv');
test_vcf_call($opts,in=>'mpileup.NA19213.NA19129',out=>'mpileup.hwe.1.out',args=>'-mv -G -');
test_vcf_call($opts,in=>'mpileup.hwe',out=>'mpileup.hwe.2.out',args=>'-mv');
test_vcf_call($opts,in=>'mpileup.hwe',out=>'mpileup.hwe.3.out',args=>'-mv -G -');
test_vcf_call($opts,in=>'mpileup.hwe',out=>'mpileup.hwe.4.out',args=>'-mv -G {PATH}/mpileup.hwe.samples');
test_vcf_call_cAls($opts,in=>'mpileup',out=>'mpileup.cAls.out',tab=>'mpileup');
test_vcf_call_cAls($opts,in=>'mpileup.2',out=>'mpileup.cAls.2.out',tab=>'mpileup.2');
test_vcf_call_cAls($opts,in=>'mpileup.3',out=>'mpileup.cAls.3.out',tab=>'mpileup.3',args=>'-i');
test_vcf_call_cAls($opts,in=>'mpileup.3',out=>'mpileup.cAls.4.out',tab=>'mpileup.4',args=>'-i');
test_vcf_call_cAls($opts,in=>'mpileup.3',out=>'mpileup.cAls.5.out',tab=>'mpileup.5',args=>'-i');
test_vcf_call_cAls($opts,in=>'mpileup.4',out=>'mpileup.cAls.6.out',tab=>'mpileup.6',args=>'-i');
test_vcf_call_cAls($opts,in=>'mpileup.5',out=>'mpileup.cAls.7.out',tab=>'mpileup.7',args=>'-i');
test_vcf_call_cAls($opts,in=>'mpileup.cals.1',out=>'mpileup.cals.8.out',tab=>'mpileup.cals.1',args=>'');
test_vcf_call_cAls($opts,in=>'mpileup.cals.2',out=>'mpileup.cals.9.out',tab=>'mpileup.cals.2',args=>'');
test_vcf_call($opts,in=>'mpileup.c',out=>'mpileup.c.1.out',args=>'-cv');
# test_vcf_call($opts,in=>'mpileup.c',out=>'mpileup.c.2.out',args=>'-cg0');
test_vcf_call($opts,in=>'mpileup.c.X',out=>'mpileup.c.X.out',args=>'-cv --ploidy-file {PATH}/mpileup.ploidy -S {PATH}/mpileup.samples');
test_vcf_call($opts,in=>'mpileup.c.X',out=>'mpileup.c.X.out',args=>'-cv --ploidy-file {PATH}/mpileup.ploidy -S {PATH}/mpileup.ped');
test_vcf_call($opts,in=>'mpileup.c.X',out=>'mpileup.c.X.2.out',args=>'-cv --ploidy-file {PATH}/mpileup.ploidy -S {PATH}/mpileup.2.samples');
test_vcf_filter($opts,in=>'filter.1',out=>'filter.1.out',args=>'-mx -g2 -G2');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.2.out',args=>q[-e'QUAL==59.2 || (INDEL=0 & (FMT/GQ=25 | FMT/DP=10))' -sModified -S.]);
test_vcf_filter($opts,in=>'filter.3',out=>'filter.3.out',args=>q[-e'INFO/DP=19'],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.4.out',args=>q[-e'INFO/DP=19' -s XX],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.5.out',args=>q[-e'INFO/DP=19' -s XX -m+],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.6.out',args=>q[-e'INFO/DP=19' -s XX -mx],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.7.out',args=>q[-e'INFO/DP=19' -s XX -m+x],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.3.out',args=>q[-e'FMT/GT="0/2"'],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.4.out',args=>q[-e'FMT/GT="0/2"' -s XX],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.5.out',args=>q[-e'FMT/GT="0/2"' -s XX -m+],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.6.out',args=>q[-e'FMT/GT="0/2"' -s XX -mx],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.3',out=>'filter.7.out',args=>q[-e'FMT/GT="0/2"' -s XX -m+x],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.8.out',args=>q[-i'FMT/GT="0/0" && AC[*]=2'],fmt=>'%POS\\t%AC[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.8.out',args=>q[-i'AC[*]=2 && FMT/GT="0/0"'],fmt=>'%POS\\t%AC[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.9.out',args=>q[-i'ALT="."'],fmt=>'%POS\\t%AC[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.4',out=>'filter.10.out',args=>q[-S . -i 'FORMAT/TEST3<25']);
test_vcf_filter($opts,in=>'filter.4',out=>'filter.10.out',args=>q[-S . -i 'FORMAT/TEST4<25']);
test_vcf_filter($opts,in=>'filter.2',out=>'filter.12.out',args=>q[-i'GT="A"'],fmt=>'%POS[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.13.out',args=>q[-i'GT="RR"'],fmt=>'%POS[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.14.out',args=>q[-i'GT="RA"'],fmt=>'%POS[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.14.out',args=>q[-i'GT="AR"'],fmt=>'%POS[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.15.out',args=>q[-i'GT="AA"'],fmt=>'%POS[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.16.out',args=>q[-i'GT="aA"'],fmt=>'%POS[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.16.out',args=>q[-i'GT="Aa"'],fmt=>'%POS[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.17.out',args=>q[-i'GT="HOM"'],fmt=>'%POS[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.18.out',args=>q[-i'GT="HET"'],fmt=>'%POS[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.19.out',args=>q[-i'GT="HAP"'],fmt=>'%POS[\\t%GT]\\n');
test_vcf_filter($opts,in=>'filter.5',out=>'filter.20.out',args=>q[-i'AD[:1]=11'],fmt=>'%POS[\\t%AD]\\n');
test_vcf_filter($opts,in=>'filter.5',out=>'filter.21.out',args=>q[-i'AD[1:]=11'],fmt=>'%POS[\\t%AD]\\n');
test_vcf_filter($opts,in=>'filter.5',out=>'filter.22.out',args=>q[-i'FR[0:1]=11'],fmt=>'%POS[\\t%FR]\\n');
test_vcf_filter($opts,in=>'filter.5',out=>'filter.23.out',args=>q[-i'AD[*]="."'],fmt=>'%POS[\\t%AD]\\n');
test_vcf_filter($opts,in=>'filter.5',out=>'filter.24.out',args=>q[-i'AD[0:0]=="."'],fmt=>'%POS[\\t%AD]\\n');
test_vcf_filter($opts,in=>'filter.5',out=>'filter.25.out',args=>q[-i'AD[0:0]!="."'],fmt=>'%POS[\\t%AD]\\n');
test_vcf_filter($opts,in=>'filter.5',out=>'filter.26.out',args=>q[-i'QUAL=="."'],fmt=>'%POS\\t%QUAL\\n');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.27.out',args=>q[-i'N_PASS(DP>32)=1'],fmt=>'[%POS\\t%SAMPLE\\t%DP\\n]');
test_vcf_filter($opts,in=>'filter.2',out=>'filter.27.out',args=>q[-i'F_PASS(DP>32)=0.5'],fmt=>'[%POS\\t%SAMPLE\\t%DP\\n]');
test_vcf_filter($opts,in=>'filter.6',out=>'filter.28.out',args=>q[-i'F_MISSING>=1/5'],fmt=>'%POS\\n');
test_vcf_filter($opts,in=>'filter.6',out=>'filter.28.out',args=>q[-i'F_MISSING>=0.2'],fmt=>'%POS\\n');
test_vcf_filter($opts,in=>'filter.6',out=>'filter.28.out',args=>q[-i'F_PASS(GT=="mis")>=1/5'],fmt=>'%POS\\n');
test_vcf_filter($opts,in=>'filter.6',out=>'filter.28.out',args=>q[-i'F_PASS(GT=="mis")>=0.2'],fmt=>'%POS\\n');
test_vcf_sort($opts,in=>'sort',out=>'sort.out',args=>q[-m 0],fmt=>'%CHROM\\t%POS\\t%REF,%ALT\\n');
test_vcf_sort($opts,in=>'sort',out=>'sort.out',args=>q[-m 1000],fmt=>'%CHROM\\t%POS\\t%REF,%ALT\\n');
test_vcf_regions($opts,in=>'regions');
test_vcf_annotate($opts,in=>'annotate',tab=>'annotate',out=>'annotate.out',args=>'-c CHROM,POS,REF,ALT,ID,QUAL,INFO/T_INT,INFO/T_FLOAT,INDEL');
test_vcf_annotate($opts,in=>'annotate',tab=>'annotate2',out=>'annotate2.out',args=>'-c CHROM,POS,-,T_STR');
test_vcf_annotate($opts,in=>'annotate',tab=>'annotate2',out=>'annotate22.out',args=>'-c CHROM,FROM,TO,T_STR');
test_vcf_annotate($opts,in=>'annotate',vcf=>'annots',out=>'annotate3.out',args=>'-c STR,ID,QUAL,FILTER');
test_vcf_annotate($opts,in=>'annotate2',vcf=>'annots2',out=>'annotate4.out',args=>'-c ID,QUAL,FILTER,INFO,FMT');
test_vcf_annotate($opts,in=>'annotate2',vcf=>'annots2',out=>'annotate5.out',args=>'-c ID,QUAL,+FILTER,+INFO,FMT/GT -s A');
test_vcf_annotate($opts,in=>'annotate2',vcf=>'annots2',out=>'annotate18.out',args=>'-c ID,QUAL,+FILTER,+INFO,FMT/GT -s "A B"');
test_vcf_annotate($opts,in=>'annotate2',vcf=>'annots2',out=>'annotate19.out',args=>'-c ID,QUAL,+FILTER,+INFO,FMT/GT -s "A C"');
test_vcf_annotate($opts,in=>'annotate2',vcf=>'annots2',out=>'annotate20.out',args=>'-c ID,QUAL,+FILTER,+INFO,FMT/GT -s "B C"');
test_vcf_annotate($opts,in=>'annotate3',out=>'annotate6.out',args=>'-x ID,QUAL,^FILTER/fltA,FILTER/fltB,^INFO/AA,INFO/BB,^FMT/GT,FMT/PL');
test_vcf_annotate($opts,in=>'annotate3',out=>'annotate7.out',args=>'-x FORMAT');
test_vcf_annotate($opts,in=>'annotate4',vcf=>'annots4',out=>'annotate8.out',args=>'-c +INFO');
test_vcf_annotate($opts,in=>'annotate4',tab=>'annots4',out=>'annotate8.out',args=>'-c CHROM,POS,REF,ALT,+FA,+FR,+IA,+IR,+SA,+SR');
test_vcf_annotate($opts,in=>'annotate10',tab=>'annots10',out=>'annotate10.out',args=>'-c CHROM,POS,FMT/FINT,FMT/FFLT,FMT/FSTR');
test_vcf_annotate($opts,in=>'annotate2',vcf=>'annots2',out=>'annotate11.out',args=>'-c CHROM,POS,FMT/FINT,FMT/FFLT,FMT/FSTR -s A');
test_vcf_annotate($opts,in=>'annotate2',tab=>'annots11',out=>'annotate11.out',args=>'-c CHROM,POS,FMT/FINT,FMT/FFLT,FMT/FSTR -s A');
test_vcf_annotate($opts,in=>'annotate2',vcf=>'annots2',out=>'annotate12.out',args=>'-c AAA:=IINT,FMT/BBB:=FMT/FINT');
test_vcf_annotate($opts,in=>'annotate2',vcf=>'annots2',out=>'annotate13.out',args=>'-x INFO -c INFO/IINT');
test_vcf_annotate($opts,in=>'annotate2',vcf=>'annots2',out=>'annotate14.out',args=>q[-x INFO -c INFO/IINT -e'POS=3000001' -k]);
test_vcf_annotate($opts,in=>'annotate11',vcf=>'annots11',out=>'annotate15.out',args=>q[-c FMT]);
test_vcf_annotate($opts,in=>'annotate2',vcf=>'annots2',out=>'annotate16.out',args=>'-c FMT/newGT:=GT');
test_vcf_annotate($opts,in=>'annotate2',vcf=>'annots12',out=>'annotate17.out',args=>'-c FMT/GT:=newGT');
test_vcf_annotate($opts,in=>'annotate13',tab=>'annots13',out=>'annotate21.out',args=>'-c CHROM,BEG,END,ABC');
test_vcf_annotate($opts,in=>'annotate13',tab=>'annots13',out=>'annotate23.out',args=>'-c CHROM,BEG,END,ABC -l ABC:append');
test_vcf_annotate($opts,in=>'annotate13',tab=>'annots13',out=>'annotate24.out',args=>'-c CHROM,BEG,END,ABC -l ABC:unique');
test_vcf_annotate($opts,in=>'annotate14',out=>'annotate25.out',args=>'-x FILTER/XX,INFO/XX --force');
test_vcf_annotate($opts,in=>'annotate15',tab=>'annotate15',out=>'annotate26.out',args=>'-s SAMPLE1 -c CHROM,FROM,TO,FMT/FOO,BAR');
test_vcf_annotate($opts,in=>'annotate15',tab=>'annotate15',out=>'annotate27.out',args=>'-s SAMPLE2 -c CHROM,FROM,TO,FMT/FOO,BAR');
test_vcf_annotate($opts,in=>'annotate16',out=>'annotate28.out',args=>'-x FILTER');
test_vcf_annotate($opts,in=>'annotate17.1',tab=>'annotate17.1',out=>'annotate17.1.out',args=>'-c CHROM,BEG,END,A,B -l A:append,B:append');
test_vcf_annotate($opts,in=>'annotate17.2',tab=>'annotate17.1',out=>'annotate17.2.out',args=>'-c CHROM,BEG,END,A,B -l A:append,B:append');
test_vcf_annotate($opts,in=>'annotate17.3',tab=>'annotate17.3',out=>'annotate17.3.out',args=>'-c CHROM,BEG,END,A,B -l A:append,B:append');
test_vcf_annotate($opts,in=>'annotate18.1',tab=>'annotate18.1',out=>'annotate18.1.out',args=>'-c CHROM,BEG,END,A,B,C,D,E -l A:sum,B:avg,C:min,D:max,E:append');
test_vcf_annotate($opts,in=>'annotate18.2',tab=>'annotate18.2',out=>'annotate18.2.out',args=>'-c CHROM,BEG,END,A,B,C,D,E -l A:sum,B:avg,C:min,D:max,E:append');
test_vcf_plugin($opts,in=>'plugin1',out=>'missing2ref.out',cmd=>'+missing2ref --no-version');
test_vcf_plugin($opts,in=>'plugin1',out=>'missing2ref.out',cmd=>'+setGT --no-version',args=>'-- -t . -n 0');
test_vcf_plugin($opts,in=>'setGT',out=>'setGT.1.out',cmd=>'+setGT --no-version',args=>'-- -t q -n 0 -i \'GT~"." && FMT/DP=30 && GQ=150\'');
test_vcf_plugin($opts,in=>'setGT.2',out=>'setGT.2.out',cmd=>'+setGT --no-version',args=>'-- -t q -n . -i \'GT[@{QPATH}/setGT.samples.txt]="het"\'');
test_vcf_plugin($opts,in=>'setGT.2',out=>'setGT.3.out',cmd=>'+setGT --no-version',args=>'-- -t q -n . -i \'GT[@{QPATH}/setGT.samples.txt]="het" & binom(AD[@{QPATH}/setGT.samples.txt])<0.1\'');
test_vcf_annotate($opts,in=>'annotate9',tab=>'annots9',out=>'annotate9.out',args=>'-c CHROM,POS,REF,ALT,+ID');
test_vcf_plugin($opts,in=>'plugin1',out=>'fill-AN-AC.out',cmd=>'+fill-AN-AC --no-version');
test_vcf_plugin($opts,in=>'dosage',out=>'dosage.1.out',cmd=>'+dosage',args=>'-- -t PL');
test_vcf_plugin($opts,in=>'dosage',out=>'dosage.2.out',cmd=>'+dosage',args=>'-- -t GL');
test_vcf_plugin($opts,in=>'dosage',out=>'dosage.3.out',cmd=>'+dosage',args=>'-- -t GT');
test_vcf_plugin($opts,in=>'fixploidy',out=>'fixploidy.out',cmd=>'+fixploidy --no-version',args=>'-- -s {PATH}/fixploidy.samples -p {PATH}/fixploidy.ploidy');
test_vcf_plugin($opts,in=>'view.PL',out=>'guess-ploidy.PL.out',cmd=>'+guess-ploidy',args=>'-vrX | grep -v bcftools');
test_vcf_plugin($opts,in=>'view.GL',out=>'guess-ploidy.GL.out',cmd=>'+guess-ploidy',args=>'-vrX | grep -v bcftools');
test_vcf_plugin($opts,in=>'view.GL',out=>'view.PL.vcf',cmd=>'+tag2tag --no-version',args=>'-- -r --gl-to-pl');
test_vcf_plugin($opts,in=>'view.GP',out=>'view.GT.vcf',cmd=>'+tag2tag --no-version',args=>'-- -r --gp-to-gt -t 0.2');
test_vcf_plugin($opts,in=>'query.variantkey',out=>'query.add-variantkey.vcf',cmd=>'+add-variantkey',args=>'');
test_vcf_plugin($opts,in=>'query.variantkey',out=>'variantkey-hex.out',cmd=>'+variantkey-hex',args=>'test/');
test_vcf_plugin($opts,in=>'query.nucleotide',out=>'query.allele-length.tsv',cmd=>'+allele-length',args=>'');
test_vcf_plugin($opts,in=>'merge.a',out=>'fill-tags.out',cmd=>'+fill-tags --no-version',args=>'-- -t AN,AC,AC_Hom,AC_Het,AC_Hemi');
test_vcf_plugin($opts,in=>'view',out=>'fill-tags.2.out',cmd=>'+fill-tags --no-version',args=>'-- -t AC,AN,AF,MAF,NS');
test_vcf_plugin($opts,in=>'view',out=>'fill-tags.3.out',cmd=>'+fill-tags --no-version',args=>'-- -t AC -S {PATH}/fill-tags.3.smpl');
test_vcf_plugin($opts,in=>'view',out=>'fill-tags.5.out',cmd=>'+fill-tags --no-version',args=>'-- -t "DP=sum(DP)" -S {PATH}/fill-tags.3.smpl');
test_vcf_plugin($opts,in=>'many-alts',out=>'fill-tags.4.out',cmd=>'+fill-tags --no-version',args=>'-- -t AN,AC');
test_vcf_plugin($opts,in=>'fill-tags-hemi',out=>'fill-tags-hemi.1.out',cmd=>'+fill-tags --no-version');
test_vcf_plugin($opts,in=>'fill-tags-hemi',out=>'fill-tags-hemi.2.out',cmd=>'+fill-tags --no-version',args=>'-- -d');
test_vcf_plugin($opts,in=>'fill-tags-hwe',out=>'fill-tags-hwe.out',cmd=>'+fill-tags --no-version');
test_vcf_plugin($opts,in=>'fill-tags-AN0',out=>'fill-tags-AN0.out',cmd=>'+fill-tags --no-version');
test_vcf_plugin($opts,in=>'view',out=>'view.GTisec.out',cmd=>'+GTisec',args=>' | grep -v bcftools');
test_vcf_plugin($opts,in=>'view',out=>'view.GTisec.H.out',cmd=>'+GTisec',args=>'-- -H | grep -v bcftools');
test_vcf_plugin($opts,in=>'view',out=>'view.GTisec.Hm.out',cmd=>'+GTisec',args=>'-- -Hm | grep -v bcftools');
test_vcf_plugin($opts,in=>'view',out=>'view.GTisec.Hmv.out',cmd=>'+GTisec',args=>'-- -Hmv | grep -v bcftools');
test_vcf_plugin($opts,in=>'view',out=>'view.GTisec.Hv.out',cmd=>'+GTisec',args=>'-- -Hv | grep -v bcftools');
test_vcf_plugin($opts,in=>'view',out=>'view.GTisec.m.out',cmd=>'+GTisec',args=>'-- -m | grep -v bcftools');
test_vcf_plugin($opts,in=>'view',out=>'view.GTisec.mv.out',cmd=>'+GTisec',args=>'-- -mv | grep -v bcftools');
test_vcf_plugin($opts,in=>'view',out=>'view.GTisec.v.out',cmd=>'+GTisec',args=>'-- -v | grep -v bcftools');
test_vcf_plugin($opts,in=>'trio',out=>'trio.out',cmd=>'+trio-switch-rate',args=>'-- -p {PATH}/trio.ped | grep -v bcftools');
test_vcf_plugin($opts,in=>'trio-stats',out=>'trio-stats.out',cmd=>'+trio-stats',args=>'-a 1 -p {PATH}/trio-stats.ped -d mendel-errors,transmitted | grep -v ^CMD');
test_vcf_plugin($opts,in=>'trio-stats',out=>'trio-stats.2.out',cmd=>'+trio-stats',args=>'-p {PATH}/trio-stats.ped -d mendel-errors,transmitted | grep -v ^CMD');
test_vcf_plugin($opts,in=>'indel-stats',out=>'smpl-stats.1.out',cmd=>'+smpl-stats',args=>'| grep -v ^CMD');
test_vcf_plugin($opts,in=>'indel-stats',out=>'indel-stats.1.out',cmd=>'+indel-stats',args=>'| grep -v ^CMD');
test_vcf_plugin($opts,in=>'indel-stats',out=>'indel-stats.2.out',cmd=>'+indel-stats',args=>' -p {PATH}/trio-stats.ped | grep -v ^CMD');
test_vcf_plugin($opts,in=>'indel-stats',out=>'indel-stats.3.out',cmd=>'+indel-stats',args=>' -p {PATH}/trio-stats.2.ped | grep -v ^CMD');
test_vcf_plugin($opts,in=>'ad-bias',out=>'ad-bias.out',cmd=>'+ad-bias',args=>'-- -s {PATH}/ad-bias.samples | grep -v bcftools');
test_vcf_plugin($opts,in=>'af-dist',out=>'af-dist.out',cmd=>'+af-dist',args=>' | grep -v bcftools');
test_vcf_plugin($opts,in=>'fixref.2a',out=>'fixref.2.out',index=>['fixref.2b'],cmd=>'+fixref',args=>'-- -f {PATH}/norm.fa -i {TMP}/fixref.2b.vcf.gz');
test_vcf_plugin($opts,in=>'fixref.3',out=>'fixref.3.out',cmd=>'+fixref',args=>'-- -f {PATH}/fixref.3.fa -m top');
test_vcf_plugin($opts,in=>'fixref.2a',out=>'fixref.4.out',index=>['fixref.2b'],cmd=>'+fixref',args=>'-- -f {PATH}/norm.fa -m ref-alt');
test_vcf_plugin($opts,in=>'fixref.2a',out=>'fixref.5.out',index=>['fixref.2b'],cmd=>'+fixref',args=>'-- -f {PATH}/norm.fa -m flip');
test_vcf_plugin($opts,in=>'aa',out=>'aa.out',cmd=>'+fill-from-fasta',args=>'-- -f {PATH}/aa.fa -c AA -h {PATH}/aa.hdr -i \'TYPE="snp"\'');
test_vcf_plugin($opts,in=>'aa',out=>'aa.2.out',cmd=>'+fill-from-fasta',args=>'-- -f {PATH}/aa.fa -c REF -N');
test_vcf_plugin($opts,in=>'ref',out=>'ref.out',cmd=>'+fill-from-fasta',args=>'-- -f {PATH}/norm.fa -c REF');
test_vcf_plugin($opts,in=>'view',out=>'view.GTsubset.NA1.out',cmd=>'+GTsubset --no-version',args=>'-- -s NA00001');
test_vcf_plugin($opts,in=>'view',out=>'view.GTsubset.NA1NA2.out',cmd=>'+GTsubset --no-version',args=>'-- -s NA00001,NA00002');
test_vcf_plugin($opts,in=>'view',out=>'view.GTsubset.NA1NA2NA3.out',cmd=>'+GTsubset --no-version',args=>'-- -s NA00001,NA00002,NA00003');
test_vcf_plugin($opts,in=>'mendelian',out=>'mendelian.1.out',cmd=>'+mendelian',args=>'-t mom1,dad1,child1 -d');
test_vcf_plugin($opts,in=>'mendelian',out=>'mendelian.2.out',cmd=>'+mendelian',args=>'-t mom1,dad1,child1 -l+');
test_vcf_plugin($opts,in=>'mendelian',out=>'mendelian.3.out',cmd=>'+mendelian',args=>'-t mom1,dad1,child1 -lx');
test_vcf_plugin($opts,in=>'contrast',out=>'contrast.out',cmd=>'+contrast',args=>'-0 a,b -1 c');
test_vcf_plugin($opts,in=>'contrast',out=>'contrast.out',cmd=>'+contrast',args=>'-0 {PATH}/contrast0.txt -1 {PATH}/contrast1.txt');
test_vcf_plugin($opts,in=>'trio-dnm.1',out=>'trio-dnm.1.out',cmd=>'+trio-dnm',args=>"-p proband,father,mother | $$opts{bin}/bcftools query -f'%CHROM[\\t%DNM]\\t[\\t%VAF]\\n'");
test_vcf_plugin($opts,in=>'trio-dnm.2',out=>'trio-dnm.1.out',cmd=>'+trio-dnm',args=>"-p proband,father,mother --force-AD | $$opts{bin}/bcftools query -f'%CHROM[\\t%DNM]\\t[\\t%VAF]\\n'");
test_vcf_plugin($opts,in=>'trio-dnm.2',out=>'trio-dnm.2.out',cmd=>'+trio-dnm',args=>"-p proband,father,mother | $$opts{bin}/bcftools query -f'%CHROM[\\t%DNM]\\t[\\t%VAF]\\n'");
test_vcf_plugin($opts,in=>'gvcfz',out=>'gvcfz.1.out',cmd=>'+gvcfz',args=>qq[-g 'PASS:GT!="alt"' -a | $$opts{bin}/bcftools query -f'%POS\\t%REF\\t%ALT\\t%END[\\t%GT][\\t%DP][\\t%GQ][\\t%RGQ]\\n']);
test_vcf_plugin($opts,in=>'gvcfz',out=>'gvcfz.2.out',cmd=>'+gvcfz',args=>qq[-g 'PASS:GQ>10; FLT:-' -a | $$opts{bin}/bcftools query -f'%POS\\t%REF\\t%ALT\\t%FILTER\\t%END[\\t%GT][\\t%DP][\\t%GQ][\\t%RGQ]\\n']);
test_vcf_plugin($opts,in=>'gvcfz.2',out=>'gvcfz.2.1.out',cmd=>'+gvcfz',args=>qq[-g 'PASS:GT!="alt"' -a | $$opts{bin}/bcftools query -f'%POS\\t%REF\\t%ALT\\t%FILTER\\t%END[\\t%GT][\\t%DP]\\n']);
test_vcf_plugin($opts,in=>'remove-overlaps',out=>'remove-overlaps.1.out',cmd=>'+remove-overlaps',args=>'');
test_vcf_plugin($opts,in=>'remove-overlaps',out=>'remove-overlaps.2.out',cmd=>'+remove-overlaps',args=>'-d');
test_vcf_plugin($opts,in=>'split-vep',out=>'split-vep.1.out',cmd=>'+split-vep',args=>qq[-c Consequence -s worst:missense+ | $$opts{bin}/bcftools query -f'%POS\\t%Consequence\\n']);
test_vcf_plugin($opts,in=>'split-vep',out=>'split-vep.2.out',cmd=>'+split-vep',args=>qq[-c Consequence -s worst:missense+ | $$opts{bin}/bcftools query -f'%POS\\t%Consequence\\n' -i'Consequence!="."']);
test_vcf_plugin($opts,in=>'split-vep',out=>'split-vep.2.out',cmd=>'+split-vep',args=>qq[-s worst:missense+ -f'%POS\\t%Consequence\\n']);
test_vcf_plugin($opts,in=>'split-vep',out=>'split-vep.3.out',cmd=>'+split-vep',args=>qq[-s primary:missense+ -f'%POS\\t%Consequence\\n']);
test_vcf_plugin($opts,in=>'split-vep',out=>'split-vep.4.out',cmd=>'+split-vep',args=>qq[-s primary:missense+ -f'%POS\\n']);
test_vcf_plugin($opts,in=>'split-vep.2',out=>'split-vep.5.out',cmd=>'+split-vep',args=>qq[-s worst -f'%POS\\t%AF\\n']);
test_vcf_plugin($opts,in=>'split-vep.2',out=>'split-vep.6.out',cmd=>'+split-vep',args=>qq[-s worst -f'%POS\\t%INFO/AF\\n']);
test_vcf_plugin($opts,in=>'split-vep.3',out=>'split-vep.7.out',cmd=>'+split-vep',args=>qq[-s worst -f'%POS\\t%Consequence\\n']);
test_vcf_plugin($opts,in=>'split-vep.3',out=>'split-vep.8.out',cmd=>'+split-vep',args=>qq[-s worst -f'[%POS\\t%SAMPLE\\t%GT\\t%Consequence\\n]' -i'GT="alt"']);
test_vcf_plugin($opts,in=>'split-vep',out=>'split-vep.9.out', cmd=>'+split-vep',args=>qq[-t 1:14464 -f '%POS\\t%CANONICAL\\t%Consequence\\n']);
test_vcf_plugin($opts,in=>'split-vep',out=>'split-vep.10.out',cmd=>'+split-vep',args=>qq[-t 1:14464 -f '%POS\\t%CANONICAL\\t%Consequence\\n' -d]);
test_vcf_plugin($opts,in=>'split-vep',out=>'split-vep.11.out',cmd=>'+split-vep',args=>qq[-t 1:14464 -f '%POS\\t%CSQ\\n' -A tab]);
test_vcf_plugin($opts,in=>'split-vep',out=>'split-vep.12.out',cmd=>'+split-vep',args=>qq[-t 1:14464 -f '%POS\\t%CSQ\\n' -A tab -d]);
test_vcf_plugin($opts,in=>'split-vep.4',out=>'split-vep.13.out',cmd=>'+split-vep',args=>qq[-f '%POS\\t%BCSQ\\n' -a BCSQ -A tab -d]);
test_vcf_plugin($opts,in=>'split-vep',out=>'split-vep.14.out',cmd=>'+split-vep',args=>qq[-c gnomAD_NFE_AF:real,ALLELE_NUM:int | $$opts{bin}/bcftools query -f'%POS\\t%gnomAD_NFE_AF\\t%ALLELE_NUM\\n']);
test_vcf_plugin($opts,in=>'split-vep',out=>'split-vep.14.out',cmd=>'+split-vep',args=>qq[-c gnomAD_NFE_AF,ALLELE_NUM          | $$opts{bin}/bcftools query -f'%POS\\t%gnomAD_NFE_AF\\t%ALLELE_NUM\\n']);
test_vcf_plugin($opts,in=>'split-vep.5',out=>'split-vep.15.out',cmd=>'+split-vep',args=>qq[-s :synonymous    -c Consequence | $$opts{bin}/bcftools query -f'%POS\\t%Consequence\\n']);
test_vcf_plugin($opts,in=>'split-vep.5',out=>'split-vep.16.out',cmd=>'+split-vep',args=>qq[-s :synonymous -x -c Consequence | $$opts{bin}/bcftools query -f'%POS\\t%Consequence\\n']);
test_vcf_plugin($opts,in=>'parental-origin',out=>'parental-origin.1.out',cmd=>'+parental-origin',args=>qq[-r 20:100 -p proband,father,mother -t del | grep -v ^#]);
test_vcf_plugin($opts,in=>'parental-origin',out=>'parental-origin.2.out',cmd=>'+parental-origin',args=>qq[-r 20:101 -p proband,father,mother -t del | grep -v ^#]);
test_vcf_plugin($opts,in=>'parental-origin',out=>'parental-origin.3.out',cmd=>'+parental-origin',args=>qq[-r 20:102 -p proband,father,mother -t del | grep -v ^#]);
test_vcf_plugin($opts,in=>'parental-origin',out=>'parental-origin.4.out',cmd=>'+parental-origin',args=>qq[-r 20:103 -p proband,father,mother -t dup | grep -v ^#]);
test_vcf_plugin($opts,in=>'parental-origin',out=>'parental-origin.5.out',cmd=>'+parental-origin',args=>qq[-r 20:104 -p proband,father,mother -t dup | grep -v ^#]);
test_vcf_concat($opts,in=>['concat.1.a','concat.1.b'],out=>'concat.1.vcf.out',do_bcf=>0,args=>'');
test_vcf_concat($opts,in=>['concat.1.a','concat.1.b'],out=>'concat.1.bcf.out',do_bcf=>1,args=>'');
test_vcf_concat($opts,in=>['concat.2.a','concat.2.b'],out=>'concat.2.vcf.out',do_bcf=>0,args=>'-a');
test_vcf_concat($opts,in=>['concat.2.a','concat.2.b'],out=>'concat.2.bcf.out',do_bcf=>1,args=>'-a');
test_vcf_concat($opts,in=>['concat.2.a','concat.2.b'],out=>'concat.4.vcf.out',do_bcf=>0,args=>'-aD');
test_vcf_concat($opts,in=>['concat.2.a','concat.2.b'],out=>'concat.4.bcf.out',do_bcf=>1,args=>'-aD');
test_vcf_concat($opts,in=>['concat.3.a','concat.3.b','concat.3.0','concat.3.c','concat.3.d','concat.3.e','concat.3.f'],out=>'concat.3.vcf.out',do_bcf=>0,args=>'-l');
test_vcf_concat($opts,in=>['concat.3.a','concat.3.b','concat.3.0','concat.3.c','concat.3.d','concat.3.e','concat.3.f'],out=>'concat.3.bcf.out',do_bcf=>1,args=>'-l');
test_naive_concat($opts,name=>'naive_concat',max_hdr_lines=>10000,max_body_lines=>10000,nfiles=>10);
test_vcf_reheader($opts,in=>'reheader',out=>'reheader.1.out',header=>'reheader.hdr');
test_vcf_reheader($opts,in=>'reheader',out=>'reheader.2.out',samples=>'reheader.samples');
test_vcf_reheader($opts,in=>'reheader',out=>'reheader.2.out',samples=>'reheader.samples2');
test_vcf_reheader($opts,in=>'reheader',out=>'reheader.3.out',samples=>'reheader.samples3');
test_vcf_reheader($opts,in=>'reheader',out=>'reheader.4.out',samples=>'reheader.samples4');
test_vcf_reheader($opts,in=>'empty',out=>'reheader.empty.out',header=>'reheader.empty.hdr');
test_vcf_reheader($opts,in=>'reheader.2',out=>'reheader.5.out',args=>'-f {PATH}/reheader.fai',nostdin=>1);
test_vcf_reheader($opts,in=>'reheader.2',out=>'reheader.5.out',args=>'-h {PATH}/reheader.2.hdr -f {PATH}/reheader.fai',nostdin=>1);
test_rename_chrs($opts,in=>'annotate');
test_vcf_convert($opts,in=>'convert',out=>'convert.gs.gt.gen',args=>'-g -,.');
test_vcf_convert($opts,in=>'convert',out=>'convert.gs.gt.samples',args=>'-g .,-');
test_vcf_convert($opts,in=>'convert',out=>'convert.gs.pl.gen',args=>'-g -,. --tag PL');
test_vcf_convert($opts,in=>'convert',out=>'convert.gs.pl.samples',args=>'-g .,- --tag PL');
test_vcf_convert($opts,in=>'check',out=>'check.gs.vcfids.gen',args=>'-g -,. --vcf-ids');
test_vcf_convert($opts,in=>'check',out=>'check.gs.vcfids.samples',args=>'-g .,- --vcf-ids');
test_vcf_convert($opts,in=>'check',out=>'check.gs.chrom.gen',args=>'-g -,. --chrom');
test_vcf_convert($opts,in=>'check',out=>'check.gs.chrom.samples',args=>'-g .,- --chrom');
test_vcf_convert($opts,in=>'check',out=>'check.gs.vcfids_chrom.gen',args=>'-g -,. --chrom --vcf-ids');
test_vcf_convert($opts,in=>'check',out=>'check.gs.vcfids_chrom.samples',args=>'-g .,- --chrom --vcf-ids');
test_vcf_convert($opts,in=>'convert',out=>'convert.hls.haps',args=>'-h -,.,.');
test_vcf_convert($opts,in=>'convert',out=>'convert.hls.legend',args=>'-h .,-,.');
test_vcf_convert($opts,in=>'convert',out=>'convert.hls.samples',args=>'-h .,.,-');
test_vcf_convert_hls2vcf($opts,h=>'convert.hls.gt.hap',l=>'convert.hls.gt.legend',s=>'convert.hls.gt.samples',out=>'convert.gt.noHead.vcf',args=>'-H');
test_vcf_convert_hs2vcf($opts,h=>'convert.hs.gt.hap',s=>'convert.hs.gt.samples',out=>'convert.gt.noHead.vcf',args=>'--hapsample2vcf');
test_vcf_convert($opts,in=>'convert',out=>'convert.hs.hap',args=>'--hapsample -,.');
test_vcf_convert($opts,in=>'convert',out=>'convert.hs.sample',args=>'--hapsample .,-');
test_vcf_convert($opts,in=>'convert.hap-missing',out=>'convert.hap-missing.haps',args=>'--haplegendsample -,.,.');
test_vcf_convert_gvcf($opts,in=>'convert.gvcf',out=>'convert.gvcf.out',fa=>'gvcf.fa',args=>'--gvcf2vcf -i\'FILTER="PASS"\'');
test_vcf_convert_tsv2vcf($opts,in=>'convert.23andme',out=>'convert.23andme.vcf',args=>'-c ID,CHROM,POS,AA -s SAMPLE1',fai=>'23andme');
test_vcf_consensus($opts,in=>'consensus',out=>'consensus.1.out',fa=>'consensus.fa',mask=>'consensus.tab',args=>'');
test_vcf_consensus_chain($opts,in=>'consensus',out=>'consensus.1.chain',chain=>'consensus.1.chain',fa=>'consensus.fa',mask=>'consensus.tab',args=>'');
test_vcf_consensus($opts,in=>'consensus',out=>'consensus.2.out',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-H 1');
test_vcf_consensus_chain($opts,in=>'consensus',out=>'consensus.2.chain',chain=>'consensus.2.chain',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-H 1');
test_vcf_consensus($opts,in=>'consensus',out=>'consensus.3.out',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-I');
test_vcf_consensus_chain($opts,in=>'consensus',out=>'consensus.3.chain',chain=>'consensus.3.chain',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-I');
test_vcf_consensus($opts,in=>'consensus',out=>'consensus.4.out',fa=>'consensus.fa',args=>'-H 1');
test_vcf_consensus_chain($opts,in=>'consensus',out=>'consensus.4.chain',chain=>'consensus.4.chain',fa=>'consensus.fa',args=>'-H 1');
test_vcf_consensus($opts,in=>'consensus2',out=>'consensus2.1.out',fa=>'consensus2.fa',args=>'-H 1');
test_vcf_consensus($opts,in=>'consensus2',out=>'consensus2.2.out',fa=>'consensus2.fa',args=>'-H 2');
test_vcf_consensus($opts,in=>'empty',out=>'consensus.5.out',fa=>'consensus.fa',args=>'');
test_vcf_consensus($opts,in=>'consensus3',out=>'consensus3.out',fa=>'consensus2.fa',args=>'-H 2 -M "?"');
test_vcf_consensus($opts,in=>'consensus3',out=>'consensus3.2.out',fa=>'consensus2.fa',args=>'-H 2 -M "?" -p xx_');
test_vcf_consensus($opts,in=>'consensus4',out=>'consensus4.out',fa=>'consensus2.fa',args=>'');
test_vcf_consensus($opts,in=>'consensus5',out=>'consensus5.out',fa=>'consensus5.fa',args=>'--haplotype LA');
test_vcf_consensus($opts,in=>'consensus6',out=>'consensus6.out',fa=>'consensus6.fa',args=>'');
test_vcf_consensus($opts,in=>'consensus7',out=>'consensus7a.out',fa=>'consensus7.fa',args=>'-H 2');
test_vcf_consensus($opts,in=>'consensus7',out=>'consensus7b.out',fa=>'consensus7.fa',args=>'-H 2pIu');
test_vcf_consensus($opts,in=>'consensus7',out=>'consensus7c.out',fa=>'consensus7.fa',args=>'-H 1');
test_vcf_consensus($opts,in=>'consensus7',out=>'consensus7d.out',fa=>'consensus7.fa',args=>'-H 1pIu');
test_mpileup($opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.1.out',args=>q[-r17:100-150],test_list=>1);
test_mpileup($opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.2.out',args=>q[-a DP,DV -r17:100-600]); # test files from samtools mpileup test suite
test_mpileup($opts,in=>[qw(mpileup.1)],out=>'mpileup/mpileup.3.out',args=>q[-B --ff 0x14 -r17:1050-1060]); # test file converted to vcf from samtools mpileup test suite
test_mpileup($opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.4.out',args=>q[-a DP,DPR,DV,DP4,INFO/DPR,SP -r17:100-600]); #test files from samtools mpileup test suite
test_mpileup($opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.5.out',args=>q[-a DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -r17:100-600]);
test_mpileup($opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.6.out',args=>q[-a DP,DV -r17:100-600 --gvcf 0,2,5]);
test_mpileup($opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.6.out',args=>q[-a DP,DV -r17:100-200,17:201-300,17:301-400,17:401-500,17:501-600 --gvcf 0,2,5]);
test_mpileup($opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.7.out',args=>q[-r17:100-150 -s HG00101,HG00102]);
test_mpileup($opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.7.out',args=>q[-r17:100-150 -S {PATH}/mplp.samples]);
test_mpileup($opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.8.out',args=>q[-r17:100-150 -s ^HG00101,HG00102]);
test_mpileup($opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.8.out',args=>q[-r17:100-150 -S ^{PATH}/mplp.samples]);
test_mpileup($opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.9.out',args=>q[-t17:100-150 -S {PATH}/mplp.9.samples]);
test_mpileup($opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.10.out',args=>q[-t17:100-150 -G {PATH}/mplp.10.samples]);
test_mpileup($opts,in=>[qw(mpileup.3)],out=>'mpileup/mpileup.11.out',args=>q[]);
test_mpileup($opts,in=>[qw(mpileup.3 mpileup.4)],out=>'mpileup/mpileup.11.out',args=>q[-s HG00102]);
test_mpileup($opts,in=>[qw(mpileup.3 mpileup.4)],out=>'mpileup/mpileup.11.out',args=>q[-s ^HG99999]);
test_mpileup($opts,in=>[qw(mpileup.3 mpileup.4)],out=>'mpileup/mpileup.11.out',args=>q[-G {PATH}/mplp.11.rgs]);
test_mpileup($opts,in=>[qw(mpileup.3 mpileup.4)],out=>'mpileup/mpileup.11.out',args=>q[-G {PATH}/mplp.11.rgs]);
test_mpileup($opts,in=>[qw(indel-AD.1)],out=>'mpileup/indel-AD.1.out',ref=>'indel-AD.1.fa',args=>q[-a AD]);
test_mpileup($opts,in=>[qw(mpileup-SCR)],out=>'mpileup/mpileup-SCR.out',ref=>'mpileup-SCR.fa',args=>q[-a INFO/SCR,FMT/SCR]);
test_csq($opts,in=>'csq',out=>'csq.1.out',cmd=>'-f {PATH}/csq.fa -g {PATH}/csq.gff3');
test_csq_real($opts,in=>'csq');


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

sub cygpath {
    my ($path) = @_;
    $path = `cygpath -m $path`;
    $path =~ s/\r?\n//;
    return $path
}

sub safe_tempdir
{
    my $dir = tempdir(CLEANUP=>1);
    if ($^O =~ /^msys/) {
        $dir = cygpath($dir);
    }
    return $dir;
}

sub parse_params
{
    my $opts = { bgzip=>"bgzip", keep_files=>0, nok=>0, nfailed=>0, tabix=>"tabix", plugins=>0, htsdir=>"" };
    my $help;
    Getopt::Long::Configure('bundling');
    my $ret = GetOptions (
            'e|exec=s' => sub { my ($tool, $path) = split /=/, $_[1]; $$opts{$tool} = $path if $path },
            't|temp-dir:s' => \$$opts{keep_files},
            'p|plugins' => \$$opts{test_plugins},
            'r|redo-outputs' => \$$opts{redo_outputs},
            'h|?|help' => \$help,
            'i|htsdir:s' => \$$opts{htsdir}
            );
    if ( !$ret or $help ) { error(); }
    if ( $$opts{htsdir} ) {
        if ($^O eq 'cygwin' || $^O =~ /^msys/) {
            # Set PATH so against-htslib-source builds can find the htslib dll
            $ENV{PATH} = "$$opts{htsdir}:"."$$opts{htsdir}/bin:"."$$opts{htsdir}/lib:".$ENV{PATH};
        }
    }
    $$opts{tmp} = $$opts{keep_files} ? $$opts{keep_files} : safe_tempdir;
    if ( $$opts{keep_files} ) { cmd("mkdir -p $$opts{keep_files}"); }
    $$opts{path} = $FindBin::RealBin;
    $$opts{bin}  = $FindBin::RealBin;
    $$opts{bin}  =~ s{/test/?$}{};
    if ($^O =~ /^msys/) {
        $$opts{path} = cygpath($$opts{path});
        $$opts{bin}  = cygpath($$opts{bin});
    }
    
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
sub _cmd3
{
    my ($cmd) = @_;

    my $tmp = "$$opts{tmp}/tmp";
    my $pid = fork();
    if ( !$pid )
    {
        exec('/bin/bash', '-o','pipefail','-c', "($cmd) 2>$tmp.e >$tmp.o");
    }
    waitpid($pid,0);

    my $status  = $? >> 8;
    my $signal  = $? & 127;

    my (@out,@err);
    if ( open(my $fh,'<',"$tmp.o") )
    {
        @out = <$fh>;
        close($fh) or error("Failed to close $tmp.o");
    }
    if ( open(my $fh,'<',"$tmp.e") )
    {
        @err = <$fh>;
        close($fh) or error("Failed to close $tmp.e");
    }
    unlink("$tmp.o");
    unlink("$tmp.e");

    return ($status,join('',@out),join('',@err));
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

    my ($ret,$out,$err) = _cmd3("$args{cmd}");
    if ( length($err) ) { $err =~ s/\n/\n\t\t/gs; $err = "\n\n\t\t$err\n"; }
    if ( $ret ) { failed($opts,$test,"Non-zero status $ret$err"); return; }
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
    if ( exists($args{exp}) ) { $exp = $args{exp}; }
    elsif ( open(my $fh,'<',"$$opts{path}/$args{out}") )
    {
        my @exp = <$fh>;
        $exp = join('',@exp);
        close($fh);
    }
    else
    {
        open(my $fh,'>',"$$opts{path}/$args{out}.new") or error("$$opts{path}/$args{out}.new: $!");
        print $fh $out;
        close($fh);
        if ( !$$opts{redo_outputs} ) { failed($opts,$test,"$$opts{path}/$args{out}.new"); return; }
    }

    $exp =~ s/\r\n/\n/g; # if ($args{text});  ## simplified as all tests are text?
    $out =~ s/\r\n/\n/g; # if ($args{text});  ## simplified as all tests are text?
    $out =~ s/e([-+])0(\d\d)/e$1$2/g if ($args{exp_fix});
    if ( $exp ne $out )
    {
        open(my $fh,'>',"$$opts{path}/$args{out}.new") or error("$$opts{path}/$args{out}.new");
        print $fh $out;
        close($fh);
        if ( !-e "$$opts{path}/$args{out}" && !exists($args{exp}) )
        {
            rename("$$opts{path}/$args{out}.new","$$opts{path}/$args{out}") or error("rename $$opts{path}/$args{out}.new $$opts{path}/$args{out}: $!");
            print "\tthe file with expected output does not exist, creating new one:\n";
            print "\t\t$$opts{path}/$args{out}\n";
        }
        else
        {
            if ( exists($args{exp}) )
            {
                open(my $fh,'>',"$$opts{path}/$args{out}") or error("$$opts{path}/$args{out}");
                print $fh $exp;
                close($fh);
            }
            failed($opts,$test,"The outputs differ:\n\t\t$$opts{path}/$args{out}\n\t\t$$opts{path}/$args{out}.new$err");
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

    # output path
    unlink("$$opts{tmp}/$args{in}.bcf.csi", "$$opts{tmp}/$args{in}.bcf.csi", "$$opts{tmp}/$args{in}.vcf.gz.tbi");
    cmd("$$opts{bin}/bcftools index -fo $$opts{tmp}/$args{in}.csi $$opts{tmp}/$args{in}.bcf");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -H $$opts{tmp}/$args{in}.bcf $args{reg}");

    # streaming
    cmd("$$opts{bin}/bcftools view -Oz $$opts{path}/$args{in}.vcf | tee $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools index -fo $$opts{tmp}/$args{in}.vcf.gz.csi");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -H $$opts{tmp}/$args{in}.vcf.gz $args{reg}");

    cmd("$$opts{bin}/bcftools view -Ob $$opts{path}/$args{in}.vcf | tee $$opts{tmp}/$args{in}.bcf | $$opts{bin}/bcftools index -fo $$opts{tmp}/$args{in}.bcf.csi");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -H $$opts{tmp}/$args{in}.bcf $args{reg}");
}

sub test_vcf_idxstats
{
    my ($opts,%args) = @_;
    cmd("$$opts{bin}/bcftools view -Oz $$opts{path}/$args{in}.vcf > $$opts{tmp}/$args{in}.vcf.gz");
    cmd("$$opts{bin}/bcftools index --tbi -f $$opts{tmp}/$args{in}.vcf.gz");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools index $args{args} $$opts{tmp}/$args{in}.vcf.gz");
    unlink("$$opts{tmp}/$args{in}.vcf.gz.tbi");
    cmd("$$opts{bin}/bcftools index --csi -f $$opts{tmp}/$args{in}.vcf.gz");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools index $args{args} $$opts{tmp}/$args{in}.vcf.gz");
    unlink("$$opts{tmp}/$args{in}.vcf.gz.csi");

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

sub test_vcf_check_merge
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    cmd("$$opts{bin}/bcftools stats -r 1 $$opts{tmp}/$args{in}.vcf.gz > $$opts{tmp}/$args{in}.1.chk");
    cmd("$$opts{bin}/bcftools stats -r 2 $$opts{tmp}/$args{in}.vcf.gz > $$opts{tmp}/$args{in}.2.chk");
    cmd("$$opts{bin}/bcftools stats -r 3 $$opts{tmp}/$args{in}.vcf.gz > $$opts{tmp}/$args{in}.3.chk");
    cmd("$$opts{bin}/bcftools stats -r 4 $$opts{tmp}/$args{in}.vcf.gz > $$opts{tmp}/$args{in}.4.chk");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/misc/plot-vcfstats -m $$opts{tmp}/$args{in}.1.chk $$opts{tmp}/$args{in}.2.chk $$opts{tmp}/$args{in}.3.chk $$opts{tmp}/$args{in}.4.chk | grep -v 'plot-vcfstats' | grep -v '^# The command' | grep -v '^# This' | grep -v '^ID\t'");
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
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools merge --no-version $args $files", exp_fix=>1);
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools merge -Ob $args $files | $$opts{bin}/bcftools view | grep -v ^##bcftools_", exp_fix => 1);
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
    $args{args} =~ s/{PATH}/$$opts{path}/g;
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
    $args{args} =~ s/{PATH}/$$opts{path}/g;
    bgzip_tabix($opts,file=>$args{tab_in},suffix=>'tab',args=>'-s 1 -b 2 -e 3');
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools isec --no-version  $args{args} -T $$opts{tmp}/$args{tab_in}.tab.gz $files");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools isec -Ob $args{args} -T $$opts{tmp}/$args{tab_in}.tab.gz $files | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
}
sub test_vcf_query
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools query $args{args} $$opts{tmp}/$args{in}.vcf.gz", exp_fix=>1);
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools query $args{args}", exp_fix=>1);
}
sub test_vcf_convert
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools convert $args{args} $$opts{tmp}/$args{in}.vcf.gz");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools convert $args{args}");
}
sub test_vcf_convert_hls2vcf
{
    my ($opts,%args) = @_;
    my $hls = join(',', map { "$$opts{path}/$_" }( $args{h}, $args{l}, $args{s} ) );
    if($^O =~ /^msys/) { $hls =~ s/\//\\\\/g; }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools convert $args{args} $hls | grep -v ^##");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools convert $args{args} $hls -Ou | $$opts{bin}/bcftools view | grep -v ^##");
}
sub test_vcf_convert_hs2vcf
{
    my ($opts,%args) = @_;
    my $hs = join(',', map { "$$opts{path}/$_" }( $args{h}, $args{s} ) );
    if($^O =~ /^msys/) { $hs =~ s/\//\\\\/g; }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools convert $args{args} $hs | grep -v ^##");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools convert $args{args} $hs -Ou | $$opts{bin}/bcftools view | grep -v ^##");
}
sub test_vcf_convert_gvcf
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools convert --no-version $args{args} -f $$opts{path}/$args{fa} $$opts{tmp}/$args{in}.vcf.gz");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools convert $args{args} -f $$opts{path}/$args{fa} | grep -v ^##bcftools");
}
sub test_vcf_convert_tsv2vcf
{
    my ($opts,%args) = @_;
    my $params = '';
    if ( exists($args{args}) ) { $params .= " $args{args}"; }
    if ( exists($args{fai} ) ) { $params .= " -f $$opts{path}/$args{fai}.fa"; }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools convert --no-version $params --tsv2vcf $$opts{path}/$args{in}");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools convert -Ou $params --tsv2vcf $$opts{path}/$args{in} | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
}
sub test_vcf_norm
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    my $params = '';
    if ( exists($args{args}) ) { $params .= " $args{args}"; }
    if ( exists($args{fai} ) ) { $params .= " -f $$opts{path}/$args{fai}.fa"; }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools norm --no-version $params $$opts{tmp}/$args{in}.vcf.gz",exp_fix=>1);
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools norm -Ob $params $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools view | grep -v ^##bcftools_",exp_fix=>1);
}
sub test_vcf_view
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});

    if ( !exists($args{args}) ) { $args{args} = ''; }
    if ( !exists($args{reg}) ) { $args{reg} = ''; }
    if ( exists($args{tgts}) ) { $args{args} .= "-T $$opts{path}/$args{tgts}"; }
    $args{args} =~ s/{PATH}/$$opts{path}/g;
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view --no-version $args{args} $$opts{tmp}/$args{in}.vcf.gz $args{reg}", exp_fix=>1);
    unless ($args{args} =~ /-H/) {
        test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -Ob $args{args} $$opts{tmp}/$args{in}.vcf.gz $args{reg} | $$opts{bin}/bcftools view | grep -v ^##bcftools_", exp_fix=>1);
    }
}
sub test_vcf_64bit
{
    my ($opts,%args) = @_;
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view $$opts{path}/$args{in}.vcf -H", exp_fix=>1);
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view $$opts{path}/$args{in}.vcf     | $$opts{bin}/bcftools view -H", exp_fix=>1);
    if ( $args{do_bcf} )
    {
        test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view $$opts{path}/$args{in}.vcf -Ou | $$opts{bin}/bcftools view -H", exp_fix=>1);
    }
}
sub test_vcf_call
{
    my ($opts,%args) = @_;
    $args{args} =~ s/{PATH}/$$opts{path}/g;
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools call --no-version $args{args} $$opts{path}/$args{in}.vcf");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools call -Ob $args{args} $$opts{path}/$args{in}.vcf | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
}
sub test_vcf_call_cAls
{
    my ($opts,%args) = @_;
    my $args = exists($args{args}) ? $args{args} : '';
    bgzip_tabix($opts,file=>$args{tab},suffix=>'tab',args=>'-s1 -b2 -e2');
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools call --no-version -mA -C alleles -T $$opts{tmp}/$args{tab}.tab.gz $args $$opts{path}/$args{in}.vcf");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools call -Ob -mA -C alleles -T $$opts{tmp}/$args{tab}.tab.gz $args $$opts{path}/$args{in}.vcf | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
}
sub test_vcf_filter
{
    my ($opts,%args) = @_;
    my $pipe = 'grep -v ^##bcftools_';
    if ( exists($args{fmt}) )
    {
        $pipe = "$$opts{bin}/bcftools query -f '$args{fmt}'";
    }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools filter $args{args} $$opts{path}/$args{in}.vcf | $pipe", exp_fix=>1);
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools filter -Ob $args{args} $$opts{path}/$args{in}.vcf | $$opts{bin}/bcftools view | $pipe", exp_fix=>1);
}
sub test_vcf_sort
{
    my ($opts,%args) = @_;
    my $pipe = "$$opts{bin}/bcftools query -f '$args{fmt}'";
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools sort $args{args} $$opts{path}/$args{in}.vcf | $pipe");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools sort -Ob $args{args} $$opts{path}/$args{in}.vcf | $$opts{bin}/bcftools view | $pipe");
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
    # Under msys the isatty function fails to recognise the terminal.
    # Skip these tests for now.
    return if ($^O =~ /^msys/);

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
        $hdr = -e "$$opts{path}/$args{in}.hdr" ? "-h $$opts{path}/$args{in}.hdr" : '';
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
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools annotate $annot_fname $hdr $args{args} $in_fname | $$opts{bin}/bcftools view | grep -v ^##bcftools_", exp_fix=>1);
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools annotate -Ob $annot_fname $hdr $args{args} $in_fname | $$opts{bin}/bcftools view | grep -v ^##bcftools_", exp_fix=>1);
}
sub test_vcf_plugin
{
    my ($opts,%args) = @_;
    if ( !$$opts{test_plugins} ) { return; }
    $ENV{BCFTOOLS_PLUGINS} = "$$opts{bin}/plugins";
    if ( !exists($args{args}) ) { $args{args} = ''; }
    my $wpath = $$opts{path}; 
    if ($^O =~ /^msys/) {
        $wpath = `cygpath -w $$opts{path}`;
        $wpath =~ s/\r?\n//;
    }
    $args{args} =~ s/{QPATH}/$wpath/g;
    $args{args} =~ s/{PATH}/$$opts{path}/g;
    $args{cmd}  =~ s/{PATH}/$$opts{path}/g;
    $args{args} =~ s/{TMP}/$$opts{tmp}/g;
    $args{cmd}  =~ s/{TMP}/$$opts{tmp}/g;
    bgzip_tabix_vcf($opts,"$args{in}");
    if ( exists($args{index}) )
    {
        for my $file (@{$args{index}}) { bgzip_tabix_vcf($opts,$file); }
    }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools $args{cmd} $$opts{tmp}/$args{in}.vcf.gz $args{args} | grep -v ^##bcftools_");

    cmd("$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz > $$opts{tmp}/$args{in}.bcf");
    cmd("$$opts{bin}/bcftools index -f $$opts{tmp}/$args{in}.bcf");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools $args{cmd} $$opts{tmp}/$args{in}.bcf $args{args} | grep -v ^##bcftools_", exp_fix=>1);
}
sub test_vcf_concat
{
    my ($opts,%args) = @_;
    my $files;
    for my $file (@{$args{in}})
    {
        if ( $args{do_bcf} )
        {
            cmd("$$opts{bin}/bcftools view --no-version -Ob $$opts{tmp}/$file.vcf.gz > $$opts{tmp}/$file.bcf");
            cmd("$$opts{bin}/bcftools index -f $$opts{tmp}/$file.bcf");
            $files .= " $$opts{tmp}/$file.bcf";
        }
        else
        {
            bgzip_tabix_vcf($opts,$file);
            $files .= " $$opts{tmp}/$file.vcf.gz";
        }
    }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools concat --no-version $args{args} $files");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools concat -Ob $args{args} $files | $$opts{bin}/bcftools view | grep -v ^##bcftools_");
}
sub test_vcf_reheader
{
    my ($opts,%args) = @_;
    cmd("$$opts{bin}/bcftools view --no-version -Ob $$opts{path}/$args{in}.vcf > $$opts{tmp}/$args{in}.bcf");
    cmd("$$opts{bin}/bcftools view --no-version -Oz $$opts{path}/$args{in}.vcf > $$opts{tmp}/$args{in}.vcf.gz");

    my $arg;
    if ( exists($args{args}) )
    {
        $arg = $args{args};
        $arg =~ s/{PATH}/$$opts{path}/g;
    }
    elsif ( exists($args{header}) )
    {
        $arg = "-h $$opts{path}/$args{header}";
    }
    else
    {
        $arg = "-s $$opts{path}/$args{samples}";
    }
    for my $file ("$$opts{path}/$args{in}.vcf","$$opts{tmp}/$args{in}.bcf","$$opts{tmp}/$args{in}.vcf.gz")
    {
        # bcf header lines can come in different order
        my %bcf_args = ();
        if ( $file=~/\.bcf$/ && -e "$$opts{path}/$args{out}.bcf" ) { %bcf_args = ( out=>"$args{out}.bcf" ); }
        test_cmd($opts,%args,%bcf_args,cmd=>"$$opts{bin}/bcftools reheader $arg $file | $$opts{bin}/bcftools view --no-version");
        test_cmd($opts,%args,%bcf_args,cmd=>"cat $file | $$opts{bin}/bcftools reheader $arg | $$opts{bin}/bcftools view --no-version") unless $args{nostdin};
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
    my $chain = $args{chain} ? "-c $$opts{tmp}/$args{chain}" : '';
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools consensus $$opts{tmp}/$args{in}.vcf.gz -f $$opts{path}/$args{fa} $args{args} $mask $chain");
}
sub test_vcf_consensus_chain
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    my $mask = $args{mask} ? "-m $$opts{path}/$args{mask}" : '';
    my $chain = $args{chain} ? "-c $$opts{tmp}/$args{chain}.new" : '';
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools consensus $$opts{tmp}/$args{in}.vcf.gz -f $$opts{path}/$args{fa} $args{args} $mask $chain >/dev/null; cat $$opts{tmp}/$args{chain}.new");
}

sub test_naive_concat
{
    my ($opts,%args) = @_;

    my $seed = srand();
    print STDERR "Random seed for test_naive_concat: $seed\n";

    my @files = ();
    my $exp   = '';
    for (my $n=0; $n<$args{nfiles}; $n++)
    {
        my $nhdr = 1 + int(rand($args{max_hdr_lines}));
        my $nbdy = int(rand($args{max_body_lines}));
        my $file = "$$opts{tmp}/$args{name}.$n";
        push @files,$file;

        open(my $fh,'>',"$file.vcf") or error("$file.vcf: $!");
        print $fh "##fileformat=VCFv4.0\n";
        print $fh "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
        print $fh "##contig=<ID=1,length=62435964>\n";
        for (my $i=0; $i<$nhdr; $i++)
        {
            my $x = rand;
            print $fh "##INFO=<ID=XX$i,Number=1,Type=Integer,Description=\"Test Tag $x\">\n";
        }
        print $fh join("\t",'#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')."\n";

        # let one of the files have no body
        if ( $n!=3 )
        {
            for (my $i=1; $i<=$nbdy; $i++)
            {
                my $x = int(rand(1000));
                my $line = join("\t",'1',$i,'.','A','C','.','.',"DP=$x")."\n";
                print $fh  $line;
                $exp .= $line;
            }
        }
        close($fh) or error("close failed: $file.vcf");
    }

    for my $file (@files)
    {
        cmd("$$opts{bin}/bcftools view -Ob -o $file.bcf $file.vcf");
        cmd("$$opts{bin}/bcftools view -Oz -o $file.vcf.gz $file.vcf");
    }

    my $bcfs = join('.bcf ',@files).'.bcf';
    test_cmd($opts,exp=>$exp,out=>"concat.naive.bcf.out",cmd=>"$$opts{bin}/bcftools concat --naive-force $bcfs | $$opts{bin}/bcftools view -H");

    my $vcfs = join('.vcf.gz ',@files).'.vcf.gz';
    test_cmd($opts,exp=>$exp,out=>"concat.naive.vcf.out",cmd=>"$$opts{bin}/bcftools concat --naive $vcfs | $$opts{bin}/bcftools view -H");
}

sub test_mpileup
{
    my ($opts,%args) = @_;

    if ($args{test_list})
    {
        # make a local copy, create bams, index the bams and the reference
        open(my $fh1,'>',"$$opts{tmp}/mpileup.bam.list") or error("$$opts{tmp}/mpileup.bam.list: $!");
        open(my $fh2,'>',"$$opts{tmp}/mpileup.cram.list") or error("$$opts{tmp}/mpileup.cram.list: $!");
        open(my $fh3,'>',"$$opts{tmp}/mpileup.bam.urllist") or error("$$opts{tmp}/mpileup.bam.urllist: $!");
        open(my $fh4,'>',"$$opts{tmp}/mpileup.cram.urllist") or error("$$opts{tmp}/mpileup.cram.urllist: $!");
        for my $file (@{$args{in}})
        {
            print $fh1 "$$opts{path}/mpileup/$file.bam\n";
            print $fh2 "$$opts{path}/mpileup/$file.cram\n";
            my $atmp = $^O =~ /^msys/ ? cygpath("$$opts{path}/mpileup") : abs_path("$$opts{path}/mpileup");
            unless ($atmp =~ /^\//) { $atmp = "/$atmp"; }
            print $fh3 "file://$atmp/$file.bam\n";
            print $fh4 "file://$atmp/$file.cram\n";
        }
        close($fh1);
        close($fh2);
        close($fh3);
        close($fh4);
    }
    my $ref = exists($args{ref}) ? $args{ref} : "mpileup.ref.fa";

    $args{args} =~ s/{PATH}/$$opts{path}/g;
    for my $fmt ('bam','cram')
    {
        my @files = ();
        for my $file (@{$args{in}}) { push @files, "$$opts{path}/mpileup/$file.$fmt"; }
        my $files = join(' ',@files);
        my $grep_hdr = "grep -v ^##bcftools | grep -v ^##reference";
        test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools mpileup $args{args} -f $$opts{path}/mpileup/$ref $files | $grep_hdr");
        test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools mpileup $args{args} -f $$opts{path}/mpileup/$ref -Ob $files | $$opts{bin}/bcftools view  | $grep_hdr");
        if ($args{test_list})
        {
            test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools mpileup $args{args} -f $$opts{path}/mpileup/$ref -b $$opts{tmp}/mpileup.$fmt.list --no-version | grep -v ^##reference");
            test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools mpileup $args{args} -f $$opts{path}/mpileup/$ref -Ob -b $$opts{tmp}/mpileup.$fmt.urllist | $$opts{bin}/bcftools view  | $grep_hdr");
        }
    }
}

sub test_csq
{
    my ($opts,%args) = @_;
    $args{cmd}  =~ s/{PATH}/$$opts{path}/g;
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools csq $args{cmd} $$opts{path}/$args{in}.vcf | $$opts{bin}/test/csq/sort-csq | $$opts{bin}/bcftools query -f'%POS\\t%REF\\t%ALT\\t%EXP\\n%POS\\t%REF\\t%ALT\\t%BCSQ\\n\\n'");
}
sub test_csq_real
{
    my ($opts,%args) = @_;

    my $dirname = "$$opts{path}/$args{in}";
    opendir(my $dh,$dirname) or error("opendir $dirname: $!");
    while (my $dir=readdir($dh))
    {
        if ( !($dir=~/^E/) or !-d "$dirname/$dir" ) { next; }
        my $gff = "$dirname/$dir/$dir.gff";
        my $ref = "$dirname/$dir/$dir.fa";
        opendir(my $tmp,"$dirname/$dir") or error("opendir: $dirname/$dir: $!");
        while (my $file=readdir($tmp))
        {
            my $cmd = undef;
            
            if ( $file=~/\.vcf$/ )
            {
                my $bname = $`;
                my $vcf   = "$dirname/$dir/$file";
                my $out   = "$args{in}/$dir/$bname.txt";
                my $outl  = "$args{in}/$dir/$bname.txt-l";

                my @nsmpl = `$$opts{bin}/bcftools query -l $vcf`;
                if ( !@nsmpl )
                {
                    $cmd = "| $$opts{bin}/test/csq/sort-csq | $$opts{bin}/bcftools query -f'%POS\\t%REF\\t%ALT\\t%EXP\\n%POS\\t%REF\\t%ALT\\t%BCSQ\\n\\n'";
                }
                else
                {
                    $cmd = "| $$opts{bin}/bcftools query -f'[%POS\\t%REF\\t%ALT\\t%TBCSQ\\n]\\n'";
                }
                if ( -e "$$opts{path}/$out")
                {
                    test_cmd($opts,%args,out=>$out,cmd=>"$$opts{bin}/bcftools csq -f $ref -g $gff $vcf $cmd");
                }
                if ( -e "$$opts{path}/$outl" )
                {
                    test_cmd($opts,%args,out=>$outl,cmd=>"$$opts{bin}/bcftools csq -l -f $ref -g $gff $vcf $cmd");
                }
                next;
            }
            if ( $file=~/\.cmd$/ )
            {
                chomp(my $cmd = (`cat $dirname/$dir/$file`)[0]);
                $cmd =~ s/{bin}/$$opts{bin}/g;
                test_cmd($opts,%args,out=>"$args{in}/$dir/$file.out",cmd=>"cd $dirname/$dir && $cmd");
            }
        }
        closedir($tmp);
    }
    closedir($dh);
}
