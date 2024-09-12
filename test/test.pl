#!/usr/bin/env perl
#
#   Copyright (C) 2012-2024 Genome Research Ltd.
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
run_test(\&test_usage,$opts,cmd=>'bcftools');
run_test(\&test_tabix,$opts,in=>'merge.a',reg=>'2:3199812-3199812',out=>'tabix.2.3199812.out');
run_test(\&test_tabix,$opts,in=>'merge.a',reg=>'1:3000151-3000151',out=>'tabix.1.3000151.out');
run_test(\&test_index,$opts,in=>'large_chrom_csi_limit',reg=>'chr20:1-2147483647',out=>'large_chrom_csi_limit.20.1.2147483647.out'); # 2147483647 (1<<31-1) is the current chrom limit for csi. bcf conversion and indexing fail above this
run_test(\&test_index,$opts,in=>'large_chrom_csi_limit',reg=>'chr20',out=>'large_chrom.20.1.2147483647.out'); # this fails until bug resolved
run_test(\&test_vcf_idxstats,$opts,in=>'idx',args=>'-s',out=>'idx.out');
run_test(\&test_vcf_idxstats,$opts,in=>'idx',args=>'-n',out=>'idx_count.out');
run_test(\&test_vcf_idxstats,$opts,in=>'empty',args=>'-s',out=>'empty.idx.out');
run_test(\&test_vcf_idxstats,$opts,in=>'empty',args=>'-n',out=>'empty.idx_count.out');
run_test(\&test_vcf_check,$opts,in=>'check',out=>'check.chk');
run_test(\&test_vcf_check_merge,$opts,in=>'check',out=>'check_merge.chk');
run_test(\&test_vcf_stats,$opts,in=>['stats.a','stats.b'],out=>'stats.chk',args=>'-s -');
run_test(\&test_vcf_stats,$opts,in=>['stats.a','stats.b'],out=>'stats.B.chk',args=>'-s B');
run_test(\&test_vcf_stats,$opts,in=>['stats.counts'],out=>'stats.counts.chk',args=>'-s -');
run_test(\&test_vcf_stats,$opts,in=>['stats.counts'],out=>'stats.counts.2.chk',args=>q[-s - -i 'type="snp"']);
run_test(\&test_vcf_stats,$opts,in=>['stats.vaf'],out=>'stats.vaf.1.chk',args=>q[-s -]);
run_test(\&test_vcf_isec,$opts,in=>['isec.a','isec.b'],out=>'isec.ab.out',args=>'-n =2');
run_test(\&test_vcf_isec,$opts,in=>['isec.a','isec.b'],out=>'isec.ab.flt.out',args=>'-n =2 -i"STRLEN(REF)==2"');
run_test(\&test_vcf_isec,$opts,in=>['isec.a','isec.b'],out=>'isec.ab.both.out',args=>'-n =2 -c both');
run_test(\&test_vcf_isec,$opts,in=>['isec.a','isec.b'],out=>'isec.ab.any.out',args=>'-n =2 -c any');
run_test(\&test_vcf_isec,$opts,in=>['isec.a','isec.b'],out=>'isec.ab.C.out',args=>'-C -c any');
run_test(\&test_vcf_isec2,$opts,vcf_in=>['isec.a'],tab_in=>'isec',out=>'isec.tab.out',args=>'');
run_test(\&test_vcf_isec,$opts,in=>['isec-miss.1.1','isec-miss.1.2','isec-miss.1.3'],out=>'isec-miss.1.1.out',args=>'-n +1 -r 20:100,20:140,12:55,20:140,20:100');
run_test(\&test_vcf_isec,$opts,in=>['isec-miss.1.1','isec-miss.1.2','isec-miss.1.3'],out=>'isec-miss.1.1.out',args=>'-R {PATH}/isec-miss.1.regs.txt -n +1');
run_test(\&test_vcf_isec,$opts,in=>['isec-miss.2.1','isec-miss.2.2','isec-miss.2.3'],out=>'isec-miss.2.1.out',args=>'-n +1 -r 20:100,20:140,12:55,20:140,20:100');
run_test(\&test_vcf_isec,$opts,in=>['isec-miss.2.1','isec-miss.2.2','isec-miss.2.3'],out=>'isec-miss.2.1.out',args=>'-R {PATH}/isec-miss.1.regs.txt -n +1');
run_test(\&test_vcf_merge,$opts,in=>['merge.11.a','merge.11.b'],out=>'merge.11.1.out',args=>'');
run_test(\&test_vcf_merge,$opts,in=>['merge.join.a','merge.join.b'],out=>'merge.join.1.out',args=>'-i AF:join');
run_test(\&test_vcf_merge,$opts,in=>['merge.LPL.a'],out=>'merge.LPL.0.out',args=>'--force-single');
run_test(\&test_vcf_merge,$opts,in=>['merge.LPL.a','merge.LPL.b','merge.LPL.c'],out=>'merge.LPL.1.out',args=>'--force-samples');
run_test(\&test_vcf_merge,$opts,in=>['merge.LPL.a','merge.LPL.b','merge.LPL.c'],out=>'merge.LPL.2.out',args=>'--force-samples -L 1');
run_test(\&test_vcf_merge,$opts,in=>['merge.LPL.a','merge.LPL.b','merge.LPL.c'],out=>'merge.LPL.3.out',args=>'--force-samples -L 2');
run_test(\&test_vcf_merge,$opts,in=>['merge.LPL.a','merge.LPL.b','merge.LPL.c'],out=>'merge.LPL.4.out',args=>'--force-samples -L 3');
run_test(\&test_vcf_merge,$opts,in=>['merge.LPL.a','merge.LPL.b','merge.LPL.c'],out=>'merge.LPL.5.out',args=>'--force-samples -L 4');
run_test(\&test_vcf_merge,$opts,in=>['merge.LPL.a','merge.LPL.b','merge.LPL.c'],out=>'merge.LPL.6.out',args=>'--force-samples -L 5');
run_test(\&test_vcf_merge,$opts,in=>['merge.a','merge.b','merge.c'],out=>'merge.abc.out',args=>'--force-samples');
run_test(\&test_vcf_merge,$opts,in=>['merge.a','merge.b','merge.c'],out=>'merge.abc.2.out',args=>'--force-samples -Fx');
run_test(\&test_vcf_merge,$opts,in=>['merge.a','merge.b','merge.c'],out=>'merge.abc.3.out',args=>'--force-samples -0');
run_test(\&test_vcf_merge,$opts,in=>['merge.2.a','merge.2.b'],out=>'merge.2.none.out',args=>'--force-samples -m none');
run_test(\&test_vcf_merge,$opts,in=>['merge.2.a','merge.2.b'],out=>'merge.2.both.out',args=>'--force-samples -m both');
run_test(\&test_vcf_merge,$opts,in=>['merge.2.a','merge.2.b'],out=>'merge.2.all.out',args=>'--force-samples -m all');
run_test(\&test_vcf_merge,$opts,in=>['merge.3.a','merge.3.b'],out=>'merge.3.out',args=>'--force-samples -i TR:sum,TA:sum,TG:sum');
run_test(\&test_vcf_merge,$opts,in=>['merge.4.a','merge.4.b'],out=>'merge.4.out',args=>'--force-samples -m id');
run_test(\&test_vcf_merge,$opts,in=>['gvcf.merge.1','gvcf.merge.2','gvcf.merge.3'],out=>'gvcf.merge.1.out',args=>'--gvcf -');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.2.a','merge.gvcf.2.b','merge.gvcf.2.c'],out=>'merge.gvcf.2.out',args=>'--gvcf -',types=>['vcf']);
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.2.a','merge.gvcf.2.b','merge.gvcf.2.c'],out=>'merge.gvcf.2.1.out',args=>'--gvcf - -m both,*',types=>['vcf']);
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.2.a','merge.gvcf.2.b','merge.gvcf.2.c'],out=>'merge.gvcf.2.2.out',args=>'--gvcf - -m both,**',types=>['vcf']);
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.3.a','merge.gvcf.3.b'],out=>'merge.gvcf.3.out',args=>'--gvcf - -i SRC:join',types=>['vcf']);
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.4.a','merge.gvcf.4.b'],out=>'merge.gvcf.4.out',args=>'--gvcf -');
run_test(\&test_vcf_merge,$opts,in=>['merge.5.a','merge.5.b'],out=>'merge.5.out');
run_test(\&test_vcf_merge,$opts,in=>['merge.6.a','merge.6.b'],out=>'merge.6.out');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.7.a','merge.gvcf.7.b'],out=>'merge.gvcf.7.out',args=>'--gvcf -');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.8.a','merge.gvcf.8.b'],out=>'merge.gvcf.8.out',args=>'--gvcf -');
run_test(\&test_vcf_merge,$opts,in=>['merge.7.a','merge.7.b'],out=>'merge.9.out',args=>'--force-samples');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.9a','merge.gvcf.9b','merge.gvcf.9c','merge.gvcf.9d'],out=>'merge.gvcf.9.1.out',args=>'--gvcf {PATH}/gvcf.fa');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.9a','merge.gvcf.9b','merge.gvcf.9c','merge.gvcf.9d'],out=>'merge.gvcf.9.2.out',args=>'--gvcf {PATH}/gvcf.fa -r 22:21-23');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.9a','merge.gvcf.9b','merge.gvcf.9c','merge.gvcf.9d','merge.gvcf.9e'],out=>'merge.gvcf.9.3.out',args=>'--gvcf {PATH}/gvcf.fa');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.9a','merge.gvcf.9b','merge.gvcf.9c','merge.gvcf.9d','merge.gvcf.9e'],out=>'merge.gvcf.9.4.out',args=>'--gvcf {PATH}/gvcf.fa -r 22:21-23');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.10.a','merge.gvcf.10.b'],out=>'merge.gvcf.10.1.out',args=>'');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.10.a','merge.gvcf.10.b'],out=>'merge.gvcf.10.2.out',args=>'-m none');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.10.a','merge.gvcf.10.b'],out=>'merge.gvcf.10.3.out',args=>'-g {PATH}/merge.gvcf.10.fa');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.10.a','merge.gvcf.10.b'],out=>'merge.gvcf.10.4.out',args=>'-g {PATH}/merge.gvcf.10.fa -m none');
run_test(\&test_vcf_merge,$opts,in=>['merge.noidx.a','merge.noidx.b','merge.noidx.c'],out=>'merge.noidx.abc.out',args=>'');
run_test(\&test_vcf_merge,$opts,in=>['merge.noidx.a','merge.noidx.b','merge.noidx.c'],out=>'merge.noidx.abc.out',args=>'--no-index',noidx=>1);
run_test(\&test_vcf_merge,$opts,in=>['merge.8.a','merge.8.b'],out=>'merge.8.out',args=>'');
run_test(\&test_vcf_merge,$opts,in=>['merge.8.a','merge.8.b'],out=>'merge.8.out',args=>'-i AN:sum,AC:sum');
run_test(\&test_vcf_merge,$opts,in=>['merge.9.a','merge.9.b'],out=>'merge.9.1.out',args=>'');
run_test(\&test_vcf_merge,$opts,in=>['merge.9.a','merge.9.b'],out=>'merge.9.2.out',args=>'-i AN:sum,AC:sum');
run_test(\&test_vcf_merge,$opts,in=>['merge.10.a','merge.10.b'],out=>'merge.10.1.out',args=>'-m none');
run_test(\&test_vcf_merge,$opts,in=>['merge.10.a','merge.10.b'],out=>'merge.10.2.out',args=>'-m both');
run_test(\&test_vcf_merge,$opts,in=>['merge.10.a','merge.10.b'],out=>'merge.10.3.out',args=>'-m snp-ins-del');
run_test(\&test_vcf_merge,$opts,in=>['merge.mrules.1.a','merge.mrules.1.b'],out=>'merge.mrules.1.1.out',args=>'--gvcf -');
run_test(\&test_vcf_merge,$opts,in=>['merge.mrules.1.a','merge.mrules.1.b'],out=>'merge.mrules.1.2.out',args=>'--gvcf - -M AD:.,PL:.');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.5.a','merge.gvcf.5.b'],out=>'merge.gvcf.5.1.out',args=>'--gvcf -');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.5.a','merge.gvcf.5.b'],out=>'merge.gvcf.5.1.out',args=>'--gvcf - --merge none');
run_test(\&test_vcf_merge,$opts,in=>['merge.gvcf.11.a','merge.gvcf.11.b','merge.gvcf.11.c'],out=>'merge.gvcf.11.1.out',args=>'--gvcf -');
# run_test(\&test_vcf_merge_big,$opts,in=>'merge_big.1',out=>'merge_big.1.1',nsmpl=>79000,nfiles=>79,nalts=>486,args=>'');   # commented out for speed
run_test(\&test_vcf_query,$opts,in=>'query.string',out=>'query.string.1.out',args=>q[-f '%CHROM\\t%POS\\t%CLNREVSTAT\\n' -i'CLNREVSTAT="criteria_provided,_conflicting_interpretations"']);
run_test(\&test_vcf_query,$opts,in=>'query.string',out=>'query.string.1.out',args=>q[-f '%CHROM\\t%POS\\t%CLNREVSTAT\\n' -i'CLNREVSTAT="criteria_provided" || CLNREVSTAT="_conflicting_interpretations"']);
run_test(\&test_vcf_query,$opts,in=>'query.string',out=>'query.string.2.out',args=>q[-f '%CHROM\\t%POS\\t%CLNREVSTAT\\n' -i'CLNREVSTAT="criteria_provided" && CLNREVSTAT="_conflicting_interpretations"']);
run_test(\&test_vcf_query,$opts,in=>'query.string.2',out=>'query.string.2.1.out',args=>q[-f '%CHROM\\t%POS\\t%INFO/STR\\n' -i'INFO/STR=@{PATH}/query.string.2.1.txt']);
run_test(\&test_vcf_query,$opts,in=>'query.string.2',out=>'query.string.2.2.out',args=>q[-f '%CHROM\\t%POS[\\t%STR]\\n' -i'FMT/STR=@{PATH}/query.string.2.2.txt']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.out',args=>q[-f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%DP4\\t%AN[\\t%GT\\t%TGT]\\n']);
run_test(\&test_vcf_query,$opts,in=>'query.variantkey',out=>'query.variantkey.hex.out',args=>q[-f '%RSX\\t%VKX\\n']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.2.out',args=>q[-f'%XRI\\n' -i'XRI[*]>1111']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.3.out',args=>q[-f'%XRF\\n' -i'XRF[*]=2e6']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[5]="PQR"']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[*]="GHI"']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[2]~"H"']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[3]!~"H" && XGS!="."']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.4.out',args=>q[-f'%XGS\\n' -i'XGS[*]~"H"']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.5.out',args=>q[-f'%POS %REF %ALT\\n' -i'REF~"C" && ALT[*]~"CT"']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.6.out',args=>q[-f'%POS %REF %ALT\\n' -i'N_ALT=2']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.7.out',args=>q[-f'%POS %AN\\n' -i'AN!=2*N_SAMPLES']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.8.out',args=>q[-f'%POS[ %GL]\\n' -i'min(abs(GL[*:0]))=10']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.9.out',args=>q[-f'%POS %CIGAR\\n' -i'strlen(CIGAR[*])=4']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.10.out',args=>q[-f'%POS[ %GT]\\n' -i'AC[0]=3']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.10.out',args=>q[-f'%POS[ %GT]\\n' -i'AF[0]=3/4']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.11.out',args=>q[-f'%POS[ %GT]\\n' -i'MAC[0]=1']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.11.out',args=>q[-f'%POS[ %GT]\\n' -i'MAF[0]=1/4']);
run_test(\&test_vcf_query,$opts,in=>'view.vectors',out=>'query.12.out',args=>q[-f'I8=%I8 I16=%I16 I32=%I32 IF=%IF IA8=%IA8 IA16=%IA16 IA32=%IA32 IAF=%IAF IA8=%IA8{1} IA16=%IA16{1} IA32=%IA32{1} IAF=%IAF{1} [ %F8:%F16:%F32:%FF]\\n']);
run_test(\&test_vcf_query,$opts,in=>'query.filter',out=>'query.13.out',args=>q[-f'[%POS %SAMPLE %GT\\n]' -i'GT ="1"']);
run_test(\&test_vcf_query,$opts,in=>'query.filter',out=>'query.14.out',args=>q[-f'[%POS %SAMPLE %GT\\n]' -i'GT!="1"']);
run_test(\&test_vcf_query,$opts,in=>'query.filter',out=>'query.15.out',args=>q[-f'[%POS %SAMPLE %GT\\n]' -e'GT ="1"']);
run_test(\&test_vcf_query,$opts,in=>'query.filter',out=>'query.16.out',args=>q[-f'[%POS %SAMPLE %GT\\n]' -e'GT!="1"']);
run_test(\&test_vcf_query,$opts,in=>'query.2',out=>'query.17.out',args=>q[-f'%XX_A %XX.A %XX.A0 %xx.a0\\n']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.18.out',args=>q[-i'IINT="."'  -f'%POS %IINT\\n']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.19.out',args=>q[-i'IINT!="."' -f'%POS %IINT\\n']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.18.out',args=>q[-e'IINT!="."' -f'%POS %IINT\\n']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.19.out',args=>q[-e'IINT="."'  -f'%POS %IINT\\n']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.20.out',args=>q[-i'IFLT="."'  -f'%POS %IFLT\\n']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.21.out',args=>q[-i'IFLT!="."' -f'%POS %IFLT\\n']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.20.out',args=>q[-e'IFLT!="."' -f'%POS %IFLT\\n']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.21.out',args=>q[-e'IFLT="."'  -f'%POS %IFLT\\n']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.22.out',args=>q[-i'ISTR="."'  -f'%POS %ISTR\\n']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.23.out',args=>q[-i'ISTR!="."' -f'%POS %ISTR\\n']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.23.out',args=>q[-e'ISTR="."'  -f'%POS %ISTR\\n']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.22.out',args=>q[-e'ISTR!="."' -f'%POS %ISTR\\n']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.24.out',args=>q[-i'FILTER="q11"' -f'%POS %ISTR\\n']);
run_test(\&test_vcf_query,$opts,in=>'filter.11',out=>'query.76.out',args=>q[-i'FILTER="A"' -f'%POS %FILTER\\n']);
run_test(\&test_vcf_query,$opts,in=>'filter.11',out=>'query.77.out',args=>q[-i'FILTER!="A"' -f'%POS %FILTER\\n']);
run_test(\&test_vcf_query,$opts,in=>'filter.11',out=>'query.78.out',args=>q[-i'FILTER~"A"' -f'%POS %FILTER\\n']);
run_test(\&test_vcf_query,$opts,in=>'filter.11',out=>'query.79.out',args=>q[-i'FILTER!~"A"' -f'%POS %FILTER\\n']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.25.out',args=>q[-f'%LINE']);
run_test(\&test_vcf_query,$opts,in=>'query.filter-type',out=>'query.26.out',args=>q[-f'%POS\\t%REF\\t%ALT\\n' -i'type="snp"']);
run_test(\&test_vcf_query,$opts,in=>'query.filter-type',out=>'query.27.out',args=>q[-f'%POS\\t%REF\\t%ALT\\n' -i'type~"snp"']);
run_test(\&test_vcf_query,$opts,in=>'query.filter-type',out=>'query.28.out',args=>q[-f'%POS\\t%REF\\t%ALT\\n' -i'type!="snp"']);
run_test(\&test_vcf_query,$opts,in=>'query.filter-type',out=>'query.29.out',args=>q[-f'%POS\\t%REF\\t%ALT\\n' -i'type!~"snp"']);
run_test(\&test_vcf_query,$opts,in=>'query.filter-type',out=>'query.67.out',args=>q[-f'%POS\\t%REF\\t%ALT\\n' -i'INFO/TYPE="xxx"']);
run_test(\&test_vcf_query,$opts,in=>'filter-missing-floats',out=>'query.30.out',args=>q[-f'%POS\\t%A_AF\\t%B_AF\\t%C_AF\\n' -i'A_AF>=0.0001 || B_AF >= 0.0001 || C_AF >= 0.0001']);
run_test(\&test_vcf_query,$opts,in=>'filter-missing-floats',out=>'query.31.out',args=>q[-f'%POS\\t%A_AF\\t%B_AF\\t%C_AF\\n' -e'A_AF>=0.0001 || B_AF >= 0.0001 || C_AF >= 0.0001']);
run_test(\&test_vcf_query,$opts,in=>'missing',out=>'query.32.out',args=>q[-i'FMT/FINT!="."' -f'[\t%FINT]\\n']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.2',out=>'query.33.out',args=>q[-f'[%GT]\\n' -i'GT~"0/[1-9]" || GT~"[1-9]/0"']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.34.out',args=>q[-f'[%POS %SAMPLE %FGS\\n]\\n' -i'FGS[*:1-2,4]="EE"']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.34.out',args=>q[-f'[%POS %SAMPLE %FGS\\n]\\n' -i'"EE"=FGS[*:1-2,4]']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.35.out',args=>q[-f'[%POS %SAMPLE %FGI\\n]\\n' -i'FGI[*:1-2,5]=6']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.35.out',args=>q[-f'[%POS %SAMPLE %FGI\\n]\\n' -i'6=FGI[*:1-2,5]']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.36.out',args=>q[-f'[%POS %SAMPLE %FGF\\n]\\n' -i'FGF[*:1-3,4]=5e-5']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.37.out',args=>q[-f'[%POS %SAMPLE %FGS\\n]\\n' -i'FGS[*:4]!="EE"']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.38.out',args=>q[-f'[%POS %SAMPLE %FGS\\n]\\n' -i'FGS[*:4]="."']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.39.out',args=>q[-f'[%POS %SAMPLE %FGS\\n]\\n' -i'FGS[*:4]!="."']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.40.out',args=>q[-f'[%POS %SAMPLE %FGI\\n]\\n' -i'FGI[*:1-2,5]="."']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.41.out',args=>q[-f'[%POS %SAMPLE %FGI\\n]\\n' -i'FGI[*:1-2,5]!="."']);
run_test(\&test_vcf_query,$opts,in=>'filter.5',out=>'query.42.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP=1 | FMT/DP=2']);
run_test(\&test_vcf_query,$opts,in=>'filter.5',out=>'query.43.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP=1 | FMT/DP="."']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.44.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP=19']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.45.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP=19 |  FMT/DP="."']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.46.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP=19 || FMT/DP="."']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.47.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP[0]=19']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.48.out',args=>q[-f'[%POS %SAMPLE %DP\\n]\\n' -i'FMT/DP[0]!=19']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.49.out',args=>q[-f'[%POS %SAMPLE %DP %GQ\\n]\\n' -i'FMT/DP=19 &  FMT/GQ=589']);
run_test(\&test_vcf_query,$opts,in=>'view.filter',out=>'query.50.out',args=>q[-f'[%POS %SAMPLE %DP %GQ\\n]\\n' -i'FMT/DP=19 && FMT/GQ=1']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.3',out=>'query.51.out',args=>q[-f'[\\t%GT\\n]\\n' -i'GT~"1" && GT~"2"']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.3',out=>'query.52.out',args=>q[-f'[\\t%GT\\n]\\n' -i'GT~"1" &  GT~"2"']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.3',out=>'query.53.out',args=>q[-f'%POS[\\t%GT]\\n' -i'COUNT(GT="het")=1']);
run_test(\&test_vcf_query,$opts,in=>'filter.5',out=>'query.54.out',args=>q[-f'[%POS  %SAMPLE  %AD\\n]\\n' -i'AD[:0]+AD[:1] > 12']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.4',out=>'query.55.out',args=>q[-f'%POS\\t%REF\\t%ALT[\\t%GT]\\n' -e'TYPE!="snp" || ALT="*"']);
run_test(\&test_vcf_query,$opts,in=>'view',out=>'query.56.out',args=>q[-f'%ID\\n' -i 'ID=@].$$opts{path}.q[/query.56.out']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.5',out=>'query.57.out',args=>q[-f'[%POS\\t%SAMPLE\\t%GT\\t%AD\\n]' -i'GT="het" & binom(FMT/AD)>0.01']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.5',out=>'query.58.out',args=>q[-f'[%POS\\t%SAMPLE\\t%GT\\t%AD\\n]' -i'GT="het" & binom(FMT/AD[:0],FMT/AD[:1])>0.01']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.5',out=>'query.59.out',args=>q[-f'%POS\\t%AD\\n' -i'binom(INFO/AD[0],INFO/AD[1])>0.01']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.5',out=>'query.92.out',args=>q[-f'[%POS\\t%GT\\t%SAMPLE\\t%AD\\n]' -i'FMT/AD[GT]==10']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.5',out=>'query.93.out',args=>q[-f'[%POS\\t%GT\\t%SAMPLE\\t%AD\\n]' -i'sSUM(FMT/AD[GT])==210']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.5',out=>'query.94.out',args=>q[-f'[%POS\\t%GT\\t%SAMPLE\\t%AD\\n]' -i'FMT/AD[0:GT]==30']);
run_test(\&test_vcf_query,$opts,in=>'query.gtidx.1',out=>'query.gtidx.1.out',args=>q[-f'[%GT %TC %TS %TI %TF\\n]' -i'FMT/TF[:GT]==1 && FMT/TF[:GT]==2']);
run_test(\&test_vcf_query,$opts,in=>'query.gtidx.1',out=>'query.gtidx.1.out',args=>q[-f'[%GT %TC %TS %TI %TF\\n]' -i'FMT/TI[:GT]==1 && FMT/TF[:GT]==2']);
run_test(\&test_vcf_query,$opts,in=>'query.gtidx.1',out=>'query.gtidx.1.out',args=>q[-f'[%GT %TC %TS %TI %TF\\n]' -i'FMT/TS[:GT]=="BB" && FMT/TS[:GT]=="CC"']);
run_test(\&test_vcf_query,$opts,in=>'query.gtidx.1',out=>'query.gtidx.1.out',args=>q[-f'[%GT %TC %TS %TI %TF\\n]' -i'FMT/TC[:GT]=="B"  && FMT/TC[:GT]=="C"']);
run_test(\&test_vcf_query,$opts,in=>'query.gtidx.1',out=>'query.gtidx.2.out',args=>q[-f'[%GT %TI\\n]' -i'FMT/TI[:GT]=="."']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.60.out',args=>q[-f'%CHROM %POS\\n' -i'CHROM="4"']);
run_test(\&test_vcf_query,$opts,in=>'query.negative',out=>'query.61.out',args=>q[-f'%POS\\t%TAG1\\n' -i'(TAG1>=-129 && TAG1<=-120) || (TAG1>=-32769 && TAG1<=-32760)']);
run_test(\&test_vcf_query,$opts,in=>'query.negative',out=>'query.61.out',args=>q[-f'%POS\\t%TAGV1\\n' -i'(TAGV1>=-129 && TAGV1<=-120) || (TAGV1>=-32769 && TAGV1<=-32760)']);
run_test(\&test_vcf_query,$opts,in=>'query.negative',out=>'query.62.out',args=>q[-f'%POS\\t%TAG2\\n' -i'(TAG2>=-129 && TAG2<=-120) || (TAG2>=-32769 && TAG2<=-32760)']);
run_test(\&test_vcf_query,$opts,in=>'query.negative',out=>'query.62.out',args=>q[-f'%POS\\t%TAGV2\\n' -i'(TAGV2>=-129 && TAGV2<=-120) || (TAGV2>=-32769 && TAGV2<=-32760)']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.63.out',args=>q[-f'[%POS\\t%SAMPLE\\t%GQ\\n]' -i'N_PASS(GQ<20)==1']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.11',out=>'query.80.out',args=>q[-f'[%POS\\t%SAMPLE\\t%GT\\n]' -i'N_PASS(GT="alt")==1']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.11',out=>'query.81.out',args=>q[-f'[%POS\\t%SAMPLE\\t%GT\\n]' -i'N_PASS(GT="alt")!=1']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.11',out=>'query.81.out',args=>q[-f'[%POS\\t%SAMPLE\\t%GT\\n]' -e'N_PASS(GT="alt")==1']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.11',out=>'query.80.out',args=>q[-f'[%POS\\t%SAMPLE\\t%GT\\n]' -e'N_PASS(GT="alt")!=1']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.64.out',args=>q[-f'%CHROM\\t%POS\\t%INFO\\t%FORMAT\\n' -s D,C]);
run_test(\&test_vcf_query,$opts,in=>'query.pbinom.1',out=>'query.65.out',args=>q[-f'[%POS %SAMPLE %GT %AD %PBINOM(AD)\\n]' -i'phred(binom(FMT/AD))>=0']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.6',out=>'query.66.out',args=>q[-f'%POS\\n' -i'POS==16777217 || POS==33554432 || POS=118673904']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'sum(II)==6']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'sum(FORMAT/FI)==7']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'avg(II)==2']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'avg(FORMAT/FI)==1.75']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'mean(II)==2']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'mean(FORMAT/FI)==1.75']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'median(II)==2']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.7',out=>'query.68.out',args=>q[-f'%POS\\t%II[\\t%FI]\\n' -i'median(FORMAT/FI)==1.5']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.8',out=>'query.69.out',args=>q[-f'%POS\\t%REF\\t%ALT\\t%ILEN\\n' -i'%ILEN==1']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.8',out=>'query.70.out',args=>q[-f'%POS\\t%REF\\t%ALT\\t%ILEN\\n' -i'ILEN==1']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.9',out=>'query.71.out',args=>q[-f'[%POS  %SAMPLE  %AD\\n]' -i'FMT/AD[:0] < FMT/AD[:1]']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.9',out=>'query.72.out',args=>q[-f'[%POS  %SAMPLE  %AD\\n]' -i'FMT/AD[:0] > FMT/AD[:1]']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.13',out=>'query.84.out',args=>q[-f'[ %AD\\n]' -i'AD[:1] / sum(AD[*]) > 0.5']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.10',out=>'query.73.out',args=>q[-f'%POS  %NUM_TAG\\n' -i'COUNT(INFO/NUM_TAG)=2']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.10',out=>'query.74.out',args=>q[-f'%POS  %STR_TAG\\n' -i'COUNT(INFO/STR_TAG)=2']);
run_test(\&test_vcf_query,$opts,in=>'query',out=>'query.75.out',args=>q[-f '%CHROM:%POS\\t%N_PASS(GT="alt" & GQ>110)\\t[\\t%GT]\\t[\\t%GQ]\n']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.12',out=>'query.82.out',args=>q[-f '%CHROM:%POS[\\t%SAMPLE=%GT]\\n' -e 'GT="mis"' -s 1,3,0]);
run_test(\&test_vcf_query,$opts,in=>'query.filter.12',out=>'query.83.out',args=>q[-f '%CHROM:%POS[\\t%SAMPLE=%GT]\\n' -e 'GT="mis"' -s 0,1,3]);
run_test(\&test_vcf_query,$opts,in=>'query.smpl',out=>'query.smpl.1.out',args=>q[-f "[%SAMPLE %GT\n]" -S {PATH}/query.smpl.11.txt]);
run_test(\&test_vcf_query,$opts,in=>'query.smpl',out=>'query.smpl.2.out',args=>q[-l -S {PATH}/query.smpl.11.txt]);
run_test(\&test_vcf_query,$opts,in=>'query.smpl',out=>'query.smpl.3.out',args=>q[-f "[%SAMPLE %GT\n]" -S ^{PATH}/query.smpl.11.txt]);
run_test(\&test_vcf_query,$opts,in=>'query.smpl',out=>'query.smpl.4.out',args=>q[-l -S ^{PATH}/query.smpl.11.txt]);
run_test(\&test_vcf_query,$opts,in=>'query.smpl',out=>'query.smpl.5.out',args=>q[-f "[%SAMPLE %GT\n]" -S {PATH}/query.smpl.txt]);
run_test(\&test_vcf_query,$opts,in=>'query.smpl',out=>'query.smpl.6.out',args=>q[-l -S {PATH}/query.smpl.txt]);
run_test(\&test_vcf_query,$opts,in=>'query.filter.id',out=>'query.filter.id.1.out',args=>q[-f'%ID\\n' -i'ID~"s12"']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.id',out=>'query.filter.id.2.out',args=>q[-f'%ID\\n' -i'ID="rs123"']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.id',out=>'query.filter.id.3.out',args=>q[-f'%ID\\n' -i'ID="abc"']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.id',out=>'query.filter.id.3.out',args=>q[-f'%ID\\n' -i'ID=@].$$opts{path}.q[/query.filter.id.3.txt']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.id',out=>'query.filter.id.4.out',args=>q[-f'%ID\\n' -i'ID!="abc"']);
run_test(\&test_vcf_query,$opts,in=>'query.filter.id',out=>'query.filter.id.4.out',args=>q[-f'%ID\\n' -i'ID!=@].$$opts{path}.q[/query.filter.id.3.txt']);
run_test(\&test_vcf_query,$opts,in=>'filter.12',out=>'query.85.out',args=>q[-i'FILTER="A"' -f'%FILTER\\n']);
run_test(\&test_vcf_query,$opts,in=>'filter.12',out=>'query.86.out',args=>q[-i'FILTER~"A"' -f'%FILTER\\n']);
run_test(\&test_vcf_query,$opts,in=>'filter.12',out=>'query.87.out',args=>q[-i'FILTER="A;B"' -f'%FILTER\\n']);
run_test(\&test_vcf_query,$opts,in=>'filter.12',out=>'query.88.out',args=>q[-i'FILTER!="A;B"' -f'%FILTER\\n']);
run_test(\&test_vcf_query,$opts,in=>'filter.12',out=>'query.89.out',args=>q[-i'FILTER~"A;B"' -f'%FILTER\\n']);
run_test(\&test_vcf_query,$opts,in=>'filter.12',out=>'query.90.out',args=>q[-i'FILTER!~"A;B"' -f'%FILTER\\n']);
run_test(\&test_vcf_query,$opts,in=>'filter.10',out=>'query.91.out',args=>q[-i'DP%10==2' -f'[ %DP]\\n']);
run_test(\&test_vcf_query,$opts,in=>'filter.13',out=>'query.99.out',args=>q[-i'REF="N"' -f'%CHROM %POS %REF %ALT %QUAL\\n']);
run_test(\&test_vcf_query,$opts,in=>'query.header',out=>'query.95.out',args=>q[-H -f'[%CHROM %POS  %SAMPLE %DP %GT\\n]']);
run_test(\&test_vcf_query,$opts,in=>'query.header',out=>'query.96.out',args=>q[-H -f'[%CHROM %POS  %SAMPLE %DP %GT\\t]\\n']);
run_test(\&test_vcf_query,$opts,in=>'query.header',out=>'query.97.out',args=>q[-H -f'%CHROM %POS[ %SAMPLE %DP %GT]\\n']);
run_test(\&test_vcf_query,$opts,in=>'query.header',out=>'query.97.out',args=>q[-H -f'%CHROM %POS[ %SAMPLE %DP %GT]']);
run_test(\&test_vcf_query,$opts,in=>'query.header',out=>'query.98.out',args=>q[-H -f'%CHROM %POS[ %SAMPLE][ %DP][ %GT]\\n']);
run_test(\&test_vcf_query,$opts,in=>'query.header',out=>'query.98.out',args=>q[-H -f'%CHROM %POS[ %SAMPLE][ %DP][ %GT]']);
run_test(\&test_vcf_query,$opts,in=>'query.header',out=>'query.98.2.out',args=>q[-HH -f'%CHROM %POS[ %SAMPLE][ %DP][ %GT]']);
run_test(\&test_vcf_query,$opts,in=>'query.filter-or',out=>'query.filter-or.1.out',args=>q[-f'[%SAMPLE %DP\\n]' -i'DP=1 || DP=2']);
run_test(\&test_vcf_query,$opts,in=>'query.filter-or',out=>'query.filter-or.2.out',args=>q[-f'[%SAMPLE %DP\\n]' -i'DP=1 |  DP=2']);
run_test(\&test_vcf_norm,$opts,in=>'norm.sort',out=>'norm.sort.1.out',args=>'-m -');
run_test(\&test_vcf_norm,$opts,in=>'norm.sort',out=>'norm.sort.2.out',args=>'-m - -S lex');
run_test(\&test_vcf_norm,$opts,in=>'norm.join-missing-ploidy',out=>'norm.join-missing-ploidy.1.out',args=>'-m +both');
run_test(\&test_vcf_norm,$opts,in=>'norm.split.5',out=>'norm.split.5.1.out',args=>'-m - --multi-overlaps .');
run_test(\&test_vcf_norm,$opts,in=>'norm.symbolic.3',out=>'norm.symbolic.3.1.out',fai=>'norm.symbolic.3',args=>'');
run_test(\&test_vcf_norm,$opts,in=>'norm',out=>'norm.out',fai=>'norm',args=>'-cx');
run_test(\&test_vcf_norm,$opts,in=>'norm.split',out=>'norm.split.out',args=>'-m-');
run_test(\&test_vcf_norm,$opts,in=>'norm.split.2',out=>'norm.split.2.out',args=>'-m-');
run_test(\&test_vcf_norm,$opts,in=>'norm.split.3',out=>'norm.split.3.out',args=>'-m- --force');
run_test(\&test_vcf_norm,$opts,in=>'norm.split.4',out=>'norm.split.4.1.out',args=>'-m-');
run_test(\&test_vcf_norm,$opts,in=>'norm.split.4',out=>'norm.split.4.2.out',args=>'-m- --keep-sum AD');
run_test(\&test_vcf_norm,$opts,in=>'norm.split',fai=>'norm',out=>'norm.split.and.norm.out',args=>'-m-');
run_test(\&test_vcf_norm,$opts,in=>'norm.merge',out=>'norm.merge.out',args=>'-m+');
run_test(\&test_vcf_norm,$opts,in=>'norm.merge.2',out=>'norm.merge.2.out',args=>'-m+');
run_test(\&test_vcf_norm,$opts,in=>'norm.merge.3',out=>'norm.merge.3.out',args=>'-m+');
run_test(\&test_vcf_norm,$opts,in=>'norm.merge',out=>'norm.merge.strict.out',args=>'-m+ -s');
run_test(\&test_vcf_norm,$opts,in=>'norm.setref',out=>'norm.setref.out',args=>'-Nc s',fai=>'norm');
run_test(\&test_vcf_norm,$opts,in=>'norm.telomere',out=>'norm.telomere.out',fai=>'norm');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup',out=>'norm.rmdup.1.out',args=>'-d snps');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup',out=>'norm.rmdup.2.out',args=>'-d indels');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup',out=>'norm.rmdup.3.out',args=>'-d both');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup',out=>'norm.rmdup.4.out',args=>'-d all');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup',out=>'norm.rmdup.5.out',args=>'-d none');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup',out=>'norm.rmdup.5.out',args=>'-d exact');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup.2',out=>'norm.rmdup.2.1.out',args=>'-d none');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup.2',out=>'norm.rmdup.2.1.out',args=>'-d exact');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup.2',out=>'norm.rmdup.2.1.out',args=>'-d indels');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup.2',out=>'norm.rmdup.2.2.out',args=>'-d any');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup.2',out=>'norm.rmdup.2.2.out',args=>'-d both');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup.2',out=>'norm.rmdup.2.2.out',args=>'-d snps');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup.3',fai=>'norm.rmdup.3',out=>'norm.rmdup.3.1.out',args=>'-d exact');
run_test(\&test_vcf_norm,$opts,in=>'norm.rmdup.3',fai=>'norm.rmdup.3',out=>'norm.rmdup.3.2.out',args=>'-d all');
run_test(\&test_vcf_norm,$opts,in=>'norm.2',fai=>'norm.2',out=>'norm.2.out',args=>'-c s -a');
run_test(\&test_vcf_norm,$opts,in=>'norm.iupac',fai=>'norm.iupac',out=>'norm.iupac.out',args=>'-c s');
run_test(\&test_vcf_norm,$opts,in=>'norm.3',fai=>'norm.3',out=>'norm.3.out',args=>'-c s');
run_test(\&test_vcf_norm,$opts,in=>'norm.3',fai=>'norm.3',out=>'norm.3.2.out',args=>q[-c s -i'alt="N"']);
run_test(\&test_vcf_norm,$opts,in=>'atomize.split.1',out=>'atomize.split.1.0.out',args=>['-a --old-rec-tag OLD_REC','-m -any --force']);
run_test(\&test_vcf_norm,$opts,in=>'atomize.split.1',out=>'atomize.split.1.1.out',args=>['-m -any --old-rec-tag OLD_REC --force','-a']);
run_test(\&test_vcf_norm,$opts,in=>'atomize.split.1',out=>'atomize.split.1.1.out',args=>'-m -any --old-rec-tag OLD_REC --force -a');
run_test(\&test_vcf_norm,$opts,in=>'atomize.split.1',out=>'atomize.split.1.2.out',args=>'--atomize --atom-overlaps . --old-rec-tag OLD_REC');
run_test(\&test_vcf_norm,$opts,in=>'atomize.split.1',out=>'atomize.split.1.3.out',args=>'--atomize --old-rec-tag OLD_REC');
run_test(\&test_vcf_norm,$opts,in=>'atomize.split.2',out=>'atomize.split.2.1.out',args=>'--atomize --old-rec-tag OLD_REC');
run_test(\&test_vcf_norm,$opts,in=>'atomize.split.2',out=>'atomize.split.2.2.out',args=>'--atomize --atom-overlaps . --old-rec-tag OLD_REC');
run_test(\&test_vcf_norm,$opts,in=>'atomize.split.3',out=>'atomize.split.3.1.out',args=>'--atomize --atom-overlaps .');
run_test(\&test_vcf_norm,$opts,in=>'atomize.split.4',out=>'atomize.split.4.1.out',args=>'--atomize --atom-overlaps . --old-rec-tag OLD_REC');
run_test(\&test_vcf_norm,$opts,in=>'atomize.split.4',out=>'atomize.split.4.2.out',args=>q[--atomize --atom-overlaps . --old-rec-tag OLD_REC -i 'ILEN=0']);
run_test(\&test_vcf_norm,$opts,in=>'atomize.split.5',out=>'atomize.split.5.1.out',args=>q[--atomize --old-rec-tag OLD_REC --atom-overlaps .]);
run_test(\&test_vcf_norm,$opts,in=>'atomize.split.5',out=>'atomize.split.5.2.out',args=>q[--atomize --old-rec-tag OLD_REC]);
run_test(\&test_vcf_norm,$opts,in=>'norm.4',out=>'norm.4.1.out',args=>'-m +both');
run_test(\&test_vcf_norm,$opts,in=>'norm.4',out=>'norm.4.2.out',args=>'-m +any');
run_test(\&test_vcf_norm,$opts,in=>'norm.5',out=>'norm.5.1.out',args=>'-m - --multi-overlaps 0');
run_test(\&test_vcf_norm,$opts,in=>'norm.5',out=>'norm.5.2.out',args=>'-m - --multi-overlaps .');
run_test(\&test_vcf_norm,$opts,in=>'norm.m-any',out=>'norm.m-any.1.out',args=>'-m -any');
run_test(\&test_vcf_norm,$opts,in=>'norm.phased-split',out=>'norm.phased-split.1.out',args=>'-m -any');
run_test(\&test_vcf_norm,$opts,in=>'norm.phased-join',out=>'norm.phased-join.1.out',args=>'-m +any');
run_test(\&test_vcf_norm,$opts,in=>'norm.symbolic',fai=>'norm.symbolic',out=>'norm.symbolic.1.out',args=>'--old-rec-tag ORI');
run_test(\&test_vcf_norm,$opts,in=>'norm.symbolic.2',fai=>'norm.symbolic',out=>'norm.symbolic.2.out',args=>'--old-rec-tag ORI');
run_test(\&test_vcf_norm,$opts,in=>'norm.right-align',fai=>'norm.right-align',out=>'norm.right-align.1.out',args=>'--old-rec-tag ORI');
run_test(\&test_vcf_norm,$opts,in=>'norm.right-align',fai=>'norm.right-align',out=>'norm.right-align.2.out',args=>'--old-rec-tag ORI -g {PATH}/norm.right-align.gff');
run_test(\&test_vcf_norm,$opts,in=>'norm.atom-split-norm',fai=>'norm.atom-split-norm',out=>'norm.atom-split-norm.1.out',args=>'--old-rec-tag ORI -a -m -any');
run_test(\&test_vcf_norm,$opts,in=>'norm.string-tags',out=>'norm.string-tags.1.out',args=>'-m -any');
run_test(\&test_vcf_norm,$opts,in=>'norm.split.merge',out=>'norm.split.merge.1.out',args=>['-m -','-m +both']);
run_test(\&test_vcf_norm,$opts,in=>'norm.split.merge',out=>'norm.split.merge.2.out',args=>['-m -','-m +indels']);
run_test(\&test_vcf_norm,$opts,in=>'norm.split.merge',out=>'norm.split.merge.3.out',args=>['-m -','-m +snps']);
run_test(\&test_vcf_norm,$opts,in=>'norm.split.merge',out=>'norm.split.merge.3.out',args=>['-m -','-m +snps']);
run_test(\&test_vcf_norm,$opts,in=>'norm.split.merge',out=>'norm.split.merge.4.out',args=>['-m -','-m +any']);
run_test(\&test_vcf_norm,$opts,in=>'norm.split.merge',out=>'norm.split.merge.5.out',args=>q[-m - -i 'type="snp"']);
run_test(\&test_vcf_norm,$opts,in=>'norm.merge.4',out=>'norm.merge.4.1.out',args=>'-m +any');
run_test(\&test_vcf_norm,$opts,in=>'norm.merge.4',out=>'norm.merge.4.2.out',args=>'-m +both');
run_test(\&test_vcf_view,$opts,in=>'filter.string.1',out=>'filter.string.1.1.out',args=>q[-i 'INFO/TAG=@{PATH}/filter.string.1.txt']);
run_test(\&test_vcf_view,$opts,in=>'filter.string.1',out=>'filter.string.1.1.out',args=>q[-i 'FMT/TAG=@{PATH}/filter.string.1.txt']);
run_test(\&test_vcf_view,$opts,in=>'merge.gvcf.2.a',out=>'merge.gvcf.2.a.1.out',args=>'-HA');
run_test(\&test_vcf_view,$opts,in=>'merge.gvcf.2.a',out=>'merge.gvcf.2.a.2.out',args=>'-HAA');
run_test(\&test_vcf_view,$opts,in=>'weird-chr-names',out=>'weird-chr-names.1.out',args=>'',reg=>'-r 1');
run_test(\&test_vcf_view,$opts,in=>'weird-chr-names',out=>'weird-chr-names.1.out',args=>'',reg=>'-r 1:1-2');
run_test(\&test_vcf_view,$opts,in=>'weird-chr-names',out=>'weird-chr-names.1.out',args=>'',reg=>'-r 1:1,1:2');
run_test(\&test_vcf_view,$opts,in=>'weird-chr-names',out=>'weird-chr-names.2.out',args=>'',reg=>'-r 1:1-1');
run_test(\&test_vcf_view,$opts,in=>'weird-chr-names',out=>'weird-chr-names.3.out',args=>'',reg=>'-r {1:1}');
run_test(\&test_vcf_view,$opts,in=>'weird-chr-names',out=>'weird-chr-names.3.out',args=>'',reg=>'-r {1:1}:1-2');
run_test(\&test_vcf_view,$opts,in=>'weird-chr-names',out=>'weird-chr-names.3.out',args=>'',reg=>'-r {1:1}:1,{1:1}:2');
run_test(\&test_vcf_view,$opts,in=>'weird-chr-names',out=>'weird-chr-names.4.out',args=>'',reg=>'-r {1:1}:1-1');
run_test(\&test_vcf_view,$opts,in=>'weird-chr-names',out=>'weird-chr-names.5.out',args=>'',reg=>'-r {1:1-1}');
run_test(\&test_vcf_view,$opts,in=>'weird-chr-names',out=>'weird-chr-names.5.out',args=>'',reg=>'-r {1:1-1}:1-2');
run_test(\&test_vcf_view,$opts,in=>'weird-chr-names',out=>'weird-chr-names.5.out',args=>'',reg=>'-r {1:1-1}:1,{1:1-1}:2');
run_test(\&test_vcf_view,$opts,in=>'weird-chr-names',out=>'weird-chr-names.6.out',args=>'',reg=>'-r {1:1-1}:1-1');
run_test(\&test_vcf_view,$opts,in=>'weird-chr-names',args=>'',reg=>'-r {1:1-1}-2',expected_failure=>1);
run_test(\&test_vcf_view,$opts,in=>'view',out=>'view.1.out',args=>'-aUc1 -C1 -s NA00002 -v snps',reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view',out=>'view.2.out',args=>'-f PASS -Xks NA00003',reg=>'-r20,Y');
run_test(\&test_vcf_view,$opts,in=>'view',out=>'view.3.out',args=>'-xs NA00003',reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view',out=>'view.4.out',args=>q[-i 'QUAL==999 && (FS<20 || FS>=41.02) && ICF>-0.1 && HWE*2>1.2'],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view',out=>'view.5.out',args=>q[-p],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view',out=>'view.6.out',args=>q[-P],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view',out=>'view.7.out',args=>q[-hm2 -M2 -q0.3 -Q0.7],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view',out=>'view.8.out',args=>q[-Hu],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view',out=>'view.9.out',args=>q[-GVsnps],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view',out=>'view.10.out',args=>q[-ne 'INDEL=1 || PV4[0]<0.006'],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view',out=>'view.exclude.out',args=>'-s ^NA00003',reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.omitgenotypes',out=>'view.omitgenotypes.out',args=>'',reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.omitgenotypes',out=>'view.dropgenotypes.out',args=>'-G',reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.omitgenotypes',out=>'view.dropgenotypes.noheader.out',args=>'-HG',reg=>'');
run_test(\&test_vcf_view,$opts,in=>'many.alleles',out=>'many.alleles.trim.out',args=>'-a',reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.vectors',out=>'view.vectors.A.out',args=>'-asA',reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.vectors',out=>'view.vectors.B.out',args=>'-asB',reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.vectors.2',out=>'view.vectors.C.out',args=>'-asA',reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.filter',out=>'view.filter.1.out',args=>q[-H -i'FMT/FGS[*:0]="AAAAAA"'],reg=>'');    # test expressions
run_test(\&test_vcf_view,$opts,in=>'view.filter',out=>'view.filter.2.out',args=>q[-H -i'FMT/FGS[*:2]="C"'],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.filter',out=>'view.filter.3.out',args=>q[-H -i'FMT/FGS[*:4]="EE"'],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.filter',out=>'view.filter.4.out',args=>q[-H -i'FMT/FRS[*:1]="BB"'],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.filter',out=>'view.filter.5.out',args=>q[-H -i'TXT0="text"'],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.chrs',out=>'view.chrs.out',args=>'',reg=>'',tgts=>'view.chrs.tab');
run_test(\&test_vcf_view,$opts,in=>'filter.2',out=>'filter.11.out',args=>q[-i 'POS>=3062917'],reg=>'1:3062917-3157410');
run_test(\&test_vcf_view,$opts,in=>'idx.1',out=>'idx.1.out',args=>q[-H -r 1:10,1:12,1:10]);
run_test(\&test_vcf_view,$opts,in=>'idx.2',out=>'idx.2.out',args=>q[-H -r 1:1172777-1172804,1:1172806-1172808]);
run_test(\&test_vcf_view,$opts,in=>'idx.2',out=>'idx.2.out',args=>q[-H -R {PATH}/idx.2.bed]);
run_test(\&test_vcf_view,$opts,in=>'idx.3',out=>'idx.3.out',args=>q[-H -R {PATH}/idx.3.bed]);
run_test(\&test_vcf_view,$opts,in=>'idx.4',out=>'idx.4.out',args=>q[-H -R {PATH}/idx.4.bed]);
run_test(\&test_vcf_view,$opts,in=>'view-t',out=>'view-t.1.out',args=>'-Ht 2',reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view-t',out=>'view-t.2.out',args=>'-Ht 3',reg=>'');
run_test(\&test_vcf_view,$opts,in=>'overlap',out=>'overlap.0.out',args=>q[-H -r chr1:100-200 --regions-overlap 0]);
run_test(\&test_vcf_view,$opts,in=>'overlap',out=>'overlap.1.out',args=>q[-H -r chr1:100-200 --regions-overlap 1]);
run_test(\&test_vcf_view,$opts,in=>'overlap',out=>'overlap.2.out',args=>q[-H -r chr1:100-200 --regions-overlap 2]);
run_test(\&test_vcf_view,$opts,in=>'overlap',out=>'overlap.0.out',args=>q[-H -t chr1:100-200 --targets-overlap 0]);
run_test(\&test_vcf_view,$opts,in=>'overlap',out=>'overlap.1.out',args=>q[-H -t chr1:100-200 --targets-overlap 1]);
run_test(\&test_vcf_view,$opts,in=>'overlap',out=>'overlap.2.out',args=>q[-H -t chr1:100-200 --targets-overlap 2]);
run_test(\&test_vcf_view,$opts,in=>'overlap',out=>'overlap.neg0.out',args=>q[-H -t ^chr1:100-200 --targets-overlap 0]);
run_test(\&test_vcf_view,$opts,in=>'overlap',out=>'overlap.neg1.out',args=>q[-H -t ^chr1:100-200 --targets-overlap 1]);
run_test(\&test_vcf_view,$opts,in=>'overlap',out=>'overlap.neg2.out',args=>q[-H -t ^chr1:100-200 --targets-overlap 2]);
run_test(\&test_vcf_64bit,$opts,in=>'view64bit.1',out=>'view64bit.1.out',do_bcf=>1);
run_test(\&test_vcf_64bit,$opts,in=>'view64bit.2',out=>'view64bit.2.out',do_bcf=>1);
run_test(\&test_vcf_64bit,$opts,in=>'view64bit.3',out=>'view64bit.3.out');     # large coordinates don't work with BCF
run_test(\&test_vcf_64bit,$opts,in=>'view64bit.4',out=>'view64bit.4.out',do_bcf=>1);
run_test(\&test_vcf_64bit,$opts,in=>'view64bit.5',out=>'view64bit.5.out',do_bcf=>1);
run_test(\&test_vcf_view,$opts,in=>'view.minmaxac',out=>'view.minmaxac.1.out',args=>q[-H -C5:nonmajor],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.minmaxac',out=>'view.minmaxac.2.out',args=>q[-H -c6:nonmajor],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.minmaxac',out=>'view.minmaxac.1.out',args=>q[-H -q0.3:major],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.filter.annovar',out=>'view.filter.annovar.1.out',args=>q[-H -i 'Gene.refGene=="RAD21L1"'],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.filter.annovar',out=>'view.filter.annovar.2.out',args=>q[-H -i 'Gene.refGene~"NOD"'],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view.filter.annovar',out=>'view.filter.annovar.3.out',args=>q[-H -i 'LJB2_MutationTaster=="0.291000"'],reg=>'');
run_test(\&test_vcf_view,$opts,in=>'view-a',out=>'view-a.1.out',args=>q[-H -a]);
run_test(\&test_vcf_view,$opts,in=>'view.sites',out=>'view.sites.1.out',args=>'',tgts=>'view.sites.txt');
run_test(\&test_vcf_view,$opts,in=>'view.sites',out=>'view.sites.1.out',args=>'',tgts=>'view.sites.txt.gz');
run_test(\&test_vcf_head,$opts,in=>'mpileup.2.vcf',in_nheaders=>22);
run_test(\&test_vcf_head2,$opts,in=>'mpileup.2',out=>'head.1.out',args=>'-s0');
run_test(\&test_vcf_head2,$opts,in=>'mpileup.2',out=>'head.2.out',args=>'-s1');
run_test(\&test_vcf_head2,$opts,in=>'mpileup.2',out=>'head.3.out',args=>'-s2 -h2');
run_test(\&test_vcf_call,$opts,in=>'mpileup',out=>'mpileup.1.out',args=>'-mv');
run_test(\&test_vcf_call,$opts,in=>'mpileup',out=>'mpileup.2.out',args=>'-mg0');
run_test(\&test_vcf_call,$opts,in=>'mpileup',out=>'mpileup.3.out',args=>'-mv -S {PATH}/mpileup.3.samples');
run_test(\&test_vcf_call,$opts,in=>'mpileup',out=>'mpileup.4.out',args=>'-mv -S {PATH}/mpileup.4.samples');
run_test(\&test_vcf_call,$opts,in=>'mpileup',out=>'mpileup.5.out',args=>'-mv -S {PATH}/mpileup.5.samples');
run_test(\&test_vcf_call,$opts,in=>'mpileup.X',out=>'mpileup.X.out',args=>'-mv --ploidy-file {PATH}/mpileup.ploidy -S {PATH}/mpileup.samples');
run_test(\&test_vcf_call,$opts,in=>'mpileup.X',out=>'mpileup.X.out',args=>'-mv --ploidy-file {PATH}/mpileup.ploidy -S {PATH}/mpileup.ped');
run_test(\&test_vcf_call,$opts,in=>'mpileup.X',out=>'mpileup.X.2.out',args=>'-mv --ploidy-file {PATH}/mpileup.ploidy -S {PATH}/mpileup.2.samples');
run_test(\&test_vcf_call,$opts,in=>'mpileup.NA19213.NA19129',out=>'mpileup.hwe.1.out',args=>'-mv');
run_test(\&test_vcf_call,$opts,in=>'mpileup.NA19213.NA19129',out=>'mpileup.hwe.1b.out',args=>'-mv -G - --group-samples-tag AD');
run_test(\&test_vcf_call,$opts,in=>'mpileup.hwe',out=>'mpileup.hwe.2.out',args=>'-mv');
run_test(\&test_vcf_call,$opts,in=>'mpileup.hwe',out=>'mpileup.hwe.3.out',args=>'-mv -G - --group-samples-tag AD');                               # 21,3,0 becomes 0/0 because of the prior -P
run_test(\&test_vcf_call,$opts,in=>'mpileup.hwe',out=>'mpileup.hwe.4.out',args=>'-mv -G {PATH}/mpileup.hwe.samples --group-samples-tag AD');
run_test(\&test_vcf_call_cAls,$opts,in=>'mpileup',out=>'mpileup.cAls.out',tab=>'mpileup');
run_test(\&test_vcf_call_cAls,$opts,in=>'mpileup.2',out=>'mpileup.cAls.2.out',tab=>'mpileup.2');
run_test(\&test_vcf_call_cAls,$opts,in=>'mpileup.3',out=>'mpileup.cAls.3.out',tab=>'mpileup.3',args=>'-i');
run_test(\&test_vcf_call_cAls,$opts,in=>'mpileup.3',out=>'mpileup.cAls.4.out',tab=>'mpileup.4',args=>'-i');
run_test(\&test_vcf_call_cAls,$opts,in=>'mpileup.3',out=>'mpileup.cAls.5.out',tab=>'mpileup.5',args=>'-i');
run_test(\&test_vcf_call_cAls,$opts,in=>'mpileup.4',out=>'mpileup.cAls.6.out',tab=>'mpileup.6',args=>'-i');
run_test(\&test_vcf_call_cAls,$opts,in=>'mpileup.5',out=>'mpileup.cAls.7.out',tab=>'mpileup.7',args=>'-i');
run_test(\&test_vcf_call_cAls,$opts,in=>'mpileup.cals.1',out=>'mpileup.cals.8.out',tab=>'mpileup.cals.1',args=>'');
run_test(\&test_vcf_call_cAls,$opts,in=>'mpileup.cals.2',out=>'mpileup.cals.9.out',tab=>'mpileup.cals.2',args=>'');
run_test(\&test_vcf_call,$opts,in=>'mpileup.c',out=>'mpileup.c.1.out',args=>'-cv');
# run_test(\&test_vcf_call,$opts,in=>'mpileup.c',out=>'mpileup.c.2.out',args=>'-cg0');
run_test(\&test_vcf_call,$opts,in=>'mpileup.c.X',out=>'mpileup.c.X.out',args=>'-cv --ploidy-file {PATH}/mpileup.ploidy -S {PATH}/mpileup.samples');
run_test(\&test_vcf_call,$opts,in=>'mpileup.c.X',out=>'mpileup.c.X.out',args=>'-cv --ploidy-file {PATH}/mpileup.ploidy -S {PATH}/mpileup.ped');
run_test(\&test_vcf_call,$opts,in=>'mpileup.c.X',out=>'mpileup.c.X.2.out',args=>'-cv --ploidy-file {PATH}/mpileup.ploidy -S {PATH}/mpileup.2.samples');
run_test(\&test_vcf_call,$opts,in=>'call-G',out=>'call-G.1.out',args=>'-mv');
run_test(\&test_vcf_call,$opts,in=>'call-G',out=>'call-G.2.out',args=>'-mv -G - --group-samples-tag AD');
run_test(\&test_vcf_call,$opts,in=>'call-G.2',out=>'call-G.2.1.out',args=>'-mv -F AN_POP,AC_POP');
run_test(\&test_vcf_call,$opts,in=>'call.af-fixation',out=>'call.af-fixation.1.out',args=>'-m');
run_test(\&test_vcf_call,$opts,in=>'call.af-fixation',out=>'call.af-fixation.2.out',args=>'-m -G {PATH}/call.af-fixation.txt');
run_test(\&test_vcf_call,$opts,in=>'call.af-fixation',out=>'call.af-fixation.3.out',args=>'-m -G {PATH}/call.af-fixation.txt -a GP,GQ');
run_test(\&test_vcf_filter,$opts,in=>'view.filter',out=>'view.filter.6.out',args=>q[-S. -e'TXT0="text"'],reg=>'');
run_test(\&test_vcf_filter,$opts,in=>'view.filter',out=>'view.filter.7.out',args=>q[-S. -e'FMT/FRS[*:1]="BB"'],reg=>'');
run_test(\&test_vcf_filter,$opts,in=>'view.filter',out=>'view.filter.8.out',args=>q[-S. -e'FMT/FGS[*:0]="AAAAAA"'],reg=>'');
run_test(\&test_vcf_filter,$opts,in=>'view.filter',out=>'view.filter.9.out',args=>q[-S. -e'FMT/FGS[*:1]="BBB"'],reg=>'');
run_test(\&test_vcf_filter,$opts,in=>'view.filter',out=>'view.filter.10.out',args=>q[-S. -e'FMT/FGS[*:4]="EE"'],reg=>'');
run_test(\&test_vcf_filter,$opts,in=>'view.filter',out=>'view.filter.11.out',args=>q[-S. -e'FMT/STR="XX"'],reg=>'');
run_test(\&test_vcf_filter,$opts,in=>'view.filter.2',out=>'view.filter.12.out',args=>q[-S. -e'FMT/FILTER="aaa"'],reg=>'');
run_test(\&test_vcf_filter,$opts,in=>'filter.1',out=>'filter.1.out',args=>'-mx -g2 -G2');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.2.out',args=>q[-e'QUAL==59.2 || (INDEL=0 & (FMT/GQ=25 | FMT/DP=10))' -sModified -S.]);
run_test(\&test_vcf_filter,$opts,in=>'filter.3',out=>'filter.3.out',args=>q[-e'INFO/DP=19'],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.3',out=>'filter.4.out',args=>q[-e'INFO/DP=19' -s XX],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.3',out=>'filter.5.out',args=>q[-e'INFO/DP=19' -s XX -m+],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.3',out=>'filter.6.out',args=>q[-e'INFO/DP=19' -s XX -mx],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.3',out=>'filter.7.out',args=>q[-e'INFO/DP=19' -s XX -m+x],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.3',out=>'filter.3.out',args=>q[-e'FMT/GT="0/2"'],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.3',out=>'filter.4.out',args=>q[-e'FMT/GT="0/2"' -s XX],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.3',out=>'filter.5.out',args=>q[-e'FMT/GT="0/2"' -s XX -m+],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.3',out=>'filter.6.out',args=>q[-e'FMT/GT="0/2"' -s XX -mx],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.3',out=>'filter.7.out',args=>q[-e'FMT/GT="0/2"' -s XX -m+x],fmt=>'%POS\\t%FILTER\\t%DP[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.8.out',args=>q[-i'FMT/GT="0/0" && AC[*]=2'],fmt=>'%POS\\t%AC[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.8.out',args=>q[-i'AC[*]=2 && FMT/GT="0/0"'],fmt=>'%POS\\t%AC[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.9.out',args=>q[-i'ALT="."'],fmt=>'%POS\\t%AC[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.4',out=>'filter.10.out',args=>q[-S . -i 'FORMAT/TEST3<25']);
run_test(\&test_vcf_filter,$opts,in=>'filter.4',out=>'filter.10.out',args=>q[-S . -i 'FORMAT/TEST4<25']);
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.12.out',args=>q[-i'GT="A"'],fmt=>'%POS[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.13.out',args=>q[-i'GT="RR"'],fmt=>'%POS[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.14.out',args=>q[-i'GT="RA"'],fmt=>'%POS[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.14.out',args=>q[-i'GT="AR"'],fmt=>'%POS[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.15.out',args=>q[-i'GT="AA"'],fmt=>'%POS[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.16.out',args=>q[-i'GT="aA"'],fmt=>'%POS[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.16.out',args=>q[-i'GT="Aa"'],fmt=>'%POS[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.17.out',args=>q[-i'GT="HOM"'],fmt=>'%POS[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.18.out',args=>q[-i'GT="HET"'],fmt=>'%POS[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.19.out',args=>q[-i'GT="HAP"'],fmt=>'%POS[\\t%GT]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.5',out=>'filter.20.out',args=>q[-i'AD[:1]=11'],fmt=>'%POS[\\t%AD]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.5',out=>'filter.21.out',args=>q[-i'AD[1:]=11'],fmt=>'%POS[\\t%AD]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.5',out=>'filter.22.out',args=>q[-i'FR[0:1]=11'],fmt=>'%POS[\\t%FR]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.5',out=>'filter.23.out',args=>q[-i'AD[*]="."'],fmt=>'%POS[\\t%AD]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.5',out=>'filter.24.out',args=>q[-i'AD[0:0]=="."'],fmt=>'%POS[\\t%AD]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.5',out=>'filter.25.out',args=>q[-i'AD[0:0]!="."'],fmt=>'%POS[\\t%AD]\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.5',out=>'filter.26.out',args=>q[-i'QUAL=="."'],fmt=>'%POS\\t%QUAL\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.27.out',args=>q[-i'N_PASS(DP>32)=1'],fmt=>'[%POS\\t%SAMPLE\\t%DP\\n]');
run_test(\&test_vcf_filter,$opts,in=>'filter.2',out=>'filter.27.out',args=>q[-i'F_PASS(DP>32)=0.5'],fmt=>'[%POS\\t%SAMPLE\\t%DP\\n]');
run_test(\&test_vcf_filter,$opts,in=>'filter.6',out=>'filter.28.out',args=>q[-i'F_MISSING>=1/5'],fmt=>'%POS\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.6',out=>'filter.28.out',args=>q[-i'F_MISSING>=0.2'],fmt=>'%POS\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.6',out=>'filter.28.out',args=>q[-i'F_PASS(GT=="mis")>=1/5'],fmt=>'%POS\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.6',out=>'filter.28.out',args=>q[-i'F_PASS(GT=="mis")>=0.2'],fmt=>'%POS\\n');
run_test(\&test_vcf_filter,$opts,in=>'filter.7',out=>'filter.29.out',args=>'-mx -s + -g2:mnp,indel,other');
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.30.out',args=>q[-S . -e 'FORMAT/AO==4']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.30.out',args=>q[-S . -e 'MAX(FORMAT/AO[0:])==4']);    # although not desired, this is how MAX() works
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.31.out',args=>q[-S . -e 'MAX(FORMAT/AO)==4']);        # ( matches the desired sample but selects all samples, as opposed to SMPL_MAX)
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.30.out',args=>q[-S . -e 'MIN(FORMAT/AO[0:])==3']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.30.out',args=>q[-S . -e 'MIN(FORMAT/AO)==2']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.30.out',args=>q[-S . -e 'MIN(FORMAT/AO[0:])==3']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.30.out',args=>q[-S . -e 'AVG(FORMAT/AO[2:])==4']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.30.out',args=>q[-S . -e 'MEDIAN(FORMAT/AO[2:])==4']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.30.out',args=>q[-S . -e 'STDEV(FORMAT/AO[0:])=0.5']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.30.out',args=>q[-S . -e 'SUM(FORMAT/AO[0:])=7']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.32.out',args=>q[-S . -e 'SMPL_MAX(FORMAT/AO)==4']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.33.out',args=>q[-S . -e 'sMIN(FORMAT/AO)==2']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.33.out',args=>q[-S . -e 'ABS(sAVG(FORMAT/AO)-3.66666)<1e-5']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.34.out',args=>q[-S . -e 'sMEDIAN(FORMAT/AO)==4']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.33.out',args=>q[-S . -e 'ABS(sSTDEV(FORMAT/AO)-1.2472191)<1e-5']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.33.out',args=>q[-S . -e 'sSUM(FORMAT/AO)==11']);
run_test(\&test_vcf_filter,$opts,in=>'filter.9',out=>'filter.35.out',args=>q[-i 'QUAL/FMT/AD==55']);
run_test(\&test_vcf_filter,$opts,in=>'filter.9',out=>'filter.35.out',args=>q[-i 'QUAL/INFO/AD==10']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.36.out',args=>q[-S . -e 'ABS(SMPL_MAX(FORMAT/AO))=5']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.37.out',args=>q[-S . -e 'PHRED(AO[1:])>-4']);
run_test(\&test_vcf_filter,$opts,in=>'filter.8',out=>'filter.37.out',args=>q[-S . -e 'ABS(AO[1:])==2']);
run_test(\&test_vcf_filter,$opts,in=>'filter.10',out=>'filter.38.out',args=>q[-i 'sum(AD[*]) > FORMAT/DP']);
run_test(\&test_vcf_filter,$opts,in=>'filter.10',out=>'filter.38.out',args=>q[-i 'FORMAT/DP < sum(AD[*])']);
run_test(\&test_vcf_filter,$opts,in=>'filter.10',out=>'filter.39.out',args=>q[-i 'sum(AD[*]) < FORMAT/DP']);
run_test(\&test_vcf_filter,$opts,in=>'filter.10',out=>'filter.39.out',args=>q[-i 'FORMAT/DP > sum(AD[*])']);
run_test(\&test_vcf_filter,$opts,in=>'filter.1',out=>'filter.40.out',args=>q[--soft-filter xxx --mask 2:1005-1008 --mask-overlap 0]);
run_test(\&test_vcf_filter,$opts,in=>'filter.1',out=>'filter.41.out',args=>q[--soft-filter xxx --mask 2:1005-1008 --mask-overlap 1]);
run_test(\&test_vcf_filter,$opts,in=>'filter.1',out=>'filter.42.out',args=>q[--soft-filter xxx --mask 2:1005-1008 --mask-overlap 2]);
run_test(\&test_vcf_filter,$opts,in=>'filter.1',out=>'filter.43.out',args=>q[--soft-filter xxx --mask ^2:1005-1008]);
run_test(\&test_vcf_sort,$opts,in=>'sort',out=>'sort.out',args=>q[-m 0],fmt=>'%CHROM\\t%POS\\t%REF,%ALT\\n');
run_test(\&test_vcf_sort,$opts,in=>'sort',out=>'sort.out',args=>q[-m 1000],fmt=>'%CHROM\\t%POS\\t%REF,%ALT\\n');
run_test(\&test_vcf_regions,$opts,in=>'regions');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.escape.1',tab=>'annotate.escape.1',out=>'annotate.escape.1.1.out',args=>q[-c CHROM,POS,ISTR,FMT/FSTR]);
run_test(\&test_vcf_annotate,$opts,in=>'annotate.match.1',tab=>'annotate.match.1',out=>'annotate.match.1.1.out',args=>q[-c CHROM,POS,-,-,SCORE,~X,-,- -i'STR={X}']);
run_test(\&test_vcf_annotate,$opts,in=>'annotate.match.1',tab=>'annotate.match.1',out=>'annotate.match.1.2.out',args=>q[-c CHROM,POS,REF,ALT,SCORE,-,~X,- -i'INT={X}']);
run_test(\&test_vcf_annotate,$opts,in=>'annotate.match.1',tab=>'annotate.match.1',out=>'annotate.match.1.2.out',args=>q[-c CHROM,POS,REF,ALT,SCORE,-,-,~X -i'FLT={X}']);
run_test(\&test_vcf_annotate,$opts,in=>'annotate',tab=>'annotate',out=>'annotate.out',args=>'-c CHROM,POS,REF,ALT,ID,QUAL,INFO/T_INT,INFO/T_FLOAT,INDEL');
run_test(\&test_vcf_annotate,$opts,in=>'annotate',tab=>'annotate2',out=>'annotate2.out',args=>'-c CHROM,POS,-,T_STR');
run_test(\&test_vcf_annotate,$opts,in=>'annotate',tab=>'annotate2',out=>'annotate22.out',args=>'-c CHROM,FROM,TO,T_STR');
run_test(\&test_vcf_annotate,$opts,in=>'annotate',vcf=>'annots',out=>'annotate3.out',args=>'-c STR,ID,QUAL,FILTER');
run_test(\&test_vcf_annotate,$opts,in=>'annotate2',vcf=>'annots2',out=>'annotate4.out',args=>'-c ID,QUAL,FILTER,INFO,FMT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate2',vcf=>'annots2',out=>'annotate5.out',args=>'-c ID,QUAL,+FILTER,+INFO,FMT/GT -s A');
run_test(\&test_vcf_annotate,$opts,in=>'annotate2',vcf=>'annots2',out=>'annotate18.out',args=>'-c ID,QUAL,+FILTER,+INFO,FMT/GT -s "A B"');
run_test(\&test_vcf_annotate,$opts,in=>'annotate2',vcf=>'annots2',out=>'annotate19.out',args=>'-c ID,QUAL,+FILTER,+INFO,FMT/GT -s "A C"');
run_test(\&test_vcf_annotate,$opts,in=>'annotate2',vcf=>'annots2',out=>'annotate20.out',args=>'-c ID,QUAL,+FILTER,+INFO,FMT/GT -s "B C"');
run_test(\&test_vcf_annotate,$opts,in=>'annotate3',out=>'annotate6.out',args=>'-x ID,QUAL,^FILTER/fltA,FILTER/fltB,^INFO/AA,INFO/BB,^FMT/GT,FMT/PL');
run_test(\&test_vcf_annotate,$opts,in=>'annotate3',out=>'annotate7.out',args=>'-x FORMAT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate4',vcf=>'annots4',out=>'annotate8.out',args=>'-c +INFO');
run_test(\&test_vcf_annotate,$opts,in=>'annotate4',tab=>'annots4',out=>'annotate8.out',args=>'-c CHROM,POS,REF,ALT,+FA,+FR,+IA,+IR,+SA,+SR');
run_test(\&test_vcf_annotate,$opts,in=>'annotate10',tab=>'annots10',out=>'annotate10.out',args=>'-c CHROM,POS,FMT/FINT,FMT/FFLT,FMT/FSTR');
run_test(\&test_vcf_annotate,$opts,in=>'annotate2',vcf=>'annots2',out=>'annotate11.out',args=>'-c CHROM,POS,FMT/FINT,FMT/FFLT,FMT/FSTR -s A');
run_test(\&test_vcf_annotate,$opts,in=>'annotate2',tab=>'annots11',out=>'annotate11.out',args=>'-c CHROM,POS,FMT/FINT,FMT/FFLT,FMT/FSTR -s A');
run_test(\&test_vcf_annotate,$opts,in=>'annotate2',vcf=>'annots2',out=>'annotate12.out',args=>'-c AAA:=IINT,FMT/BBB:=FMT/FINT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate2',vcf=>'annots2',out=>'annotate13.out',args=>'-x INFO -c INFO/IINT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate2',vcf=>'annots2',out=>'annotate14.out',args=>q[-x INFO -c INFO/IINT -e'POS=3000001' -k]);
run_test(\&test_vcf_annotate,$opts,in=>'annotate11',vcf=>'annots11',out=>'annotate15.out',args=>q[-c FMT]);
run_test(\&test_vcf_annotate,$opts,in=>'annotate2',vcf=>'annots2',out=>'annotate16.out',args=>'-c FMT/newGT:=GT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate2',vcf=>'annots12',out=>'annotate17.out',args=>'-c FMT/GT:=newGT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate13',tab=>'annots13',out=>'annotate21.out',args=>'-c CHROM,BEG,END,ABC');
run_test(\&test_vcf_annotate,$opts,in=>'annotate13',tab=>'annots13',out=>'annotate23.out',args=>'-c CHROM,BEG,END,ABC -l ABC:append');
run_test(\&test_vcf_annotate,$opts,in=>'annotate13',tab=>'annots13',out=>'annotate24.out',args=>'-c CHROM,BEG,END,ABC -l ABC:unique');
run_test(\&test_vcf_annotate,$opts,in=>'annotate14',out=>'annotate25.out',args=>'-x FILTER/XX,INFO/XX --force');
run_test(\&test_vcf_annotate,$opts,in=>'annotate15',tab=>'annotate15',out=>'annotate26.out',args=>'-s SAMPLE1 -c CHROM,FROM,TO,FMT/FOO,BAR');
run_test(\&test_vcf_annotate,$opts,in=>'annotate15',tab=>'annotate15',out=>'annotate27.out',args=>'-s SAMPLE2 -c CHROM,FROM,TO,FMT/FOO,BAR');
run_test(\&test_vcf_annotate,$opts,in=>'annotate15',tab=>'annotate15',out=>'annotate27.out',args=>q[-s SAMPLE2 -c CHROM,FROM,TO,FMT/FOO,BAR -H '##FORMAT=<ID=FOO,Number=1,Type=String,Description="Some description">' -H '##INFO=<ID=BAR,Number=1,Type=Integer,Description="Some description">']);
run_test(\&test_vcf_annotate,$opts,in=>'annotate16',out=>'annotate28.out',args=>'-x FILTER');
run_test(\&test_vcf_annotate,$opts,in=>'annotate17.1',tab=>'annotate17.1',out=>'annotate17.1.out',args=>'-c CHROM,BEG,END,A,B -l A:append,B:append');
run_test(\&test_vcf_annotate,$opts,in=>'annotate17.2',tab=>'annotate17.1',out=>'annotate17.2.out',args=>'-c CHROM,BEG,END,A,B -l A:append,B:append');
run_test(\&test_vcf_annotate,$opts,in=>'annotate17.3',tab=>'annotate17.3',out=>'annotate17.3.out',args=>'-c CHROM,BEG,END,A,B -l A:append,B:append');
run_test(\&test_vcf_annotate,$opts,in=>'annotate18.1',tab=>'annotate18.1',out=>'annotate18.1.out',args=>'-c CHROM,BEG,END,A,B,C,D,E -l A:sum,B:avg,C:min,D:max,E:append');
run_test(\&test_vcf_annotate,$opts,in=>'annotate18.2',tab=>'annotate18.2',out=>'annotate18.2.out',args=>'-c CHROM,BEG,END,A,B,C,D,E -l A:sum,B:avg,C:min,D:max,E:append');
run_test(\&test_vcf_annotate,$opts,in=>'annotate19.dst',vcf=>'annotate19.src',out=>'annotate19.1.out',args=>'-c INFO/ID:=ID,INFO/INFO_ID:=INFO/ID,ID,=ID:=INFO/ID');
run_test(\&test_vcf_annotate,$opts,in=>'annotate19.dst',vcf=>'annotate19.src',out=>'annotate19.2.out',args=>'-c FILTER,INFO/FILTER:=FILTER,INFO/INFO_FILTER:=INFO/FILTER');
run_test(\&test_vcf_annotate,$opts,in=>'annotate19.dst',vcf=>'annotate19.src',out=>'annotate19.3.out',args=>'-c INFO/FILTER:=FILTER,INFO/INFO_FILTER:=INFO/FILTER');
run_test(\&test_vcf_annotate,$opts,in=>'annotate19.dst',                      out=>'annotate19.4.out',args=>'-c INFO/FILTER:=FILTER');
run_test(\&test_vcf_annotate,$opts,in=>'annotate19.dst',vcf=>'annotate19.src',out=>'annotate19.5.out',args=>'-c INFO/FILTER:=FILTER');
run_test(\&test_vcf_annotate,$opts,in=>'annotate19.dst',vcf=>'annotate19.src',out=>'annotate19.4.out',args=>'-c INFO/FILTER:=./FILTER');
run_test(\&test_vcf_annotate,$opts,in=>'annotate19.dst',vcf=>'annotate19.src',out=>'annotate19.6.out',args=>'-c INFO/FILTER:=./FILTER,FILTER');
run_test(\&test_vcf_annotate,$opts,in=>'annotate19.dst',vcf=>'annotate19.src',out=>'annotate19.7.out',args=>'-c FILTER,INFO/FILTER:=./FILTER');
run_test(\&test_vcf_annotate,$opts,in=>'annotate20.dst',vcf=>'annotate20.src',out=>'annotate20.1.out',args=>'-c  FMT/GT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate20.dst',vcf=>'annotate20.src',out=>'annotate20.2.out',args=>'-c +FMT/GT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate20.dst',vcf=>'annotate20.src',out=>'annotate20.3.out',args=>'-c -FMT/GT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.multi',tab=>'annotate.multi',out=>'annotate.multi.1.out',args=>'-c CHROM,POS,REF,ALT,ANN -l ANN:append');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.missing-append',tab=>'annotate.missing-append',out=>'annotate.missing-append.1.out',args=>'-c CHROM,POS,REF,ALT,STR,INT,FLT -l STR:append-missing,INT:append-missing,FLT:append-missing');
run_test(\&test_vcf_annotate,$opts,in=>'annotate9',tab=>'annots9',out=>'annotate9.out',args=>'-c CHROM,POS,REF,ALT,+ID');
run_test(\&test_vcf_annotate,$opts,in=>'annotate21',out=>'annotate29.out',args=>'--rename-annots {PATH}/annotate21.txt');
run_test(\&test_vcf_annotate,$opts,in=>'annotate21',out=>'annotate29.out',args=>'-c XX:=FORMAT/X-X,YY:=FORMAT/Y-Y,AA:=FORMAT/A-A,AA:=INFO/A-A,BB:=INFO/B-B,XX:=INFO/X-X,YY:=INFO/Y-Y,fltA:=FILTER/flt-A,fltB:=FILTER/flt-B,fltX:=FILTER/flt-X,fltY:=FILTER/flt-Y');
run_test(\&test_vcf_annotate,$opts,in=>'annotate22',vcf=>'annotate22',out=>'annotate30.out',args=>'-c FMT/XX,INFO/XX -x FMT/XX,INFO/XX');
run_test(\&test_vcf_annotate,$opts,in=>'annotate23',tab=>'annotate23',out=>'annotate31.out',args=>'-c CHROM,POS,~ID,REF,ALT,INFO/END');
run_test(\&test_vcf_annotate,$opts,in=>'annotate24.dst',vcf=>'annotate24.src',out=>'annotate24.1.out',args=>'-c XX');
run_test(\&test_vcf_annotate,$opts,in=>'annotate25',tab=>'annotate25',out=>'annotate25.1.out',args=>'-c CHROM,POS,ID,REF,ALT,~INFO/END');
run_test(\&test_vcf_annotate,$opts,in=>'annotate26',tab=>'annotate26',out=>'annotate26.1.out',args=>'-c CHROM,POS,-POS');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.missing',tab=>'annotate.missing',out=>'annotate.missing.1.out',args=>'-c CHROM,POS,REF,ALT,TSTR,TFLT,TINT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.missing',tab=>'annotate.missing',out=>'annotate.missing.2.out',args=>'-c CHROM,POS,REF,ALT,.TSTR,.TFLT,.TINT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.missing',tab=>'annotate.missing',out=>'annotate.missing.3.out',args=>'-c CHROM,POS,REF,ALT,.+TSTR,.+TFLT,.+TINT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.missing',tab=>'annotate.missing',out=>'annotate.missing.4.out',args=>'-c CHROM,POS,REF,ALT,+TSTR,+TFLT,+TINT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.missing',tab=>'annotate.missing',out=>'annotate.missing.5.out',args=>'-c CHROM,POS,REF,ALT,.=TSTR,.=TFLT,.=TINT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.missing',tab=>'annotate.missing',out=>'annotate.missing.6.out',args=>'-c CHROM,POS,REF,ALT,=TSTR,=TFLT,=TINT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.olap',tab=>'annots.olap',out=>'annotate.olap.1.out',args=>'-c CHROM,BEG,END,DB -l DB:unique');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.olap',tab=>'annots.olap',out=>'annotate.olap.2.out',args=>'-c CHROM,BEG,END,DB -l DB:unique --min-overlap 0.4:0.5 -m XXX');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.id',vcf=>'annots.id',out=>'annotate.id.1.out',args=>'-c ALT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.id',vcf=>'annots.id',out=>'annotate.id.2.out',args=>'-c +ALT');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.id.2',vcf=>'annots.id.2',out=>'annotate.id.2.1.out',args=>'--pair-logic some -c +ID');
run_test(\&test_vcf_annotate,$opts,in=>'annotate.id.2',vcf=>'annots.id.2',out=>'annotate.id.2.2.out',args=>'--pair-logic any -c +ID');
run_test(\&test_vcf_annotate,$opts,in=>'annotate27',tab=>'annotate27',out=>'annotate.32.out',args=>'-c CHROM,POS,REF,ALT,EVIDENCE');
run_test(\&test_vcf_annotate,$opts,in=>'annotate28',tab=>'annots28',out=>'annotate28.1.out',args=>'-c CHROM,POS,REF,ALT,FMT/TEST -s smpl1,smpl2');
run_test(\&test_vcf_annotate,$opts,in=>'annotate28',tab=>'annots28',out=>'annotate28.2.out',args=>'-c CHROM,POS,REF,ALT,FMT/TEST -s smpl2,smpl1');
run_test(\&test_vcf_annotate,$opts,in=>'annotate28',tab=>'annots28',out=>'annotate28.3.out',args=>'-c CHROM,POS,REF,ALT,FMT/TEST -s smpl1');
run_test(\&test_vcf_annotate,$opts,in=>'annotate28',tab=>'annots28',out=>'annotate28.4.out',args=>'-c CHROM,POS,REF,ALT,FMT/TEST -s smpl2');
run_test(\&test_vcf_annotate,$opts,in=>'annotate',out=>'annotate.33.out',args=>'-m XXX');
run_test(\&test_vcf_annotate,$opts,in=>'annotate34',tab=>'annots34',out=>'annotate34.out',args=>q[-c CHROM,FROM,TO,INFO/END -H '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate in reference for SV">']);
run_test(\&test_vcf_annotate,$opts,in=>'annots-mark',bed=>'annots-mark',out=>'annots-mark.1.out',args=>q[-c CHROM,FROM,TO -m TAG]);
run_test(\&test_vcf_plugin,$opts,in=>'checkploidy',out=>'checkploidy.out',cmd=>'+check-ploidy --no-version');
run_test(\&test_vcf_plugin,$opts,in=>'checkploidy.2',out=>'checkploidy.2.out',cmd=>'+check-ploidy --no-version');
run_test(\&test_vcf_plugin,$opts,in=>'checkploidy.2',out=>'checkploidy.3.out',cmd=>'+check-ploidy --no-version',args=>'-- -m');
run_test(\&test_vcf_plugin,$opts,in=>'plugin1',out=>'missing2ref.out',cmd=>'+missing2ref --no-version');
run_test(\&test_vcf_plugin,$opts,in=>'plugin1',out=>'missing2ref.out',cmd=>'+setGT --no-version',args=>'-- -t . -n 0');
run_test(\&test_vcf_plugin,$opts,in=>'setGT',out=>'setGT.1.out',cmd=>'+setGT --no-version',args=>'-- -t q -n 0 -i \'GT~"." && FMT/DP=30 && GQ=150\'');
run_test(\&test_vcf_plugin,$opts,in=>'setGT.2',out=>'setGT.2.out',cmd=>'+setGT --no-version',args=>'-- -t q -n . -i \'GT[@{QPATH}/setGT.samples.txt]="het"\'');
run_test(\&test_vcf_plugin,$opts,in=>'setGT.2',out=>'setGT.3.out',cmd=>'+setGT --no-version',args=>'-- -t q -n . -i \'GT[@{QPATH}/setGT.samples.txt]="het" & binom(AD[@{QPATH}/setGT.samples.txt])<0.1\'');
run_test(\&test_vcf_plugin,$opts,in=>'setGT.3',out=>'setGT.3.1.out',cmd=>'+setGT --no-version',args=>'-- -t a -n pM');
run_test(\&test_vcf_plugin,$opts,in=>'setGT.3',out=>'setGT.3.2.out',cmd=>'+setGT --no-version',args=>'-- -t a -n pm');
run_test(\&test_vcf_plugin,$opts,in=>'setGT.3',out=>'setGT.3.3.out',cmd=>'+setGT --no-version',args=>'-- -t a -n c:1');
run_test(\&test_vcf_plugin,$opts,in=>'setGT.3',out=>'setGT.3.4.out',cmd=>'+setGT --no-version',args=>'-- -t a -n c:"1|1"');
run_test(\&test_vcf_plugin,$opts,in=>'setGT.3',out=>'setGT.3.5.out',cmd=>'+setGT --no-version',args=>'-- -t a -n c:"m|M"');
run_test(\&test_vcf_plugin,$opts,in=>'setGT.3',out=>'setGT.3.6.out',cmd=>'+setGT --no-version',args=>'-- -t a -n c:0/1/1');
run_test(\&test_vcf_plugin,$opts,in=>'setGT.4',out=>'setGT.4.1.out',cmd=>'+setGT --no-version',args=>q[-- -t q -n . -e 'FMT/DP>90']);
run_test(\&test_vcf_plugin,$opts,in=>'setGT.4',out=>'setGT.4.2.out',cmd=>'+setGT --no-version',args=>q[-- -t q -n . -e 'FMT/DP>100']);
run_test(\&test_vcf_plugin,$opts,in=>'setGT.5',out=>'setGT.5.1.out',cmd=>'+setGT --no-version',args=>q[-- -t a -n X]);
run_test(\&test_vcf_plugin,$opts,in=>'setGT.5',out=>'setGT.5.2.out',cmd=>'+setGT --no-version',args=>q[-- -t a -n c:0/X/X]);
run_test(\&test_vcf_plugin,$opts,in=>'setGT.6',out=>'setGT.6.1.out',cmd=>'+setGT --no-version',args=>q[-- -t ./x -n .]);
run_test(\&test_vcf_plugin,$opts,in=>'setGT.6',out=>'setGT.6.1.out',cmd=>'+setGT --no-version',args=>q[-- -t . -n .]);
run_test(\&test_vcf_plugin,$opts,in=>'setGT.6',out=>'setGT.6.2.out',cmd=>'+setGT --no-version',args=>q[-- -t r:0.5 -n .]);
run_test(\&test_vcf_plugin,$opts,in=>'plugin1',out=>'fill-AN-AC.out',cmd=>'+fill-AN-AC --no-version');
run_test(\&test_vcf_plugin,$opts,in=>'dosage',out=>'dosage.1.out',cmd=>'+dosage',args=>'-- -t PL');
run_test(\&test_vcf_plugin,$opts,in=>'dosage',out=>'dosage.2.out',cmd=>'+dosage',args=>'-- -t GL');
run_test(\&test_vcf_plugin,$opts,in=>'dosage',out=>'dosage.3.out',cmd=>'+dosage',args=>'-- -t GT');
run_test(\&test_vcf_plugin,$opts,in=>'fixploidy',out=>'fixploidy.out',cmd=>'+fixploidy --no-version',args=>'-- -s {PATH}/fixploidy.samples -p {PATH}/fixploidy.ploidy');
run_test(\&test_vcf_plugin,$opts,in=>'view.PL',out=>'guess-ploidy.PL.out',cmd=>'+guess-ploidy',args=>'-vrX | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'view.GL',out=>'guess-ploidy.GL.out',cmd=>'+guess-ploidy',args=>'-vrX | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'view.GL',out=>'view.PL.vcf',cmd=>'+tag2tag --no-version',args=>'-- -r --gl-to-pl');
run_test(\&test_vcf_plugin,$opts,in=>'view.GL',out=>'view.GL-GP.vcf',cmd=>'+tag2tag --no-version',args=>'-- --gl-to-gp');
run_test(\&test_vcf_plugin,$opts,in=>'view.GP',out=>'view.GT.vcf',cmd=>'+tag2tag --no-version',args=>'-- -r --gp-to-gt -t 0.2');
run_test(\&test_vcf_plugin,$opts,in=>'tag2tag.LPL.1',out=>'tag2tag.LPL.1.1.vcf',cmd=>'+tag2tag --no-version',args=>'-- --LXX-to-XX');
run_test(\&test_vcf_plugin,$opts,in=>'tag2tag.LPL.1',out=>'tag2tag.LPL.1.2.vcf',cmd=>'+tag2tag --no-version',args=>'-- --LXX-to-XX -r');
run_test(\&test_vcf_plugin,$opts,in=>'tag2tag.LPL.1',out=>'tag2tag.LPL.1.3.vcf',cmd=>'+tag2tag --no-version',args=>'-- --LXX-to-XX -r -d AD:0,PL:255 -s 3');
run_test(\&test_vcf_plugin,$opts,in=>'query.variantkey',out=>'query.add-variantkey.vcf',cmd=>'+add-variantkey',args=>'');
run_test(\&test_vcf_plugin,$opts,in=>'query.variantkey',out=>'variantkey-hex.out',cmd=>'+variantkey-hex',args=>'test/');
run_test(\&test_vcf_plugin,$opts,in=>'query.nucleotide',out=>'query.allele-length.tsv',cmd=>'+allele-length',args=>'');
run_test(\&test_vcf_plugin,$opts,in=>'merge.a',out=>'fill-tags.out',cmd=>'+fill-tags --no-version',args=>'-- -t AN,AC,AC_Hom,AC_Het,AC_Hemi');
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'fill-tags.2.out',cmd=>'+fill-tags --no-version',args=>'-- -t AC,AN,AF,MAF,NS');
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'fill-tags.3.out',cmd=>'+fill-tags --no-version',args=>'-- -t AC -S {PATH}/fill-tags.3.smpl');
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'fill-tags.5.out',cmd=>'+fill-tags --no-version',args=>'-- -t "DP:1=int(sum(FORMAT/DP))" -S {PATH}/fill-tags.3.smpl');
run_test(\&test_vcf_plugin,$opts,in=>'many-alts',out=>'fill-tags.4.out',cmd=>'+fill-tags --no-version',args=>'-- -t AN,AC');
run_test(\&test_vcf_plugin,$opts,in=>'fill-tags-hemi',out=>'fill-tags-hemi.1.out',cmd=>'+fill-tags --no-version');
run_test(\&test_vcf_plugin,$opts,in=>'fill-tags-hemi',out=>'fill-tags-hemi.2.out',cmd=>'+fill-tags --no-version',args=>'-- -d');
run_test(\&test_vcf_plugin,$opts,in=>'fill-tags-hwe',out=>'fill-tags-hwe.out',cmd=>'+fill-tags --no-version');
run_test(\&test_vcf_plugin,$opts,in=>'fill-tags-hwe',out=>'fill-tags-func.out',cmd=>'+fill-tags --no-version',args=>q[-- -t 'XX:1=F_PASS(GT="alt")']);
run_test(\&test_vcf_plugin,$opts,in=>'fill-tags-AN0',out=>'fill-tags-AN0.out',cmd=>'+fill-tags --no-version',args=>'-- -t all,END,TYPE,F_MISSING');
run_test(\&test_vcf_plugin,$opts,in=>'fill-tags-VAF',out=>'fill-tags-VAF.out',cmd=>'+fill-tags --no-version',args=>'-- -t VAF,VAF1');
run_test(\&test_vcf_plugin,$opts,in=>'fill-tags-AD',out=>'fill-tags-AD.1.out',cmd=>'+fill-tags --no-version',args=>q[-- -t 'INFO/DP:1=int(sum(FMT/AD))']);
run_test(\&test_vcf_plugin,$opts,in=>'fill-tags-AD',out=>'fill-tags-AD.2.out',cmd=>'+fill-tags --no-version',args=>q[-- -t 'INFO/DP:1=int(sum(INFO/AD))']);
run_test(\&test_vcf_plugin,$opts,in=>'fill-tags-AD',out=>'fill-tags-AD.3.out',cmd=>'+fill-tags --no-version',args=>q[-- -t 'FORMAT/DP:1=int(smpl_sum(FMT/AD))']);
run_test(\&test_vcf_plugin,$opts,in=>'fill-tags-AD',out=>'fill-tags-AD.4.out',cmd=>'+fill-tags --no-version',args=>q[-- -t 'XX=N_PASS(FMT/AD[:0]<=10)','YY=N_PASS(FMT/AD[:0]>10)']);
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'view.GTisec.out',cmd=>'+GTisec',args=>' | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'view.GTisec.H.out',cmd=>'+GTisec',args=>'-- -H | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'view.GTisec.Hm.out',cmd=>'+GTisec',args=>'-- -Hm | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'view.GTisec.Hmv.out',cmd=>'+GTisec',args=>'-- -Hmv | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'view.GTisec.Hv.out',cmd=>'+GTisec',args=>'-- -Hv | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'view.GTisec.m.out',cmd=>'+GTisec',args=>'-- -m | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'view.GTisec.mv.out',cmd=>'+GTisec',args=>'-- -mv | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'view.GTisec.v.out',cmd=>'+GTisec',args=>'-- -v | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'trio',out=>'trio.out',cmd=>'+trio-switch-rate',args=>'-- -p {PATH}/trio.ped | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'trio-stats',out=>'trio-stats.out',cmd=>'+trio-stats',args=>'-a 1 -p {PATH}/trio-stats.ped -d mendel-errors,transmitted | grep -v ^CMD');
run_test(\&test_vcf_plugin,$opts,in=>'trio-stats',out=>'trio-stats.2.out',cmd=>'+trio-stats',args=>'-p {PATH}/trio-stats.ped -d mendel-errors,transmitted | grep -v ^CMD');
run_test(\&test_vcf_plugin,$opts,in=>'indel-stats',out=>'smpl-stats.1.out',cmd=>'+smpl-stats',args=>'| grep -v ^CMD');
run_test(\&test_vcf_plugin,$opts,in=>'indel-stats',out=>'indel-stats.1.out',cmd=>'+indel-stats',args=>'| grep -v ^CMD');
run_test(\&test_vcf_plugin,$opts,in=>'indel-stats',out=>'indel-stats.2.out',cmd=>'+indel-stats',args=>' -p {PATH}/trio-stats.ped | grep -v ^CMD');
run_test(\&test_vcf_plugin,$opts,in=>'indel-stats',out=>'indel-stats.3.out',cmd=>'+indel-stats',args=>' -p {PATH}/trio-stats.2.ped | grep -v ^CMD');
run_test(\&test_vcf_plugin,$opts,in=>'ad-bias',out=>'ad-bias.out',cmd=>'+ad-bias',args=>'-- -s {PATH}/ad-bias.samples | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'ad-bias.2',out=>'ad-bias.out',cmd=>'+ad-bias',args=>'-- -s {PATH}/ad-bias.samples | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'ad-bias',out=>'ad-bias.2.out',cmd=>'+ad-bias',args=>'--no-version -- -s {PATH}/ad-bias.samples -c | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'ad-bias.2',out=>'ad-bias.2.out',cmd=>'+ad-bias',args=>'--no-version -- -s {PATH}/ad-bias.samples -c | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'af-dist',out=>'af-dist.out',cmd=>'+af-dist',args=>' | grep -v bcftools');
run_test(\&test_vcf_plugin,$opts,in=>'fixref.2a',out=>'fixref.2.out',index=>['fixref.2b'],cmd=>'+fixref',args=>'-- -f {PATH}/norm.fa -i {TMP}/fixref.2b.vcf.gz');
run_test(\&test_vcf_plugin,$opts,in=>'fixref.3',out=>'fixref.3.out',cmd=>'+fixref',args=>'-- -f {PATH}/fixref.3.fa -m top');
run_test(\&test_vcf_plugin,$opts,in=>'fixref.2a',out=>'fixref.4.out',index=>['fixref.2b'],cmd=>'+fixref',args=>'-- -f {PATH}/norm.fa -m ref-alt');
run_test(\&test_vcf_plugin,$opts,in=>'fixref.2a',out=>'fixref.5.out',index=>['fixref.2b'],cmd=>'+fixref',args=>'-- -f {PATH}/norm.fa -m flip');
run_test(\&test_vcf_plugin,$opts,in=>'fixref.2a',out=>'fixref.6.out',cmd=>'+fixref',args=>'-- -f {PATH}/norm.fa -m flip-all');
run_test(\&test_vcf_plugin,$opts,in=>'fixref.2a',out=>'fixref.7.out',cmd=>'+fixref',args=>'-- -f {PATH}/norm.fa -m swap');
run_test(\&test_vcf_plugin,$opts,in=>'aa',out=>'aa.out',cmd=>'+fill-from-fasta',args=>'-- -f {PATH}/aa.fa -c AA -h {PATH}/aa.hdr -i \'TYPE="snp"\'');
run_test(\&test_vcf_plugin,$opts,in=>'aa',out=>'aa.2.out',cmd=>'+fill-from-fasta',args=>'-- -f {PATH}/aa.fa -c REF -N');
run_test(\&test_vcf_plugin,$opts,in=>'ref',out=>'ref.out',cmd=>'+fill-from-fasta',args=>'-- -f {PATH}/norm.fa -c REF');
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'view.GTsubset.NA1.out',cmd=>'+GTsubset --no-version',args=>'-- -s NA00001');
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'view.GTsubset.NA1NA2.out',cmd=>'+GTsubset --no-version',args=>'-- -s NA00001,NA00002');
run_test(\&test_vcf_plugin,$opts,in=>'view',out=>'view.GTsubset.NA1NA2NA3.out',cmd=>'+GTsubset --no-version',args=>'-- -s NA00001,NA00002,NA00003');
run_test(\&test_vcf_plugin,$opts,in=>'mendelian',out=>'mendelian.1.out',cmd=>'+mendelian2',args=>'-p child1,dad1,mom1 -md');
run_test(\&test_vcf_plugin,$opts,in=>'mendelian',out=>'mendelian.6.out',cmd=>'+mendelian2',args=>'-p child1,dad1,mom1 -mg');
run_test(\&test_vcf_plugin,$opts,in=>'mendelian',out=>'mendelian.3.out',cmd=>'+mendelian2',args=>'-p child1,dad1,mom1 -me');
run_test(\&test_vcf_plugin,$opts,in=>'mendelian',out=>'mendelian.4.out',cmd=>'+mendelian2',args=>'-p child1,dad1,mom1 -ma');
run_test(\&test_vcf_plugin,$opts,in=>'mendelian',out=>'mendelian.7.out',cmd=>'+mendelian2',args=>'-p child1,dad1,mom1 -mm');
run_test(\&test_vcf_plugin,$opts,in=>'mendelian',out=>'mendelian.8.out',cmd=>'+mendelian2',args=>'-p child1,dad1,mom1 | grep -v ^#');
run_test(\&test_vcf_plugin,$opts,in=>'contrast',out=>'contrast.out',cmd=>'+contrast',args=>'-a PASSOC,FASSOC,NOVELAL,NOVELGT -0 a,b -1 c');
run_test(\&test_vcf_plugin,$opts,in=>'contrast',out=>'contrast.out',cmd=>'+contrast',args=>'-a PASSOC,FASSOC,NOVELAL,NOVELGT -0 {PATH}/contrast0.txt -1 {PATH}/contrast1.txt');
run_test(\&test_vcf_plugin,$opts,in=>'contrast',out=>'contrast.1.out',cmd=>'+contrast',args=>'-a NASSOC -0 a,b,c -1 d --force-samples');
run_test(\&test_vcf_plugin,$opts,in=>'contrast.1',out=>'contrast.1.1.out',cmd=>'+contrast',args=>'-a NOVELAL,NOVELGT -0 A -1 B');
run_test(\&test_vcf_plugin,$opts,in=>'contrast.1',out=>'contrast.1.2.out',cmd=>'+contrast',args=>'-a NOVELGT -0 A -1 B');
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.1',out=>'trio-dnm/trio-dnm.1.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother --ppl --dnm-tag DNM:log | $$opts{bin}/bcftools query -f'[\\t%DNM]\\t[\\t%VAF]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.2',out=>'trio-dnm/trio-dnm.1.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother --ppl --dnm-tag DNM:log --force-AD | $$opts{bin}/bcftools query -f'[\\t%DNM]\\t[\\t%VAF]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.4',out=>'trio-dnm/trio-dnm.4.1.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother --use-DNG | $$opts{bin}/bcftools query -f'[\\t%DNM]\\t[\\t%VAF]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.4',out=>'trio-dnm/trio-dnm.4.1.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother           | $$opts{bin}/bcftools query -f'[\\t%DNM]\\t[\\t%VAF]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.4',out=>'trio-dnm/trio-dnm.4.2.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother --use-DNG --dnm-tag DNM:log | $$opts{bin}/bcftools query -f'[\\t%DNM]\\t[\\t%VAF]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.4',out=>'trio-dnm/trio-dnm.4.2.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother           --dnm-tag DNM:log | $$opts{bin}/bcftools query -f'[\\t%DNM]\\t[\\t%VAF]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.5',out=>'trio-dnm/trio-dnm.5.1.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother --use-DNG --dnm-tag DNM:log | $$opts{bin}/bcftools query -f'[\\t%DNM]\\t[\\t%VAF]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.5',out=>'trio-dnm/trio-dnm.5.1.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother           --dnm-tag DNM:log | $$opts{bin}/bcftools query -f'[\\t%DNM]\\t[\\t%VAF]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.6',out=>'trio-dnm/trio-dnm.6.1.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother --use-DNG --dnm-tag DNM:log | $$opts{bin}/bcftools query -f'[\\t%DNM]\\t[\\t%VAF]\\n'"); # incorrect miss by DNG
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.6',out=>'trio-dnm/trio-dnm.6.2.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother           --dnm-tag DNM:log | $$opts{bin}/bcftools query -f'[\\t%DNM]\\t[\\t%VAF]\\t[\\t%VA]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.7',out=>'trio-dnm/trio-dnm.7.1.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother --use-DNG --dnm-tag DNM:log | $$opts{bin}/bcftools query -f'[\\t%DNM]\\t[\\t%VAF]\\n'"); # incorrect miss, low PL
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.7',out=>'trio-dnm/trio-dnm.7.1.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother        --dnm-tag DNM:log | $$opts{bin}/bcftools query -f'[\\t%DNM]\\t[\\t%VAF]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.8',out=>'trio-dnm/trio-dnm.8.1.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother  | $$opts{bin}/bcftools query -f'[\\t%DNM]\\t[\\t%VAF]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.9',out=>'trio-dnm/trio-dnm.9.1.out',cmd=>'+trio-dnm2',args=>"-p 1X:proband,father,mother --use-NAIVE | $$opts{bin}/bcftools query -f'[\\t%DNM]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.9',out=>'trio-dnm/trio-dnm.9.2.out',cmd=>'+trio-dnm2',args=>"-p 2X:proband,father,mother --use-NAIVE | $$opts{bin}/bcftools query -f'[\\t%DNM]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.10',out=>'trio-dnm/trio-dnm.10.1.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother --with-pAD | $$opts{bin}/bcftools query -f'[\\t%DNM][\\t%VAF]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.11',out=>'trio-dnm/trio-dnm.11.1.out',cmd=>'+trio-dnm2',args=>"-p proband,father,mother | $$opts{bin}/bcftools query -f'%CHROM:%POS  DNM=[%DNM ]\\tAD=[%AD ]\\tQS=[%QS ]\\tVAF=[%VAF ]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'trio-dnm/trio-dnm.11',out=>'trio-dnm/trio-dnm.11.2.out',cmd=>'+trio-dnm2',args=>"-p 1X:proband,father,mother --strictly-novel | $$opts{bin}/bcftools query -f'%CHROM:%POS  DNM=[%DNM ]\\tAD=[%AD ]\\tQS=[%QS ]\\tVAF=[%VAF ]\\n'");
run_test(\&test_vcf_plugin,$opts,in=>'gvcfz',out=>'gvcfz.1.out',cmd=>'+gvcfz',args=>qq[-g 'PASS:GT!="alt"' -a | $$opts{bin}/bcftools query -f'%POS\\t%REF\\t%ALT\\t%END[\\t%GT][\\t%DP][\\t%GQ][\\t%RGQ]\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'gvcfz',out=>'gvcfz.2.out',cmd=>'+gvcfz',args=>qq[-g 'PASS:GQ>10; FLT:-' -a | $$opts{bin}/bcftools query -f'%POS\\t%REF\\t%ALT\\t%FILTER\\t%END[\\t%GT][\\t%DP][\\t%GQ][\\t%RGQ]\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'gvcfz.2',out=>'gvcfz.2.1.out',cmd=>'+gvcfz',args=>qq[-g 'PASS:GT!="alt"' -a | $$opts{bin}/bcftools query -f'%POS\\t%REF\\t%ALT\\t%FILTER\\t%END[\\t%GT][\\t%DP]\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'remove-overlaps.1',out=>'remove-overlaps.1.1.out',cmd=>'+remove-overlaps',args=>'-m overlap');
run_test(\&test_vcf_plugin,$opts,in=>'remove-overlaps.1',out=>'remove-overlaps.1.2.out',cmd=>'+remove-overlaps',args=>'-m overlap -M overlap');
run_test(\&test_vcf_plugin,$opts,in=>'remove-overlaps.1',out=>'remove-overlaps.1.3.out',cmd=>'+remove-overlaps',args=>'-m overlap -O t');
run_test(\&test_vcf_plugin,$opts,in=>'remove-overlaps.1',out=>'remove-overlaps.1.4.out',cmd=>'+remove-overlaps',args=>'-m overlap --reverse');
run_test(\&test_vcf_plugin,$opts,in=>'remove-overlaps.1',out=>'remove-overlaps.1.5.out',cmd=>'+remove-overlaps',args=>'-m dup -M DUP');
run_test(\&test_vcf_plugin,$opts,in=>'remove-overlaps.1',out=>'remove-overlaps.1.6.out',cmd=>'+remove-overlaps',args=>'-m dup -M unique --reverse');
run_test(\&test_vcf_plugin,$opts,in=>'remove-overlaps.2',out=>'remove-overlaps.2.1.out',cmd=>'+remove-overlaps',args=>q[-m 'min(QUAL)' -M rmme]);
run_test(\&test_vcf_plugin,$opts,in=>'remove-overlaps.3',out=>'remove-overlaps.3.1.out',cmd=>'+remove-overlaps',args=>q[-m 'min(QUAL)' -M rmme]);
run_test(\&test_vcf_plugin,$opts,in=>'remove-overlaps.3',out=>'remove-overlaps.3.1.out',cmd=>'+remove-overlaps',args=>q[-m 'min(QUAL)' -M rmme --missing 0]);
run_test(\&test_vcf_plugin,$opts,in=>'remove-overlaps.3',out=>'remove-overlaps.3.2.out',cmd=>'+remove-overlaps',args=>q[-m 'min(QUAL)' -M rmme --missing DP]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.1.out',cmd=>'+split-vep',args=>qq[-c Consequence -s worst:missense+ | $$opts{bin}/bcftools query -f'%POS\\t%Consequence\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.2.out',cmd=>'+split-vep',args=>qq[-c Consequence -s worst:missense+ | $$opts{bin}/bcftools query -f'%POS\\t%Consequence\\n' -i'Consequence!="."']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.2.out',cmd=>'+split-vep',args=>qq[-s worst:missense+ -f'%POS\\t%Consequence\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.3.out',cmd=>'+split-vep',args=>qq[-s primary:missense+ -f'%POS\\t%Consequence\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.3.out',cmd=>'+split-vep',args=>qq[-s CANONICAL=YES:missense+ -f'%POS\\t%Consequence\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.4.out',cmd=>'+split-vep',args=>qq[-s primary:missense+ -f'%POS\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.4.out',cmd=>'+split-vep',args=>qq[-s CANONICAL=YES:missense+ -f'%POS\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.2',out=>'split-vep.5.out',cmd=>'+split-vep',args=>qq[-s worst -f'%POS\\t%AF\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.2',out=>'split-vep.6.out',cmd=>'+split-vep',args=>qq[-s worst -f'%POS\\t%AF\\n' -a BCSQ]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.2',out=>'split-vep.6.out',cmd=>'+split-vep',args=>qq[-s worst -f'%POS\\t%INFO/AF\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.2',out=>'split-vep.6.out',cmd=>'+split-vep',args=>qq[-s worst -f'%POS\\t%INFO/AF\\n' -a BCSQ]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.3',out=>'split-vep.7.out',cmd=>'+split-vep',args=>qq[-s worst -f'%POS\\t%Consequence\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.3',out=>'split-vep.8.out',cmd=>'+split-vep',args=>qq[-s worst -f'[%POS\\t%SAMPLE\\t%GT\\t%Consequence\\n]' -i'GT="alt"']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.9.out', cmd=>'+split-vep',args=>qq[-t 1:14464 -f '%POS\\t%CANONICAL\\t%Consequence\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.10.out',cmd=>'+split-vep',args=>qq[-t 1:14464 -f '%POS\\t%CANONICAL\\t%Consequence\\n' -d]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.11.out',cmd=>'+split-vep',args=>qq[-t 1:14464 -f '%POS\\t%CSQ\\n' -A tab]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.12.out',cmd=>'+split-vep',args=>qq[-t 1:14464 -f '%POS\\t%CSQ\\n' -A tab -d]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.12.2.out',cmd=>'+split-vep',args=>qq[-t 1:14464 -f '%POS\\t%CSQ\\n' -A tab -d -H]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.12.3.out',cmd=>'+split-vep',args=>qq[-t 1:14464 -f '%POS\\t%CSQ\\n' -A tab -d -HH]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.4',out=>'split-vep.13.out',cmd=>'+split-vep',args=>qq[-f '%POS\\t%BCSQ\\n' -a BCSQ -A tab -d]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.4',out=>'split-vep.13.out',cmd=>'+split-vep',args=>qq[-f '%POS\\t%BCSQ\\n'         -A tab -d]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.14.out',cmd=>'+split-vep',args=>qq[-c gnomAD_NFE_AF:real,ALLELE_NUM:int | $$opts{bin}/bcftools query -f'%POS\\t%gnomAD_NFE_AF\\t%ALLELE_NUM\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.14.out',cmd=>'+split-vep',args=>qq[-c gnomAD_NFE_AF,ALLELE_NUM          | $$opts{bin}/bcftools query -f'%POS\\t%gnomAD_NFE_AF\\t%ALLELE_NUM\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.5',out=>'split-vep.15.out',cmd=>'+split-vep',args=>qq[-s :synonymous    -c Consequence | $$opts{bin}/bcftools query -f'%POS\\t%Consequence\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.5',out=>'split-vep.16.out',cmd=>'+split-vep',args=>qq[-s :synonymous -x -c Consequence | $$opts{bin}/bcftools query -f'%POS\\t%Consequence\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.6',out=>'split-vep.17.out',cmd=>'+split-vep',args=>qq[-c SAS_AF | grep ID=SAS_AF]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.6',out=>'split-vep.18.out',cmd=>'+split-vep',args=>qq[-c - | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.6',out=>'split-vep.19.out',cmd=>'+split-vep',args=>qq[-c - -s worst | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.7',out=>'split-vep.20.out',cmd=>'+split-vep',args=>qq[--annotation 'ANN' -c IMPACT -i 'INFO/IMPACT[*] ~ "MODIFIER"' | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.8',out=>'split-vep.21.out',cmd=>'+split-vep',args=>qq[-a BCSQ -s worst -c Consequence,amino_acid_change,BCSQ | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.8',out=>'split-vep.21.out',cmd=>'+split-vep',args=>qq[        -s worst -c Consequence,amino_acid_change,BCSQ | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.8',out=>'split-vep.22.out',cmd=>'+split-vep',args=>qq[-a BCSQ -s worst -f '%POS\\t%Consequence\\t%amino_acid_change\\t%BCSQ\\n' | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.8',out=>'split-vep.22.out',cmd=>'+split-vep',args=>qq[        -s worst -f '%POS\\t%Consequence\\t%amino_acid_change\\t%BCSQ\\n' | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.9',out=>'split-vep.23.out',cmd=>'+split-vep',args=>qq[-a BCSQ -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%BCSQ\\n' -d -A tab | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.9',out=>'split-vep.24.out',cmd=>'+split-vep',args=>qq[-a BCSQ -f '%xPOS %xConsequence\\n' -p x | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.10',out=>'split-vep.25.out',cmd=>'+split-vep',args=>qq[-a CSQ -f '%M_CAP_pred %M_CAP_score\\n' | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.10',out=>'split-vep.25.out',cmd=>'+split-vep',args=>qq[-a CSQ -f '%xM_CAP_pred %xM_CAP_score\\n' -p x | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.10',out=>'split-vep.25.out',cmd=>'+split-vep',args=>qq[       -f '%xM_CAP_pred %xM_CAP_score\\n' -p x | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.gene-list',out=>'split-vep.gene-list.1.out',cmd=>'+split-vep',args=>qq[-d -f '%CHROM:%POS %Gene %Consequence\\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.gene-list',out=>'split-vep.gene-list.2.out',cmd=>'+split-vep',args=>qq[-d -f '%CHROM:%POS %Gene %Consequence\\n' -g {PATH}/split-vep.gene-list.txt]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.gene-list',out=>'split-vep.gene-list.2.out',cmd=>'+split-vep',args=>qq[-d -f '%CHROM:%POS %Gene %Consequence\\n' -g {PATH}/split-vep.mixed-list.txt --gene-list-fields Feature,SYMBOL]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.gene-list',out=>'split-vep.gene-list.3.out',cmd=>'+split-vep',args=>qq[-d -f '%CHROM:%POS %Gene %Consequence\\n' -g +{PATH}/split-vep.gene-list.txt]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.gene-list',out=>'split-vep.gene-list.3.out',cmd=>'+split-vep',args=>qq[-d -f '%CHROM:%POS %Gene %Consequence\\n' -g +{PATH}/split-vep.mixed-list.txt --gene-list-fields Feature,SYMBOL]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.broken-LoF',out=>'split-vep.broken-LoF.out',cmd=>'+split-vep',args=>qq[-d -f '%CHROM:%POS %Consequence %LoF_info\\n' -a vep]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.broken-LoF',out=>'split-vep.broken-LoF.2.out',cmd=>'+split-vep',args=>qq[-d -f '%CHROM:%POS %LoF_info\\n' -a vep -i 'Consequence=="frameshift_variant"']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.26.out',cmd=>'+split-vep',args=>qq[-f'%POS\\n' -i'SYMBOL~"SAMD11"']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep',out=>'split-vep.27.out',cmd=>'+split-vep',args=>qq[-f'%POS\\t%MAX_AF\\n' -i'MAX_AF>0.999' -c MAX_AF,MAX_AF:float,MAX_AF]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.filter',out=>'split-vep.filter.1.out',cmd=>'+split-vep',args=>qq[-s worst -i'CSQ~"nonsense"' -f '%POS %Consequence %Feature %BIOTYPE']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.filter',out=>'split-vep.filter.2.out',cmd=>'+split-vep',args=>qq[-s worst -i'CSQ~"nonsense"' -f '%POS %Consequence %Feature %BIOTYPE %CSQ']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.select-tr-expr',out=>'select-tr-expr.1.out',cmd=>'+split-vep',args=>qq[-f '%POS %CANONICAL %MANE_SELECT %Consequence' -s primary]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.select-tr-expr',out=>'select-tr-expr.1.out',cmd=>'+split-vep',args=>qq[-f '%POS %CANONICAL %MANE_SELECT %Consequence' -s CANONICAL=YES]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.select-tr-expr',out=>'select-tr-expr.2.out',cmd=>'+split-vep',args=>qq[-f '%POS %CANONICAL %MANE_SELECT %Consequence' -s CANONICAL!=YES]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.select-tr-expr',out=>'select-tr-expr.2.out',cmd=>'+split-vep',args=>qq[-f '%POS %CANONICAL %MANE_SELECT %Consequence' -s CANONICAL=]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.select-tr-expr',out=>'select-tr-expr.1.out',cmd=>'+split-vep',args=>qq[-f '%POS %CANONICAL %MANE_SELECT %Consequence' -s mane]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.select-tr-expr',out=>'select-tr-expr.1.out',cmd=>'+split-vep',args=>qq[-f '%POS %CANONICAL %MANE_SELECT %Consequence' -s MANE_SELECT!='']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.select-tr-expr',out=>'select-tr-expr.3.out',cmd=>'+split-vep',args=>qq[-s mane]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.select-tr-expr',out=>'select-tr-expr.3.out',cmd=>'+split-vep',args=>qq[-s primary]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.select-tr-expr',out=>'select-tr-expr.3.out',cmd=>'+split-vep',args=>qq[-s PolyPhen~damaging]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.types',out=>'split-vep.types.1.out',cmd=>'+split-vep',args=>qq[-f'%cDNA_position %CDS_position %Protein_position']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.types',out=>'split-vep.types.2.out',cmd=>'+split-vep',args=>qq[-f'%cDNA_position %CDS_position %Protein_position' -c cDNA_position:int,CDS_position:int,Protein_position:int]);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.filter-missing-bt2098',out=>'split-vep.filter-missing-bt2098.1.out',cmd=>'+split-vep',args=>qq[-e 'gnomad_genomes_AF >= 0.02' | $$opts{bin}/bcftools query -f'%CHROM %POS \t %gnomad_genomes_AF\n']);
run_test(\&test_vcf_plugin,$opts,in=>'split-vep.filter-missing-bt2098',out=>'split-vep.filter-missing-bt2098.1.out',cmd=>'+split-vep',args=>qq[-i 'gnomad_genomes_AF < 0.02 || gnomad_genomes_AF=="."' | $$opts{bin}/bcftools query -f'%CHROM %POS \t %gnomad_genomes_AF\n']);
run_test(\&test_vcf_plugin,$opts,in=>'parental-origin',out=>'parental-origin.1.out',cmd=>'+parental-origin',args=>qq[-r 20:100 -p proband,father,mother -t del | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'parental-origin',out=>'parental-origin.2.out',cmd=>'+parental-origin',args=>qq[-r 20:101 -p proband,father,mother -t del | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'parental-origin',out=>'parental-origin.3.out',cmd=>'+parental-origin',args=>qq[-r 20:102 -p proband,father,mother -t del | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'parental-origin',out=>'parental-origin.4.out',cmd=>'+parental-origin',args=>qq[-r 20:103 -p proband,father,mother -t dup | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'parental-origin',out=>'parental-origin.5.out',cmd=>'+parental-origin',args=>qq[-r 20:104 -p proband,father,mother -t dup | grep -v ^#]);
run_test(\&test_vcf_plugin,$opts,in=>'prune.1',out=>'prune.1.1.out',cmd=>'+prune -w 1 -a r2,LD,HD');                        # annotate with r2,LD,HD at the previous site
run_test(\&test_vcf_plugin,$opts,in=>'prune.2',out=>'prune.2.1.out',cmd=>'+prune -w 1 -a r2,LD,HD');
run_test(\&test_vcf_plugin,$opts,in=>'prune.1',out=>'prune.1.2.out',cmd=>'+prune -w 2 -a r2 -m 0.5 -f MaxR2');              # prune within 2bp, max r2=0.5, soft filter
run_test(\&test_vcf_plugin,$opts,in=>'prune.1',out=>'prune.1.3.out',cmd=>'+prune -w 2 -a r2 -m 0.5 ');                      # prune within 2bp, max r2=0.5
run_test(\&test_vcf_plugin,$opts,in=>'prune.1',out=>'prune.1.4.out',cmd=>'+prune -w 2bp -n 1 --AF-tag AF');                 # leave 1 site within 2bp windows, prioritize by AF
run_test(\&test_vcf_plugin,$opts,in=>'prune.1',out=>'prune.1.5.out',cmd=>q[+prune -w 2bp -n 1 --AF-tag AF -i 'GT="alt"']);  # same as above but first discard REF-only sites
run_test(\&test_vcf_plugin,$opts,in=>'prune.1',out=>'prune.1.6.out',cmd=>'+prune -w 2bp -n 1 -N 1st');
run_test(\&test_vcf_plugin,$opts,in=>'prune.1',out=>'prune.1.7.out',cmd=>'+prune -w 2bp -n 1 -N rand --random-seed 1');
run_test(\&test_vcf_plugin,$opts,in=>'variant-distance',out=>'variant-distance.1.out',cmd=>'+variant-distance');
run_test(\&test_vcf_plugin,$opts,in=>'variant-distance',out=>'variant-distance.1.out',cmd=>'+variant-distance -d nearest');
run_test(\&test_vcf_plugin,$opts,in=>'variant-distance',out=>'variant-distance.2.out',cmd=>'+variant-distance -d fwd');
run_test(\&test_vcf_plugin,$opts,in=>'variant-distance',out=>'variant-distance.3.out',cmd=>'+variant-distance -d rev');
run_test(\&test_vcf_plugin,$opts,in=>'variant-distance',out=>'variant-distance.4.out',cmd=>'+variant-distance -d both');
run_test(\&test_plugin_split,$opts,in=>'split.1',out=>'split.1.1.out',tmp=>'split.1.1');
run_test(\&test_plugin_split,$opts,in=>'split.1',out=>'split.1.2.out',tmp=>'split.1.2',args=>'-S {PATH}/split.smpl.1.2.txt');
run_test(\&test_plugin_split,$opts,in=>'split.1',out=>'split.1.3.out',tmp=>'split.1.3',args=>'-S {PATH}/split.smpl.1.3.txt');
run_test(\&test_plugin_split,$opts,in=>'split.1',out=>'split.1.4.out',tmp=>'split.1.4',args=>q[-S {PATH}/split.smpl.1.3.txt -i 'GT[0]="alt"']);
run_test(\&test_plugin_split,$opts,in=>'split.1',out=>'split.1.5.out',tmp=>'split.1.5',args=>q[-S {PATH}/split.smpl.1.3.txt -i 'GT="alt"']);
run_test(\&test_plugin_split,$opts,in=>'split.1',out=>'split.1.6.out',tmp=>'split.1.6',args=>q[-S {PATH}/split.smpl.1.4.txt -i 'GT="alt"']);
run_test(\&test_plugin_split,$opts,in=>'split.1',out=>'split.1.7.out',tmp=>'split.1.7',args=>q[-G {PATH}/split.grp.1.1.txt]);
run_test(\&test_plugin_split,$opts,in=>'split.2',out=>'split.2.1.out',tmp=>'split.2.1',args=>q[]);
run_test(\&test_plugin_scatter,$opts,in=>'scatter.1',out=>'scatter.1.1.out',tmp=>'scatter.1.1',args=>q[-n 3]);
run_test(\&test_plugin_scatter,$opts,in=>'scatter.1',out=>'scatter.1.2.out',tmp=>'scatter.1.2',args=>q[-s 21,22]);
run_test(\&test_plugin_scatter,$opts,in=>'scatter.1',out=>'scatter.1.3.out',tmp=>'scatter.1.3',args=>q[-s 21,22 -x X]);
run_test(\&test_vcf_concat,$opts,in=>['concat.1.a','concat.1.b'],out=>'concat.1.vcf.out',do_bcf=>0,args=>'');
run_test(\&test_vcf_concat,$opts,in=>['concat.1.a','concat.1.b'],out=>'concat.1.bcf.out',do_bcf=>1,args=>'');
run_test(\&test_vcf_concat,$opts,in=>['concat.2.a','concat.2.b'],out=>'concat.2.vcf.out',do_bcf=>0,args=>'-a');
run_test(\&test_vcf_concat,$opts,in=>['concat.2.a','concat.2.b'],out=>'concat.2.bcf.out',do_bcf=>1,args=>'-a');
run_test(\&test_vcf_concat,$opts,in=>['concat.2.a','concat.2.b'],out=>'concat.4.vcf.out',do_bcf=>0,args=>'-aD');
run_test(\&test_vcf_concat,$opts,in=>['concat.2.a','concat.2.b'],out=>'concat.4.bcf.out',do_bcf=>1,args=>'-aD');
run_test(\&test_vcf_concat,$opts,in=>['concat.3.a','concat.3.b','concat.3.0','concat.3.c','concat.3.d','concat.3.e','concat.3.f'],out=>'concat.3.vcf.out',do_bcf=>0,args=>'-l --ligate-warn');
run_test(\&test_vcf_concat,$opts,in=>['concat.3.a','concat.3.b','concat.3.0','concat.3.c','concat.3.d','concat.3.e','concat.3.f'],out=>'concat.3.bcf.out',do_bcf=>1,args=>'-l --ligate-warn');
run_test(\&test_naive_concat,$opts,name=>'naive_concat',max_hdr_lines=>10000,max_body_lines=>10000,nfiles=>10);
run_test(\&test_vcf_concat,$opts,in=>['concat.4.a','concat.4.b'],out=>'concat.5.out',do_bcf=>0,args=>'-l');
run_test(\&test_vcf_concat,$opts,in=>['concat.4.a','concat.4.b'],out=>'concat.5.out',do_bcf=>1,args=>'-l');
run_test(\&test_vcf_concat,$opts,in=>['concat.5.a','concat.5.b','concat.5.c'],out=>'concat.5.1.out',do_bcf=>0,args=>'-l --ligate-warn');
run_test(\&test_vcf_concat,$opts,in=>['concat.5.a','concat.5.b','concat.5.c'],out=>'concat.5.1.out',do_bcf=>1,args=>'-l --ligate-warn');
run_test(\&test_vcf_concat,$opts,in=>['concat.5.a','concat.5.b','concat.5.c'],out=>'concat.5.2.out',do_bcf=>1,args=>'-l --ligate-force');
run_test(\&test_vcf_concat,$opts,in=>['concat.5.a','concat.5.b','concat.5.c'],out=>'concat.5.3.out',do_bcf=>0,args=>'-G -a -D');
run_test(\&test_vcf_concat,$opts,in=>['concat.5.a','concat.5.b','concat.5.c'],out=>'concat.5.3.out',do_bcf=>1,args=>'-G -a -D');
run_test(\&test_vcf_reheader,$opts,in=>'reheader',out=>'reheader.1.out',header=>'reheader.hdr');
run_test(\&test_vcf_reheader,$opts,in=>'reheader',out=>'reheader.2.out',samples=>'reheader.samples');
run_test(\&test_vcf_reheader,$opts,in=>'reheader',out=>'reheader.2.out',samples=>'reheader.samples2');
run_test(\&test_vcf_reheader,$opts,in=>'reheader',out=>'reheader.3.out',samples=>'reheader.samples3');
run_test(\&test_vcf_reheader,$opts,in=>'reheader',out=>'reheader.4.out',samples=>'reheader.samples4');
run_test(\&test_vcf_reheader,$opts,in=>'empty',out=>'reheader.empty.out',header=>'reheader.empty.hdr');
run_test(\&test_vcf_reheader,$opts,in=>'reheader.2',out=>'reheader.5.out',args=>'-f {PATH}/reheader.fai',nostdin=>1);
run_test(\&test_vcf_reheader,$opts,in=>'reheader.2',out=>'reheader.5.out',args=>'-h {PATH}/reheader.2.hdr -f {PATH}/reheader.fai',nostdin=>1);
run_test(\&test_vcf_reheader,$opts,in=>'reheader.3',out=>'reheader.6.out',args=>'-f {PATH}/reheader.3.fai',nostdin=>1);
run_test(\&test_rename_chrs,$opts,in=>'annotate');
run_test(\&test_vcf_convert,$opts,in=>'convert',out=>'convert.gs.gt.gen',args=>'-g -,.');
run_test(\&test_vcf_convert,$opts,in=>'convert',out=>'convert.gs.gt.ids.gen',args=>'-g -,. --vcf-ids');
run_test(\&test_vcf_convert,$opts,in=>'convert',out=>'convert.gs.gt.ids.gen6',args=>'-g -,. --vcf-ids --3N6');
run_test(\&test_vcf_convert,$opts,in=>'convert',out=>'convert.gs.gt.samples',args=>'-g .,-');
run_test(\&test_vcf_convert_hs2vcf,$opts,h=>'convert.gs.gt.ids.gen',s=>'convert.gs.gt.samples',out=>'convert.gs.vcf',args=>' --vcf-ids -G');
run_test(\&test_vcf_convert_hs2vcf,$opts,h=>'convert.gs.gt.ids.gen',s=>'convert.gs.gt.samples',out=>'convert.gs.noids.vcf',args=>'-G');
run_test(\&test_vcf_convert_hs2vcf,$opts,h=>'convert.gs.gt.ids.3N6.gen',s=>'convert.gs.gt.samples',out=>'convert.gs.noids.vcf',args=>'--3N6 -G');
run_test(\&test_vcf_convert_hs2vcf,$opts,h=>'convert.gs.gt.ids.gen.rev',s=>'convert.gs.gt.samples',out=>'convert.gs.vcf',args=>'--vcf-ids -G');
run_test(\&test_vcf_convert_hs2vcf,$opts,h=>'convert.gs.gt.ids.gen.rev',s=>'convert.gs.gt.samples',out=>'convert.gs.noids.vcf',args=>'-G');
run_test(\&test_vcf_convert,$opts,in=>'convert',out=>'convert.gs.pl.gen',args=>'-g -,. --tag PL');
run_test(\&test_vcf_convert,$opts,in=>'convert',out=>'convert.gs.pl.samples',args=>'-g .,- --tag PL');
run_test(\&test_vcf_convert,$opts,in=>'check',out=>'check.gs.vcfids.gen',args=>'-g -,. --vcf-ids');
run_test(\&test_vcf_convert,$opts,in=>'check',out=>'check.gs.vcfids.samples',args=>'-g .,- --vcf-ids');
run_test(\&test_vcf_convert,$opts,in=>'check',out=>'check.gs.chrom.gen',args=>'-g -,. --3N6');
run_test(\&test_vcf_convert,$opts,in=>'check',out=>'check.gs.chrom.samples',args=>'-g .,- --3N6');
run_test(\&test_vcf_convert,$opts,in=>'check',out=>'check.gs.vcfids_chrom.gen',args=>'-g -,. --3N6 --vcf-ids');
run_test(\&test_vcf_convert,$opts,in=>'check',out=>'check.gs.vcfids_chrom.samples',args=>'-g .,- --3N6 --vcf-ids');
run_test(\&test_vcf_convert,$opts,in=>'convert',out=>'convert.hls.haps',args=>'-h -,.,.');
run_test(\&test_vcf_convert,$opts,in=>'convert',out=>'convert.hls.legend',args=>'-h .,-,.');
run_test(\&test_vcf_convert,$opts,in=>'convert',out=>'convert.hls.ids.legend',args=>'-h .,-,. --vcf-ids');
run_test(\&test_vcf_convert,$opts,in=>'convert',out=>'convert.hls.samples',args=>'-h .,.,-');
run_test(\&test_vcf_convert_hls2vcf,$opts,h=>'convert.hls.gt.hap',l=>'convert.hls.gt.legend',s=>'convert.hls.gt.samples',out=>'convert.gt.noHead.vcf',args=>'-H');
run_test(\&test_vcf_convert,$opts,in=>'convert.hap-missing',out=>'convert.hap-missing.haps',args=>'--haplegendsample -,.,.');
run_test(\&test_vcf_convert,$opts,in=>'convert',out=>'convert.hs.hap',args=>'--hapsample -,.');
run_test(\&test_vcf_convert,$opts,in=>'convert',out=>'convert.hs.ids.hap',args=>'--hapsample -,. --vcf-ids');
run_test(\&test_vcf_convert,$opts,in=>'convert',out=>'convert.hs.sample',args=>'--hapsample .,-');
run_test(\&test_vcf_convert_hs2vcf,$opts,h=>'convert.hs.gt.hap',s=>'convert.hs.gt.samples',out=>'convert.gt.noHead.vcf',args=>'--hapsample2vcf');
run_test(\&test_vcf_convert_hs2vcf,$opts,h=>'convert.hs.gt.ids.hap',s=>'convert.hs.gt.samples',out=>'convert.gt.noHead.ids.vcf',args=>'--vcf-ids --hapsample2vcf');
run_test(\&test_vcf_convert_gvcf,$opts,in=>'convert.gvcf',out=>'convert.gvcf.out',fa=>'gvcf.fa',args=>'--gvcf2vcf -i\'FILTER="PASS"\'');
run_test(\&test_vcf_convert_tsv2vcf,$opts,in=>'convert.23andme',out=>'convert.23andme.vcf',args=>'-c ID,CHROM,POS,AA -s SAMPLE1',fai=>'23andme');
run_test(\&test_vcf_convert_tsv2vcf,$opts,in=>'convert.tsv',out=>'convert.tsv.vcf',args=>'-c -,CHROM,POS,REF,ALT',fai=>'23andme');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.overlaps.1',out=>'consensus.overlaps.1.1.out',fa=>'consensus.overlaps.1.fa',args=>'-s A');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.overlaps.1',out=>'consensus.overlaps.1.2.out',fa=>'consensus.overlaps.1.fa',args=>'-s B');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.overlaps.1',out=>'consensus.overlaps.1.3.out',fa=>'consensus.overlaps.1.fa',args=>'-s A -a N');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.overlaps.1',out=>'consensus.overlaps.1.4.out',fa=>'consensus.overlaps.1.fa',args=>'-s B -a N');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.beyond',out=>'consensus.beyond.1.out',fa=>'consensus.beyond.fa',args=>'');
run_test(\&test_vcf_consensus,$opts,in=>'consensus',out=>'consensus.1.out',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-s -');
run_test(\&test_vcf_consensus_chain,$opts,in=>'consensus',out=>'consensus.1.chain',chain=>'consensus.1.chain',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-s -');
run_test(\&test_vcf_consensus,$opts,in=>'consensus',out=>'consensus.2.out',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-H 1');
run_test(\&test_vcf_consensus_chain,$opts,in=>'consensus',out=>'consensus.2.chain',chain=>'consensus.2.chain',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-H 1');
run_test(\&test_vcf_consensus,$opts,in=>'consensus',out=>'consensus.3.out',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-I -s -');
run_test(\&test_vcf_consensus,$opts,in=>'consensus',out=>'consensus.16.out',fa=>'consensus.fa',args=>'-s - -I -m {PATH}/consensus.tab --mask-with X -m {PATH}/consensus.tab --mask-with lc');
run_test(\&test_vcf_consensus_chain,$opts,in=>'consensus',out=>'consensus.3.chain',chain=>'consensus.3.chain',fa=>'consensus.fa',mask=>'consensus.tab',args=>'-I');
run_test(\&test_vcf_consensus,$opts,in=>'consensus',out=>'consensus.4.out',fa=>'consensus.fa',args=>'-H 1');
run_test(\&test_vcf_consensus_chain,$opts,in=>'consensus',out=>'consensus.4.chain',chain=>'consensus.4.chain',fa=>'consensus.fa',args=>'-H 1');
run_test(\&test_vcf_consensus,$opts,in=>'consensus2',out=>'consensus2.1.out',fa=>'consensus2.fa',args=>'-H 1');
run_test(\&test_vcf_consensus,$opts,in=>'consensus2',out=>'consensus2.2.out',fa=>'consensus2.fa',args=>'-H 2');
run_test(\&test_vcf_consensus,$opts,in=>'empty',out=>'consensus.5.out',fa=>'consensus.fa',args=>'-s -');
run_test(\&test_vcf_consensus,$opts,in=>'consensus3',out=>'consensus3.out',fa=>'consensus2.fa',args=>'-H 2 -M "?"');
run_test(\&test_vcf_consensus,$opts,in=>'consensus3',out=>'consensus3.2.out',fa=>'consensus2.fa',args=>'-H 2 -M "?" -p xx_');
run_test(\&test_vcf_consensus,$opts,in=>'consensus4',out=>'consensus4.out',fa=>'consensus2.fa',args=>'-s -');
run_test(\&test_vcf_consensus,$opts,in=>'consensus5',out=>'consensus5.out',fa=>'consensus5.fa',args=>'--haplotype LA');
run_test(\&test_vcf_consensus,$opts,in=>'consensus6',out=>'consensus6.out',fa=>'consensus6.fa',args=>'-s -');
run_test(\&test_vcf_consensus,$opts,in=>'consensus7',out=>'consensus7a.out',fa=>'consensus7.fa',args=>'-H 2');
run_test(\&test_vcf_consensus,$opts,in=>'consensus7',out=>'consensus7a.out',fa=>'consensus7.fa',args=>'-H 4');
run_test(\&test_vcf_consensus,$opts,in=>'consensus7',out=>'consensus7b.out',fa=>'consensus7.fa',args=>'-H 2pIu');
run_test(\&test_vcf_consensus,$opts,in=>'consensus7',out=>'consensus7b.out',fa=>'consensus7.fa',args=>'-H 4pIu');
run_test(\&test_vcf_consensus,$opts,in=>'consensus7',out=>'consensus7c.out',fa=>'consensus7.fa',args=>'-H 1');
run_test(\&test_vcf_consensus,$opts,in=>'consensus7',out=>'consensus7c.out',fa=>'consensus7.fa',args=>'-H 3');
run_test(\&test_vcf_consensus,$opts,in=>'consensus7',out=>'consensus7d.out',fa=>'consensus7.fa',args=>'-H 1pIu');
run_test(\&test_vcf_consensus,$opts,in=>'consensus7',out=>'consensus7d.out',fa=>'consensus7.fa',args=>'-H 3pIu');
run_test(\&test_vcf_consensus,$opts,in=>'consensus8',out=>'consensus.8a.out',fa=>'consensus.fa',args=>'-s -');
run_test(\&test_vcf_consensus,$opts,in=>'consensus8',out=>'consensus.8b.out',fa=>'consensus.fa',args=>'-s - -a .');
run_test(\&test_vcf_consensus,$opts,in=>'consensus8',out=>'consensus.8c.out',fa=>'consensus.fa',args=>q[-s - -a . -i 'type="snp" || type="ref"']);
run_test(\&test_vcf_consensus,$opts,in=>'consensus8',out=>'consensus.8d.out',fa=>'consensus.fa',args=>q[-s - -a . -i 'ALT!="<DEL>"']);
run_test(\&test_vcf_consensus,$opts,in=>'consensus8',out=>'consensus.8e.out',fa=>'consensus.fa',args=>q[-s - -a . -e 'MinDP>15']);
run_test(\&test_vcf_consensus,$opts,in=>'consensus8',out=>'consensus.8f.out',fa=>'consensus.fa',args=>q[-s - -a . -e 'MinDP<15']);
run_test(\&test_vcf_consensus,$opts,in=>'consensus.9',out=>'consensus.9.1.out',fa=>'consensus.9.1.fa',args=>q[-H A]);
run_test(\&test_vcf_consensus,$opts,in=>'consensus.9',out=>'consensus.9.2.out',fa=>'consensus.9.2.fa',args=>q[-H A]);
run_test(\&test_vcf_consensus,$opts,in=>'consensus.10',out=>'consensus.9.1.out',fa=>'consensus.9.1.fa',args=>q[-H A]);
run_test(\&test_vcf_consensus,$opts,in=>'consensus.10',out=>'consensus.9.2.out',fa=>'consensus.9.2.fa',args=>q[-H A]);
run_test(\&test_vcf_consensus,$opts,in=>'consensus.11',out=>'consensus.11.1.out',fa=>'consensus.11.fa',args=>q[-s smpl]);
run_test(\&test_vcf_consensus,$opts,in=>'consensus.11',out=>'consensus.11.2.out',fa=>'consensus.11.fa',args=>q[-s smpl -a N]);
run_test(\&test_vcf_consensus,$opts,in=>'consensus.12',out=>'consensus.12.out',fa=>'consensus.12.fa',args=>'-s -');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.13',out=>'consensus.13.out',fa=>'consensus.13.fa',args=>'-s -');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.14',out=>'consensus.14.out',fa=>'consensus.14.fa',args=>'-s -');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.12',out=>'consensus.15.out',fa=>'consensus.12.fa',args=>'-s - --mark-del - --mark-ins uc --mark-snv uc');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.12',out=>'consensus.19.out',fa=>'consensus.12.fa',args=>'-s - --mark-del - --mark-ins + --mark-snv :');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.15',out=>'consensus.17.out',fa=>'consensus.15.fa',args=>'-H I --mark-ins lc --mark-snv lc');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.16',out=>'consensus.18.out',fa=>'consensus.fa',args=>'-s - -I');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.16',out=>'consensus.18.out',fa=>'consensus.fa',args=>'-H I');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.17',out=>'consensus17.1.out',fa=>'consensus2.fa',mask=>'consensus.17.bed',args=>'-s - -M N');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.18',out=>'consensus18.1.out',fa=>'consensus.18.fa',args=>'-s -');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.19',out=>'consensus19.1.out',fa=>'consensus.19.fa',args=>'-s - ');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.20',out=>'consensus20.1.out',fa=>'consensus.20.fa',args=>'-s -');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.20',out=>'consensus20.2.out',fa=>'consensus.20.fa',args=>'');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.20',out=>'consensus20.3.out',fa=>'consensus.20.fa',args=>'-M . -s b');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.20',out=>'consensus20.4.out',fa=>'consensus.20.fa',args=>'-M . -s a');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.21',out=>'consensus21.1.out',fa=>'consensus.21.fa',args=>'');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.22',out=>'consensus22.1.out',fa=>'consensus.22.fa',args=>'--regions-overlap 0 --mark-del .');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.22',out=>'consensus22.2.out',fa=>'consensus.22.fa',args=>'--regions-overlap 1 --mark-del .');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.22',out=>'consensus22.2.out',fa=>'consensus.22.fa',args=>'--regions-overlap 2 --mark-del .');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.22',out=>'consensus22.1.out',fa=>'consensus.22.fa',args=>'--regions-overlap 0');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.22',out=>'consensus22.3.out',fa=>'consensus.22.fa',args=>'--regions-overlap 1');
run_test(\&test_vcf_consensus,$opts,in=>'consensus.22',out=>'consensus22.3.out',fa=>'consensus.22.fa',args=>'--regions-overlap 2');
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.1.out',args=>q[-r17:100-150 -a -AD],test_list=>1);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.2.out',args=>q[-a DP,DV -r17:100-600 -a -AD]); # test files from samtools mpileup test suite
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.1)],out=>'mpileup/mpileup.3.out',args=>q[-B --ff 0x14 -r17:1050-1060 -a -AD]); # test file converted to vcf from samtools mpileup test suite
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.4.out',args=>q[-a DP,DPR,DV,DP4,INFO/DPR,SP,-AD -r17:100-600]); #test files from samtools mpileup test suite
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.5.out',args=>q[-a DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -r17:100-600]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.6.out',args=>q[-a DP,DV,-AD -r17:100-600 --gvcf 0,2,5]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.6.out',args=>q[-a DP,DV,-AD -r17:100-200,17:201-300,17:301-400,17:401-500,17:501-600 --gvcf 0,2,5]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.7.out',args=>q[-a -AD -r17:100-150 -s HG00101,HG00102]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.7.out',args=>q[-a -AD -r17:100-150 -S {PATH}/mplp.samples]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.8.out',args=>q[-a -AD -r17:100-150 -s ^HG00101,HG00102]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.8.out',args=>q[-a -AD -r17:100-150 -S ^{PATH}/mplp.samples]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.9.out',args=>q[-a -AD -t17:100-150 -S {PATH}/mplp.9.samples]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.1 mpileup.2 mpileup.3)],out=>'mpileup/mpileup.10.out',args=>q[-a -AD -t17:100-150 -G {PATH}/mplp.10.samples]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.3)],out=>'mpileup/mpileup.11.out',args=>q[-a -AD]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.3 mpileup.4)],out=>'mpileup/mpileup.11.out',args=>q[-a -AD -s HG00102]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.3 mpileup.4)],out=>'mpileup/mpileup.11.out',args=>q[-a -AD -s ^HG99999]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.3 mpileup.4)],out=>'mpileup/mpileup.11.out',args=>q[-a -AD -G {PATH}/mplp.11.rgs]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup.3 mpileup.4)],out=>'mpileup/mpileup.11.out',args=>q[-a -AD -G {PATH}/mplp.11.rgs]);
run_test(\&test_mpileup,$opts,in=>[qw(indel-AD.1)],out=>'mpileup/indel-AD.1.out',ref=>'indel-AD.1.fa',args=>q[-a AD]);
run_test(\&test_mpileup,$opts,in=>[qw(indel-AD.1)],out=>'mpileup/indel-AD.1cns.out',ref=>'indel-AD.1.fa',args=>q[-a AD --indels-cns]);
run_test(\&test_mpileup,$opts,in=>[qw(indel-AD.2)],out=>'mpileup/indel-AD.2.out',ref=>'indel-AD.2.fa',args=>q[-a AD -r 11:75]);
run_test(\&test_mpileup,$opts,in=>[qw(indel-AD.2)],out=>'mpileup/indel-AD.3.out',ref=>'indel-AD.2.fa',args=>q[-a AD -r 11:75 --ambig-reads incAD]);
run_test(\&test_mpileup,$opts,in=>[qw(indel-AD.2)],out=>'mpileup/indel-AD.4.out',ref=>'indel-AD.2.fa',args=>q[-a AD -r 11:75 --ambig-reads incAD0]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup-SCR)],out=>'mpileup/mpileup-SCR.out',ref=>'mpileup-SCR.fa',args=>q[-a -AD,INFO/SCR,FMT/SCR]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup-filter)],out=>'mpileup/mpileup-filter.1.out',ref=>'mpileup-SCR.fa',args=>q[-a -AD -t 1:100 --skip-all-set PAIRED,PROPER_PAIR,MREVERSE]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup-filter)],out=>'mpileup/mpileup-filter.1.out',ref=>'mpileup-SCR.fa',args=>q[-a -AD -t 1:100 --skip-any-set READ1]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup-filter)],out=>'mpileup/mpileup-filter.2.out',ref=>'mpileup-SCR.fa',args=>q[-a -AD -t 1:100 --skip-all-unset READ1]);
run_test(\&test_mpileup,$opts,in=>[qw(mpileup-filter)],out=>'mpileup/mpileup-filter.2.out',ref=>'mpileup-SCR.fa',args=>q[-a -AD -t 1:100 --skip-any-unset READ1]);
run_test(\&test_mpileup,$opts,in=>[qw(annot-NMBZ.1)],ref=>'annot-NMBZ.1.fa',out=>'mpileup/annot-NMBZ.1.1.out',args=>q[-a -AD,INFO/NMBZ -r chr19:69-99]);
run_test(\&test_mpileup,$opts,in=>[qw(annot-NMBZ.2)],ref=>'annot-NMBZ.2.fa',out=>'mpileup/annot-NMBZ.2.1.out',args=>q[-a -AD,INFO/NMBZ -r chr6:75]);
run_test(\&test_mpileup,$opts,in=>[qw(annot-NMBZ.3.1 annot-NMBZ.3.2)],ref=>'annot-NMBZ.3.fa',out=>'mpileup/annot-NMBZ.3.1.out',args=>q[-a -AD,INFO/NMBZ -r chr16:75]);
run_test(\&test_csq,$opts,in=>'csq',out=>'csq.1.out',cmd=>'-f {PATH}/csq.fa -g {PATH}/csq.gff3');
run_test(\&test_csq,$opts,in=>'csq',out=>'csq.1.out',cmd=>'-f {PATH}/csq.fa -g {PATH}/csq.chr.gff3');
run_test(\&test_csq,$opts,in=>'csq.2',out=>'csq.2.out',cmd=>'-f {PATH}/csq.fa -g {PATH}/csq.2.gff',tbcsq=>1);
run_test(\&test_csq,$opts,in=>'csq.2',out=>'csq.3.out',cmd=>'-f {PATH}/csq.fa -g {PATH}/csq.2.gff --ncsq 64',tbcsq=>1);
run_test(\&test_csq,$opts,in=>'csq.nchr',out=>'csq.chr.out',cmd=>'-f {PATH}/csq.nchr.fa -g {PATH}/csq.nchr.gff',tbcsq=>1);
run_test(\&test_csq,$opts,in=>'csq.nchr',out=>'csq.chr.out',cmd=>'-f {PATH}/csq.ychr.fa -g {PATH}/csq.nchr.gff',tbcsq=>1);
run_test(\&test_csq,$opts,in=>'csq.nchr',out=>'csq.chr.out',cmd=>'-f {PATH}/csq.nchr.fa -g {PATH}/csq.ychr.gff',tbcsq=>1);
run_test(\&test_csq,$opts,in=>'csq.nchr',out=>'csq.chr.out',cmd=>'-f {PATH}/csq.ychr.fa -g {PATH}/csq.ychr.gff',tbcsq=>1);
run_test(\&test_csq,$opts,in=>'csq.ychr',out=>'csq.chr.out',cmd=>'-f {PATH}/csq.ychr.fa -g {PATH}/csq.ychr.gff',tbcsq=>1);
run_test(\&test_csq,$opts,in=>'csq.ychr',out=>'csq.chr.out',cmd=>'-f {PATH}/csq.ychr.fa -g {PATH}/csq.nchr.gff',tbcsq=>1);
run_test(\&test_csq,$opts,in=>'csq.ychr',out=>'csq.chr.out',cmd=>'-f {PATH}/csq.nchr.fa -g {PATH}/csq.ychr.gff',tbcsq=>1);
run_test(\&test_csq,$opts,in=>'csq.ychr',out=>'csq.chr.out',cmd=>'-f {PATH}/csq.nchr.fa -g {PATH}/csq.nchr.gff',tbcsq=>1);
run_test(\&test_csq_real,$opts,in=>'csq');
run_test(\&test_roh,$opts,in=>'roh.1',out=>'roh.1.1.out',args=>q[-Or -G30 --AF-dflt 0.4]);
run_test(\&test_roh,$opts,in=>'roh.1',out=>'roh.1.1.out',args=>q[-Or -G30 --AF-file {PATH}/roh.1.tab.gz]);
run_test(\&test_roh,$opts,in=>'roh.1',out=>'roh.1.1.out',args=>q[-Or -G30 --AF-file {PATH}/roh.1.tab.gz --ignore-homref]);
run_test(\&test_roh,$opts,in=>'roh.1',out=>'roh.1.2.out',args=>q[    -G30 --AF-dflt 0.4 -r 1:100174876-100318245]);
run_test(\&test_roh,$opts,in=>'roh.1',out=>'roh.1.3.out',args=>q[    -G30 --AF-dflt 0.4 -r 1:100174876-100318245 --ignore-homref]);
run_test(\&test_roh,$opts,in=>'roh.1',out=>'roh.1.3.out',args=>q[    -G30 --AF-dflt 0.4 -r 1:100174876-100318245 --ignore-homref --include-noalt]);
run_test(\&test_roh,$opts,in=>'roh.1',out=>'roh.1.4.out',args=>q[    -G30 --AF-dflt 0.4 -r 1:100174876-100318245 --include-noalt]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.1.2.out',args=>q[-e 0 --no-HWE-prob]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.1.out',args=>q[-e 0]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.1.out',args=>q[-e 0 -u GT,GT]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.1.out',args=>q[-e 0 -u GT,PL]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.1.out',args=>q[-e 0 -u PL,GT]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.1.out',args=>q[-e 0 -u PL,PL]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.1.out',args=>q[-e 0 -p s1,s1]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.2',gts=>'gtcheck.1.gts',out=>'gtcheck.2.out',args=>q[-e 0]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.3',out=>'gtcheck.3.out',args=>q[-e 0 ]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.3',out=>'gtcheck.3.out',args=>q[-e 0 -p B,A,C,A,C,B,D,A,D,B,D,C,E,A,E,B,E,C,E,D]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.3',out=>'gtcheck.3.out',args=>q[-e 0 -P {PATH}/gtcheck.3.pairs]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.3',out=>'gtcheck.3.out',args=>q[-e 0 -u PL]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.3',out=>'gtcheck.9.out',args=>q[-e 0 --n-matches 4],sort=>1);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.3',out=>'gtcheck.4.out',args=>q[-e 0 -s qry:E,D,C]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.3',out=>'gtcheck.4.out',args=>q[-e 0 -s qry:E,D,C -u PL]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.3',out=>'gtcheck.5.out',args=>q[-e 0 -s qry:B -s gt:D]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.3',out=>'gtcheck.5.out',args=>q[-e 0 -s qry:B -s gt:D -u PL]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.3',out=>'gtcheck.6.out',args=>q[-e 0 -s qry:B -s gt:D,C]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.3',out=>'gtcheck.6.out',args=>q[-e 0 -p B,C,B,D]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.7.out',args=>q[-e 0 -u GT,GT -H]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.7.out',args=>q[-e 0 -u GT,PL -H]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.7.out',args=>q[-e 0 -u PL,GT -H]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.7.out',args=>q[-e 0 -u PL,PL -H]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.4',out=>'gtcheck.8.out',args=>q[-e 0 -P {PATH}/gtcheck.4.pairs --distinctive-sites 3]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.4',out=>'gtcheck.8.out',args=>q[-e 0 -P {PATH}/gtcheck.4.pairs --distinctive-sites 3,1]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.10.out',args=>q[-u GT -e 30]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.10.out',args=>q[-u GT -e 30 -p s1,s1]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.1',gts=>'gtcheck.1.gts',out=>'gtcheck.11.out',args=>q[-u GT -e 300]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.3',out=>'gtcheck.12.out',args=>q[-u PL -e 30]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.ntop',gts=>'gtcheck.ntop.gts',out=>'gtcheck.ntop.1.out',args=>q[]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.ntop',gts=>'gtcheck.ntop.gts',out=>'gtcheck.ntop.2.out',args=>q[--n-matches 2]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.5',gts=>'gtcheck.5.gts',out=>'gtcheck.5.1.out',args=>q[],grep=>'grep -v Time');
run_test(\&test_gtcheck,$opts,in=>'gtcheck.6',out=>'gtcheck.6.1.out',args=>q[-p A,B,B,C]);
run_test(\&test_gtcheck,$opts,in=>'gtcheck.3',out=>'gtcheck.3.1.out',args=>q[-t 11:33 -p A,D,A,E,D,E -u GT -e 10]);

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
        "   -f, --function LIST             Run only the listed tests (e.g. 'test_mpileup,split-vep')\n",
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
            'i|htsdir:s' => \$$opts{htsdir},
            'f|function:s' => \$$opts{function}
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

sub run_test
{
    my ($func,$opts,@args) = @_;
    if ( $$opts{function} )
    {
        if ( !exists($$opts{run_function}) )
        {
            $$opts{run_function} = { map {$_=>1} split(/,/,$$opts{function}) };
            use B qw(svref_2object);
        }
        my $name = svref_2object($func)->GV->NAME;
        my %args = @args;
        my $run  = 0;
        if ( exists($$opts{run_function}{$name}) ) { $run = 1; }
        if ( !$run )
        {
            for my $func (keys %{$$opts{run_function}})
            {
                if ( exists($args{cmd}) && $args{cmd}=~/$func/ ) { $run = 1; last; }
                if ( $name=~/$func/ ) { $run = 1; last; }
            }
        }
        if ( !$run ) { return; }
    }
    &$func($opts,@args);
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
        exec('bash', '-o','pipefail','-c', $cmd) or error("Cannot execute the command [bash -o pipefail -c $cmd]: $!");
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
        exec('bash', '-o','pipefail','-c', "($cmd) 2>$tmp.e >$tmp.o");
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
    if ( $ret && !$args{expected_failure} ) { failed($opts,$test,"Non-zero status $ret$err"); return; }
    if ( $args{expected_failure} )
    {
        if ( !$ret ) { failed($opts,$test,"Expected failure but the test returned $ret$err"); }
        else { passed($opts,$test,"ok, expected non-zero status"); }
        return;
    }
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
            my @diff = `diff $$opts{path}/$args{out} $$opts{path}/$args{out}.new`;
            for (my $i=0; $i<@diff; $i++) { $diff[$i] = "\t\t\t".$diff[$i]; }
            chomp($diff[-1]);
            failed($opts,$test,"The outputs differ:\n\t\t$$opts{path}/$args{out}\n\t\t$$opts{path}/$args{out}.new$err\n".join('',@diff));
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
    my ($opts,$test,$reason) = @_;
    $$opts{nok}++;
    if ( !defined $reason ) { $reason = 'ok'; }
    print ".. $reason\n\n";
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
    if ( $file=~m{/[^/]+} ) { cmd("mkdir -p $$opts{tmp}/$`"); }     # create a subdirectory if necessary
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
sub bgzip_index_bcf
{
    my ($opts,$file) = @_;
    if ( !-e "$$opts{tmp}/$file.bcf" or is_file_newer("$$opts{path}/$file.vcf","$$opts{tmp}/$file.bcf") )
    {
        cmd("$$opts{bin}/bcftools view -Ob $$opts{path}/$file.vcf -o $$opts{tmp}/$file.bcf");
    }
    if ( !-e "$$opts{tmp}/$file.bcf.csi" or is_file_newer("$$opts{tmp}/$file.bcf","$$opts{tmp}/$file.bcf.csi") )
    {
        cmd("$$opts{bin}/bcftools index -f $$opts{tmp}/$file.bcf");
    }
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
    if ( $args{args} eq '-n' )
    {
        test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools index $args{args} $$opts{tmp}/$args{in}.vcf.gz.csi");
    }
    unlink("$$opts{tmp}/$args{in}.vcf.gz.csi");

    cmd("$$opts{bin}/bcftools view -Ob $$opts{path}/$args{in}.vcf > $$opts{tmp}/$args{in}.bcf");
    cmd("$$opts{bin}/bcftools index -f $$opts{tmp}/$args{in}.bcf");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools index $args{args} $$opts{tmp}/$args{in}.bcf");
    if ( $args{args} eq '-n' )
    {
        test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools index $args{args} $$opts{tmp}/$args{in}.bcf.csi");
    }
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
    my @types = exists($args{types}) ? @{$args{types}} : (qw(vcf bcf));
    for my $type (@types)
    {
        my @files;
        if ( $args{noidx} )
        {
            for my $file (@{$args{in}})
            {
                if ( $type eq 'vcf' )
                {
                    push @files, "$$opts{path}/$file.vcf";
                }
                else
                {
                    cmd("$$opts{bin}/bcftools view --no-version -Ob -o $$opts{tmp}/$file.bcf $$opts{path}/$file.vcf");
                    push @files, "$$opts{tmp}/$file.bcf";
                }
            }
        }
        else
        {
            for my $file (@{$args{in}})
            {
                if ( $type eq 'vcf' )
                {
                    bgzip_tabix_vcf($opts,$file);
                    push @files, "$$opts{tmp}/$file.vcf.gz";
                }
                else
                {
                    cmd("$$opts{bin}/bcftools view --no-version -Ob -o $$opts{tmp}/$file.bcf $$opts{path}/$file.vcf && $$opts{bin}/bcftools index -f $$opts{tmp}/$file.bcf");
                    push @files, "$$opts{tmp}/$file.bcf";
                }
            }
        }
        my $args  = exists($args{args}) ? $args{args} : '';
        $args     =~ s/{PATH}/$$opts{path}/g;
        my $files = join(' ',@files);
        test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools merge --no-version $args $files", exp_fix=>1);
        test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools merge --no-version -Ob $args $files | $$opts{bin}/bcftools view --no-version | grep -v ^##bcftools_", exp_fix => 1);
    }
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
    $args{args} =~ s/{PATH}/$$opts{path}/g;
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
    my $cmd;
    if ( !exists($args{args}) )
    {
        $cmd = "$$opts{bin}/bcftools norm --no-version $$opts{tmp}/$args{in}.vcf.gz";
        if ( exists($args{fai}) ) { $cmd .= " -f $$opts{path}/$args{fai}.fa"; }
    }
    elsif ( ref($args{args}) eq 'ARRAY' )
    {
        my @cmd = ();
        for my $args (@{$args{args}})
        {
            $args =~ s/{PATH}/$$opts{path}/g;
            push @cmd, "$$opts{bin}/bcftools norm --no-version $args";
        }
        if ( exists($args{fai}) ) { $cmd[0] .= " -f $$opts{path}/$args{fai}.fa"; }
        $cmd[0] .= " $$opts{tmp}/$args{in}.vcf.gz";
        $cmd = join(' | ',@cmd);
    }
    else
    {
        $args{args} =~ s/{PATH}/$$opts{path}/g;
        $cmd = "$$opts{bin}/bcftools norm --no-version $args{args} $$opts{tmp}/$args{in}.vcf.gz";
        if ( exists($args{fai}) ) { $cmd .= " -f $$opts{path}/$args{fai}.fa"; }
    }
    test_cmd($opts,%args,cmd=>$cmd,exp_fix=>1);
    test_cmd($opts,%args,cmd=>"$cmd -Ou | $$opts{bin}/bcftools view | grep -v ^##bcftools_",exp_fix=>1);
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
sub gen_head_output
{
    my ($h, $n, $desc, $infile) = @_;

    my $expected = "";
    open my $in, '<', $infile or die "Couldn't open $infile: $!\n";
    while (<$in>) {
        my $counter = /^#/? \$h : \$n;
        next unless $$counter > 0;
        $expected .= $_;
        $$counter--;
    }
    close $in;
    return (out => "head.$desc.out", exp => $expected);
}
sub test_vcf_head
{
    my ($opts, %args) = @_;

    my $infile = "$$opts{path}/$args{in}";

    test_cmd($opts, %args, gen_head_output(1000, 0, "all", $infile),
             cmd => "$$opts{bin}/bcftools head $infile");
    test_cmd($opts, %args, gen_head_output(0, 0, "none", $infile),
             cmd => "$$opts{bin}/bcftools head -h 0 $infile");
    test_cmd($opts, %args, gen_head_output(1, 0, "one", $infile),
             cmd => "$$opts{bin}/bcftools head -h 1 $infile");
    test_cmd($opts, %args, gen_head_output(5, 0, "five", $infile),
             cmd => "$$opts{bin}/bcftools head -h 5 $infile");

    my $nh = $args{in_nheaders};  # Test exactly the number of headers
    test_cmd($opts, %args, gen_head_output($nh, 0, "exact", $infile),
             cmd => "$$opts{bin}/bcftools head -h $nh $infile");
    $nh++;  # Also test asking for one more line than there are headers
    test_cmd($opts, %args, gen_head_output($nh, 0, "toomany", $infile),
             cmd => "$$opts{bin}/bcftools head -h $nh $infile");

    test_cmd($opts, %args, gen_head_output(1000, 0, "alln0", $infile),
             cmd => "$$opts{bin}/bcftools head -n 0 $infile");
    test_cmd($opts, %args, gen_head_output(1000, 1, "onerec", $infile),
             cmd => "$$opts{bin}/bcftools head -n 1 $infile");
    test_cmd($opts, %args, gen_head_output(1000, 5, "fiverecs", $infile),
             cmd => "$$opts{bin}/bcftools head -n 5 $infile");
    test_cmd($opts, %args, gen_head_output(5, 5, "fiveboth", $infile),
             cmd => "$$opts{bin}/bcftools head -h 5 -n 5 < $infile");
}
sub test_vcf_head2
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools head $args{args} $$opts{tmp}/$args{in}.vcf.gz");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view --no-version -Ob $$opts{tmp}/$args{in}.vcf.gz | $$opts{bin}/bcftools head $args{args}");
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
    $args{args} =~ s/{PATH}/$$opts{path}/g;
    if ( exists($args{tab}) )
    {
        bgzip_tabix($opts,file=>$args{tab},suffix=>'tab',args=>'-s1 -b2 -e2');
        $annot_fname = "-a $$opts{tmp}/$args{tab}.tab.gz";
        $in_fname = "$$opts{path}/$args{in}.vcf";
        $hdr = (-e "$$opts{path}/$args{in}.hdr" && !($args{args}=~/-H/) && !($args{args}=~/--header-line\s/)) ? "-h $$opts{path}/$args{in}.hdr" : '';
    }
    elsif ( exists($args{bed}) )
    {
        bgzip_tabix($opts,file=>$args{bed},suffix=>'bed',args=>'-p bed');
        $annot_fname = "-a $$opts{tmp}/$args{bed}.bed.gz";
        $in_fname = "$$opts{path}/$args{in}.vcf";
        $hdr = (-e "$$opts{path}/$args{in}.hdr" && !($args{args}=~/-H/) && !($args{args}=~/--header-line\s/)) ? "-h $$opts{path}/$args{in}.hdr" : '';
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
    # Sadly, this does not work:
    #   $ENV{BCFTOOLS_PLUGINS} = "$$opts{bin}/plugins";
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
    test_cmd($opts,%args,cmd=>"export BCFTOOLS_PLUGINS=$$opts{bin}/plugins; $$opts{bin}/bcftools $args{cmd} $$opts{tmp}/$args{in}.vcf.gz $args{args} | grep -v ^##bcftools_");

    cmd("$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz > $$opts{tmp}/$args{in}.bcf");
    cmd("$$opts{bin}/bcftools index -f $$opts{tmp}/$args{in}.bcf");
    test_cmd($opts,%args,cmd=>"export BCFTOOLS_PLUGINS=$$opts{bin}/plugins; $$opts{bin}/bcftools $args{cmd} $$opts{tmp}/$args{in}.bcf $args{args} | grep -v ^##bcftools_", exp_fix=>1);
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
    bgzip_index_bcf($opts,$args{in});
    $args{args} =~ s/{PATH}/$$opts{path}/g;
    my $mask = $args{mask} ? "-m $$opts{path}/$args{mask}" : '';
    my $chain = $args{chain} ? "-c $$opts{tmp}/$args{chain}" : '';
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools consensus $$opts{tmp}/$args{in}.vcf.gz -f $$opts{path}/$args{fa} $args{args} $mask $chain");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools consensus $$opts{tmp}/$args{in}.bcf    -f $$opts{path}/$args{fa} $args{args} $mask $chain");
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

    my @hdr  = ();
    my $nhdr = 1 + int(rand($args{max_hdr_lines}));
    for (my $i=0; $i<$nhdr; $i++)
    {
        my $x = rand;
        push @hdr, "##INFO=<ID=XX$i,Number=1,Type=Integer,Description=\"Test Tag $x\">\n";
    }

    my @files = ();
    my $exp   = '';
    for (my $n=0; $n<$args{nfiles}; $n++)
    {
        my $nbdy = int(rand($args{max_body_lines}));
        my $file = "$$opts{tmp}/$args{name}.$n";
        push @files,$file;

        open(my $fh,'>',"$file.vcf") or error("$file.vcf: $!");
        print $fh "##fileformat=VCFv4.0\n";
        print $fh "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
        print $fh "##contig=<ID=1,length=62435964>\n";
        print $fh join('',@hdr);
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
    for my $fmt ('sam','bam','cram')
    {
        if ( $fmt eq 'sam' && ($args{args}=~/-r/ or $args{args}=~/-R/) ) { next; }
        my @files = ();
        for my $file (@{$args{in}})
        {
            if ( !-e "$$opts{path}/mpileup/$file.$fmt" ) { next; }
            push @files, "$$opts{path}/mpileup/$file.$fmt";
        }
        if ( !@files ) { next; }
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
    if ( $args{tbcsq} )
    {
        test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools csq $args{cmd} $$opts{path}/$args{in}.vcf | $$opts{bin}/bcftools query -f'[%TBCSQ\\n]' | perl -pe 's/[\\t,]/\\n/g' | sort");
    }
    else
    {
        test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools csq $args{cmd} $$opts{path}/$args{in}.vcf | $$opts{bin}/test/csq/sort-csq | $$opts{bin}/bcftools query -f'%POS\\t%REF\\t%ALT\\t%EXP\\n%POS\\t%REF\\t%ALT\\t%BCSQ\\n\\n'");
    }
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
            if ( $file=~/\.vcf$/ )
            {
                my $bname = $`;
                my $vcf   = "$dirname/$dir/$file";
                my $out   = "$args{in}/$dir/$bname.txt";
                my $outl  = "$args{in}/$dir/$bname.txt-l";
                my $cmd   = undef;
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
                my @cmd = grep { chomp } `cat $dirname/$dir/$file | grep -v ^#`;
                my $cmd = join(' ; ', @cmd);
                $cmd =~ s/{bin}/$$opts{bin}/g;
                test_cmd($opts,%args,out=>"$args{in}/$dir/$file.out",cmd=>"cd $dirname/$dir && $cmd");
            }
        }
        closedir($tmp);
    }
    closedir($dh);
}
sub test_plugin_split
{
    my ($opts,%args) = @_;
    if ( !$$opts{test_plugins} ) { return; }

    my ($package, $filename, $line, $test)=caller(0);
    $test =~ s/^.+:://;
    if ( !exists($args{args}) ) { $args{args} = ''; }
    $args{args} =~ s/{PATH}/$$opts{path}/g;

    cmd("export BCFTOOLS_PLUGINS=$$opts{bin}/plugins; $$opts{bin}/bcftools +split $$opts{path}/$args{in}.vcf -o $$opts{tmp}/$args{tmp} $args{args}");

    opendir(my $dh,"$$opts{tmp}/$args{tmp}") or failed($opts,$test,"Cannot read $$opts{tmp}/$args{tmp}: $!");
    my @files = sort grep { !(/^\./) } readdir($dh);
    closedir($dh) or failed($opts,$test,"Close failed: $$opts{tmp}/$args{tmp}");

    my $files = join(' ',@files);
    test_cmd($opts,%args,
        cmd=>
            "export BCFTOOLS_PLUGINS=$$opts{bin}/plugins; $$opts{bin}/bcftools +split $$opts{path}/$args{in}.vcf -o $$opts{tmp}/$args{tmp} $args{args} " .
            " && cd $$opts{tmp}/$args{tmp} " .
            " && for f in $files; do echo \$f; $$opts{bin}/bcftools query -l \$f; $$opts{bin}/bcftools view -H \$f; done"
        );
}
sub test_plugin_scatter
{
    my ($opts,%args) = @_;
    if ( !$$opts{test_plugins} ) { return; }

    my ($package, $filename, $line, $test)=caller(0);
    $test =~ s/^.+:://;
    if ( !exists($args{args}) ) { $args{args} = ''; }
    $args{args} =~ s/{PATH}/$$opts{path}/g;

    cmd("export BCFTOOLS_PLUGINS=$$opts{bin}/plugins; $$opts{bin}/bcftools +scatter $$opts{path}/$args{in}.vcf -o $$opts{tmp}/$args{tmp} $args{args}");

    opendir(my $dh,"$$opts{tmp}/$args{tmp}") or failed($opts,$test,"Cannot read $$opts{tmp}/$args{tmp}: $!");
    my @files = sort grep { !(/^\./) } readdir($dh);
    closedir($dh) or failed($opts,$test,"Close failed: $$opts{tmp}/$args{tmp}");

    my $files = join(' ',@files);
    test_cmd($opts,%args,cmd=>"export BCFTOOLS_PLUGINS=$$opts{bin}/plugins; $$opts{bin}/bcftools +scatter $$opts{path}/$args{in}.vcf -o $$opts{tmp}/$args{tmp} $args{args} && cd $$opts{tmp}/$args{tmp} && cat $files | grep -v ^##");
}
sub test_roh
{
    my ($opts,%args) = @_;
    $args{args} =~ s/{PATH}/$$opts{path}/g;
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools roh $$opts{path}/$args{in}.vcf.gz $args{args}| grep -v ^#");
}
sub test_gtcheck
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    $args{args} =~ s/{PATH}/$$opts{path}/g;
    if ( exists($args{gts}) ) { bgzip_tabix_vcf($opts,$args{gts}); }
    my $sort = exists($args{sort}) ? ' | sort' : '';
    my $gts  = exists($args{gts}) ? qq[-g $$opts{tmp}/$args{gts}.vcf.gz] : '';
    my $grep = exists($args{grep}) ? $args{grep} : "grep -v ^INFO $sort";
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools gtcheck $args{args} $$opts{tmp}/$args{in}.vcf.gz $gts | grep -v ^# | $grep");
}
sub test_vcf_merge_big
{
    my ($opts,%args) = @_;
    cmd("mkdir -p $$opts{tmp}/$args{in}");

    my $ref = 'A';
    my @alts = ();
    for (my $i=0; $i<$args{nalts}; $i++)
    {
        push @alts,'A'.('T' x ($i+1));
    }

    my $nsmpl = int($args{nsmpl}/($args{nfiles}-1));
    my $nalts = int($args{nalts}/($args{nfiles}-1));

    my $seed = srand(0);

    open(my $fh,'>',"$$opts{tmp}/$args{in}/list.txt") or error("$$opts{tmp}/$args{in}/list.txt: $!");
    my @files;
    for (my $i=0; $i<$args{nfiles}; $i++)
    {
        my @hdr = qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT);
        for (my $j=0; $j<$nsmpl; $j++) { push @hdr,'S'.($i*$nsmpl+$j); }
        my @alt = ();
        for (my $j=0; $j<$nalts; $j++)
        {
            my $k = int(rand(@alts));
            push @alt,$alts[$k];
        }
        my @out = (qw(1 3000 .),$ref,join(',',@alt),qw(. . . GT:PL));
        for (my $j=0; $j<$nsmpl; $j++)
        {
            my $igt = int(rand($nalts+1));
            my $jgt = int(rand($nalts+1));
            my @pl  = ();
            for (my $a=0; $a<1+@alt; $a++)
            {
                for (my $b=0; $b<=$a; $b++) { push @pl,int(rand(2147483000)); }
            }
            my $fmt = "$igt/$jgt:".join(',',@pl);
            push @out,$fmt;
        }

        my $file = "$$opts{tmp}/$args{in}/$i.vcf.gz";
        my $cmd = qq[$$opts{bin}/bcftools view -Oz -o $file];
        open(my $vcf,"| $cmd") or error("$cmd: $!");
        print $vcf qq[##fileformat=VCFv4.3\n];
        print $vcf qq[##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n];
        print $vcf qq[##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype Likelihood">\n];
        print $vcf qq[##contig=<ID=1,assembly=b37,length=249250621>\n];
        print $vcf qq[##reference=file:///ref.fa\n];
        print $vcf '#'.join("\t",@hdr)."\n";
        print $vcf join("\t",@out)."\n";
        close($vcf) or error("close failed: $cmd");
        cmd(qq[$$opts{bin}/bcftools index $file]);

        print $fh $file."\n";
    }
    close($fh) or error("close failed: $$opts{tmp}/$args{in}/list.txt");

    my $args  = exists($args{args}) ? $args{args} : '';
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools merge --no-version $args -l $$opts{tmp}/$args{in}/list.txt");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools merge --no-version $args -l $$opts{tmp}/$args{in}/list.txt -Ou | $$opts{bin}/bcftools view --no-version");
}


