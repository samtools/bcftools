#!/usr/bin/env perl

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
test_vcf_check($opts,in=>'check',out=>'check.chk');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.out',args=>'-n =2');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.both.out',args=>'-n =2 -c both');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.any.out',args=>'-n =2 -c any');
test_vcf_isec($opts,in=>['isec.a','isec.b'],out=>'isec.ab.C.out',args=>'-C -c any');
test_vcf_isec2($opts,vcf_in=>['isec.a'],tab_in=>'isec',out=>'isec.tab.out',args=>'');
test_vcf_merge($opts,in=>['merge.a','merge.b','merge.c'],out=>'merge.abc.out');
test_vcf_query($opts,in=>'query',out=>'query.out',args=>q[-f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%DP4\\t%AN[\\t%GT\\t%TGT]\\n']);
test_vcf_norm($opts,in=>'norm',out=>'norm.out',fai=>'norm');
test_vcf_view($opts,in=>'view',out=>'view.1.out',args=>'-aRs NA00002 -v snps',reg=>'');
test_vcf_view($opts,in=>'view',out=>'view.2.out',args=>'-f PASS -k',reg=>'-r20,Y');
test_vcf_view($opts,in=>'view',out=>'view.3.out',args=>'-xs NA00003',reg=>'');
test_vcf_view($opts,in=>'view',out=>'view.4.out',args=>q[-i '%QUAL==999 && (FS<20 || FS>=41.02) && ICF>-0.1 && HWE*2>1.2'],reg=>'');
test_vcf_view($opts,in=>'view',out=>'view.5.out',args=>q[-p],reg=>'');
test_vcf_view($opts,in=>'view',out=>'view.6.out',args=>q[-P],reg=>'');
test_vcf_call($opts,in=>'mpileup',out=>'mpileup.1.out',args=>'-mv');
test_vcf_call_cAls($opts,in=>'mpileup',out=>'mpileup.cAls.out',tab=>'mpileup');
test_vcf_filter($opts,in=>'filter',out=>'filter.out',args=>'-mx -g2 -G2');
test_vcf_regions($opts,in=>'regions');
test_vcf_annotate($opts,in=>'annotate',tab=>'annotate',out=>'annotate.out',args=>'-c CHROM,POS,REF,ALT,ID,QUAL,INFO/T_INT,INFO/T_FLOAT,INDEL');
test_vcf_annotate($opts,in=>'annotate',tab=>'annotate2',out=>'annotate2.out',args=>'-c CHROM,FROM,TO,T_STR');
test_vcf_concat($opts,in=>['concat.1.a','concat.1.b'],out=>'concat.1.vcf.out',do_bcf=>0,args=>'');
test_vcf_concat($opts,in=>['concat.1.a','concat.1.b'],out=>'concat.1.bcf.out',do_bcf=>1,args=>'');
test_vcf_concat($opts,in=>['concat.2.a','concat.2.b'],out=>'concat.2.vcf.out',do_bcf=>0,args=>'-a');
test_vcf_concat($opts,in=>['concat.2.a','concat.2.b'],out=>'concat.2.bcf.out',do_bcf=>1,args=>'-a');
test_vcf_concat($opts,in=>['concat.3.a','concat.3.b','concat.3.c','concat.3.d'],out=>'concat.3.vcf.out',do_bcf=>0,args=>'-p');
test_vcf_concat($opts,in=>['concat.3.a','concat.3.b','concat.3.c','concat.3.d'],out=>'concat.3.bcf.out',do_bcf=>1,args=>'-p');

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
        "   -r, --redo-outputs              Recreate expected output files.\n",
        "   -t, --temp-dir <path>           When given, temporary files will not be removed.\n",
        "   -h, -?, --help                  This help message.\n",
        "\n";
    exit -1;
}
sub parse_params
{
    my $opts = { keep_files=>0, nok=>0, nfailed=>0 };
    my $help;
    Getopt::Long::Configure('bundling');
    my $ret = GetOptions (
            't|temp-dir:s' => \$$opts{keep_files}, 
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

    my ($ret,$out) = _cmd("$args{cmd} 2>&1");
    if ( $ret ) { failed($opts,$test); return; }
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
        cmd("cat $$opts{path}/$file | bgzip -c > $$opts{tmp}/$file.gz");
    }
    if ( $$opts{redo_outputs} or !-e "$$opts{tmp}/$file.gz.tbi" or is_file_newer("$$opts{tmp}/$file.gz","$$opts{tmp}/$file.gz.tbi") )
    {
        cmd("$$opts{bin}/bcftools tabix -f $args{args} $$opts{tmp}/$file.gz");
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
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools tabix $$opts{tmp}/$args{in}.vcf.gz $args{reg}");

    cmd("$$opts{bin}/bcftools view -Ob $$opts{tmp}/$args{in}.vcf.gz > $$opts{tmp}/$args{in}.bcf");
    cmd("$$opts{bin}/bcftools index $$opts{tmp}/$args{in}.bcf");
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view -H $$opts{tmp}/$args{in}.bcf $args{reg}");
}
sub test_vcf_check
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools stats -s - $$opts{tmp}/$args{in}.vcf.gz | grep -v '^# The command' | grep -v '^# This' | grep -v '^ID\t'");
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
    my $files = join(' ',@files);
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools merge $files | grep -v ^##bcftools_merge");
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
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools isec $args{args} -t $$opts{tmp}/$args{tab_in}.tab.gz $files | grep -v ^##bcftools_isec");
}
sub test_vcf_query
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools query $args{args} $$opts{tmp}/$args{in}.vcf.gz");
}
sub test_vcf_norm
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools norm -f $$opts{path}/$args{fai}.fa $$opts{tmp}/$args{in}.vcf.gz | grep -v ^##bcftools_norm");
}
sub test_vcf_view
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools view $args{args} $$opts{tmp}/$args{in}.vcf.gz $args{reg} | grep -v ^##bcftools_view");
}
sub test_vcf_call
{
    my ($opts,%args) = @_;
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools call $args{args} $$opts{path}/$args{in}.vcf | grep -v ^##bcftools_call");
}
sub test_vcf_call_cAls
{
    my ($opts,%args) = @_;
    bgzip_tabix($opts,file=>$args{tab},suffix=>'tab',args=>'-s1 -b2 -e2');
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools call -mA -C alleles -t $$opts{tmp}/$args{tab}.tab.gz $$opts{path}/$args{in}.vcf | grep -v ^##bcftools_call");
}
sub test_vcf_filter
{
    my ($opts,%args) = @_;
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools filter $args{args} $$opts{path}/$args{in}.vcf | grep -v ^##bcftools_filter");
}
sub test_vcf_regions
{
    my ($opts,%args) = @_;
    bgzip_tabix_vcf($opts,$args{in});

    # regions vs targets, holding tab in memory
    my $query = q[%CHROM %POS %REF,%ALT\n];
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -t $$opts{path}/$args{in}.tab $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -r $$opts{path}/$args{in}.tab $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');

    # regions vs targets, reading tabix-ed tab
    cmd(qq[cat $$opts{path}/$args{in}.tab | bgzip -c > $$opts{tmp}/$args{in}.tab.gz]);
    cmd(qq[$$opts{bin}/bcftools tabix -f -s1 -b2 -e3 $$opts{tmp}/$args{in}.tab.gz]);
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -t $$opts{tmp}/$args{in}.tab.gz $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -r $$opts{tmp}/$args{in}.tab.gz $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');

    # regions vs targets, holding bed in memory
    cmd(qq[cat $$opts{path}/$args{in}.tab | awk '{OFS="\\t"}{print \$1,\$2-1,\$3}' > $$opts{tmp}/$args{in}.bed]);
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -t $$opts{tmp}/$args{in}.bed $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -r $$opts{tmp}/$args{in}.bed $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');

    # regions vs targets, reading tabix-ed bed
    cmd(qq[cat $$opts{tmp}/$args{in}.bed | bgzip -c > $$opts{tmp}/$args{in}.bed.gz]);
    cmd(qq[$$opts{bin}/bcftools tabix -f -p bed $$opts{tmp}/$args{in}.bed.gz]);
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -t $$opts{tmp}/$args{in}.bed.gz $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');
    test_cmd($opts,cmd=>qq[$$opts{bin}/bcftools query -f'$query' -r $$opts{tmp}/$args{in}.bed.gz $$opts{tmp}/$args{in}.vcf.gz],out=>'regions.out');
}
sub test_usage
{
    my ($opts,%args) = @_;

    my $test = "test_usage";
    print "$test:\n";
    print "\t$args{cmd}\n";
    
    my $command = $args{cmd};
    my $commandpath = $$opts{bin}."/".$command;
    my ($ret,$out) = _cmd("$commandpath 2>&1");
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
	    $section =~ s/^[[:word:]]+[[:punct:][:space:]]*//;
	    $section =~ s/^[[:space:]]+//mg;
	    $section =~ s/^[[:punct:]]+.*?\n//msg;
	    @subcommands = ($section =~ m/^([[:word:]]+)[[:space:]].*/mg);
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
    my ($ret,$out) = _cmd("$commandpath $subcommand 2>&1");
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
    bgzip_tabix($opts,file=>$args{tab},suffix=>'tab',args=>'-s1 -b2 -e2');
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools annotate -a $$opts{tmp}/$args{tab}.tab.gz -h $$opts{path}/$args{in}.hdr $args{args} $$opts{path}/$args{in}.vcf | grep -v ^##bcftools_annotate");
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
            cmd("$$opts{bin}/bcftools index $$opts{tmp}/$file.bcf");
            $files .= " $$opts{tmp}/$file.bcf";
        }
        else
        {
            bgzip_tabix_vcf($opts,$file); 
            $files .= " $$opts{tmp}/$file.vcf.gz";
        }
    }
    test_cmd($opts,%args,cmd=>"$$opts{bin}/bcftools concat $args{args} $files | grep -v ^##bcftools_");
}
