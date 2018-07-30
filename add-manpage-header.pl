#!/usr/bin/env perl
#
# Author: petr.danecek@sanger
#

use strict;
use warnings;
use Carp;

my $opts = parse_params();
add_header($opts);

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print 
        "Usage: add-manpage-header.pl [OPTIONS] bcftools-man.html > bcftools.html\n",
        "Options:\n",
        "   -h, -?, --help            This help message\n",
        "\n";
    exit -1;
}
sub parse_params
{
    my $opts = {};
    while (defined(my $arg=shift(@ARGV)))
    {
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        if ( -e $arg && !exists($$opts{manpage}) ) { $$opts{manpage}=$arg; next }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    if ( !exists($$opts{manpage}) ) { error(); }
    return $opts;
}


sub add_header
{
    my ($opts) = @_;
    open(my $fh,'<',$$opts{manpage}) or error("$$opts{manpage}: $!");
    while (my $line=<$fh>)
    {
        if ( $line=~/<body>/ )
        {
            print $`.$&;
            print q[
                <div style="background-color:#FCF8E3;text-align:center;padding:1em;">
                This documentation refers to the <b>latest development version</b> of BCFtools which can
                be downloaded from github, see <a href="https://samtools.github.io/bcftools/">instructions</a>.
                <P>
                Please refer to <a href="http://www.htslib.org/doc/bcftools.html">htslib.org</a> for
                documentation for the latest <b>versioned release</b>.
                </div>
            ];
            print $';
            next;
        }
        print $line;
    }
    close($fh) or error("close failed: $$opts{manpage}");
}

