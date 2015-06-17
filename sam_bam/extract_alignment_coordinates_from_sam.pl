#!/usr/bin/env perl
 
=head1 DESCRIPTION

Reads an alignment SAM file and extracts the coordinates for both query and subject.

The output format is tab-delimited with the following columns:

1.  Query ID
2.  Query start position
3.  Query end position
4.  Query alignment length
5.  Reference ID
2.  Reference start position
3.  Reference end position
4.  Reference alignment length

Original author:

Shaun Jackman
https://gist.github.com/sjackman

=cut

use strict;
use warnings;
 
use Getopt::Long;
use File::Spec;
 
 
my $script_name = (File::Spec->splitpath($0))[2];
my $usage = <<HEREDOC;
Usage: $script_name [SAM file]...

Options:

    --header    include a header line
    --help      show this help message
HEREDOC
 
my %options = ();
my $getopt_success = GetOptions(
    \%options,
    '--header',
    '--help',
);
 
if ($options{help}) { print $usage; exit 0; }
die $usage unless $getopt_success;
 
print join(
    "\t",
    'qname',
    'qstart',
    'qend',
    'qlen',
    'rname',
    'rstart',
    'rend',
    'rlen',
) . "\n" if $options{header};
 
while (<>) {
    my @fields = split /\t/;
    next unless @fields >= 10;
    my $qname = $fields[0];
    my $rname = $fields[2];
    my $rstart = $fields[3];
    my $cigar = $fields[5];
    my $qseq  = $fields[9];
#print sprintf("qseq len: %d\n", length($qseq));
 
    # Query
    my $qstart = 1;
    $_ = $cigar;
    s/^(\d+)[SH]/$qstart += $1/e;
    my $qlen = 0;
    $_ = $cigar;
    s/(\d+)[M=XI]/$qlen += $1/eg;
    my $qend = $qstart + $qlen - 1;
 
    # Reference
    my $rlen = 0;
    $_ = $cigar;
    s/(\d+)[M=XDN]/$rlen += $1/eg;
    my $rend = $rstart + $rlen - 1;
 
    if ($rname eq '*') {
        # case: query sequence is unmapped
        print join("\t", $qname, $qstart, $qend, $qlen) . "\n";
    } else {
        print join("\t", $qname, $qstart, $qend, $qlen, $rname, $rstart, $rend, $rlen) . "\n";
    }
}
