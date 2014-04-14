#!/usr/bin/env perl

=head1 NAME

convert_RNAmmer_to_gff3.pl - convert bad gff2 format RNAmmer output to good GFF3

=head1 SYNOPSIS

USAGE: template.pl 
            --input=/path/to/some_file.out 

=head1 OPTIONS

B<--input,-i>
The raw output from RNAmmer:

##gff-version2
##source-version RNAmmer-1.2
##date 2014-03-26
##Type DNA
# seqname           source                      feature     start      end   score   +/-  frame  attribute
# ---------------------------------------------------------------------------------------------------------
tp.assembly.567497685.1 RNAmmer-1.2     rRNA    308022  312190  3003.7  -       .       28s_rRNA
tp.assembly.567468735.1 RNAmmer-1.2     rRNA    916705  922401  2985.9  -       .       28s_rRNA
tp.assembly.567492885.1 RNAmmer-1.2     rRNA    1034463 1034568 5.8     -       .       8s_rRNA
tp.assembly.567478335.1 RNAmmer-1.2     rRNA    762933  763050  69.8    -       .       8s_rRNA
tp.assembly.567478335.1 RNAmmer-1.2     rRNA    769881  769998  69.8    -       .       8s_rRNA
tp.assembly.567478335.1 RNAmmer-1.2     rRNA    775140  775257  69.8    -       .       8s_rRNA
tp.assembly.567492885.1 RNAmmer-1.2     rRNA    1036777 1036886 5.7     -       .       8s_rRNA
tp.assembly.567497685.1 RNAmmer-1.2     rRNA    312768  314511  1441.7  -       .       18s_rRNA
tp.assembly.567468735.1 RNAmmer-1.2     rRNA    922979  924722  1441.7  -       .       18s_rRNA
# ---------------------------------------------------------------------------------------------------------

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

File converter

=head1  INPUT

Input above.

=head1  OUTPUT
GFF3 to STDOUT

=head1  CONTACT

    Kyle Tretina
    kyletretina@gmail.com

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'input|i=s',
                          'log|l=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

## open the input file
my $ifh;
open($ifh, "<$options{input}") || die "can't open input file: $!";

## globals
my $i=1;
## parse the file
print "##gff-version 3\n";
foreach my $line (<$ifh>){
	if($line =~ /^#/){next;}
	my @cols = split /[\t]/, $line;
	chomp @cols;
	my $contig = $cols[0];
	my $start = $cols[3];
	my $stop = $cols[4];
	my $target = $cols[8];
	my $score = $cols[5];
	if ($start < $stop){
		print "$contig\tRNAmmer\tgene\t$start\t$stop\t$score\t+\t.\tID=$target\_$i\n";
		print "$contig\tRNAmmer\trRNA\t$start\t$stop\t$score\t+\t.\tID=$target\_$i\_rRNA;Parent=$target\_$i\n";
		print "$contig\tRNAmmer\texon\t$start\t$stop\t$score\t+\t.\tID=$target\_$i\_exon;Parent=$target\_$i\_rRNA\n";
		$i++;
	}else{
		print "$contig\tRNAmmer\tgene\t$start\t$stop\t$score\t-\t.\tID=$target\_$i\n";
                print "$contig\tRNAmmer\trRNA\t$start\t$stop\t$score\t-\t.\tID=$target\_$i\_rRNA;Parent=$target\_$i\n";
                print "$contig\tRNAmmer\texon\t$start\t$stop\t$score\t-\t.\tID=$target\_$i\_exon;Parent=$target\_$i\_rRNA\n";
		$i++;
	}
}

exit(0);


sub _log {
    my $msg = shift;
    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    my $options = shift;
    ## make sure required arguments were passed
    my @required = qw( input );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    ## handle some defaults
    $options{optional_argument2}   = 'foo'  unless ($options{optional_argument2});
}



