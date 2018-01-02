#!/usr/bin/env perl

=head1 NAME

tRNAScan_SE_to_gff3.pl - convert raw output of tRNAScan-SE to gff3

=head1 SYNOPSIS

USAGE: convert_tRNAScanSE_to_gff3.pl 
            --input=/path/to/some_file.out 

=head1 OPTIONS

B<--input,-i>
    The raw output from tRNAScan-SE:
    
    Sequence                        tRNA    Bounds  tRNA    Anti    Intron Bounds   Cove
    Name                    tRNA #  Begin   End     Type    Codon   Begin   End     Score
    --------                ------  ----    ------  ----    -----   -----   ----    ------
    tp.assembly.567468735.1         1       91820   91902   Tyr     GTA     91857   91866   66.58
    tp.assembly.567468735.1         2       171777  171849  Phe     GAA     0       0       70.28
    tp.assembly.567468735.1         3       172144  172215  His     GTG     0       0       64.04
    tp.assembly.567468735.1         4       852847  852919  Thr     AGT     0       0       75.69
    tp.assembly.567468735.1         5       877291  877362  Trp     CCA     0       0       68.97
    tp.assembly.567468735.1         6       1468229 1468300 Cys     GCA     0       0       72.10
    tp.assembly.567468735.1         7       2507459 2507530 Pro     AGG     0       0       62.33
    tp.assembly.567468735.1         8       2507198 2507127 Pro     CGG     0       0       65.73
    tp.assembly.567468735.1         9       2506317 2506246 Pro     TGG     0       0       66.60
    tp.assembly.567468735.1         10      2463785 2463713 Lys     TTT     0       0       79.47
    tp.assembly.567468735.1         11      2191149 2191069 Leu     CAG     0       0       57.47
    tp.assembly.567468735.1         12      1633307 1633237 Gly     CCC     0       0       65.52
    tp.assembly.567468735.1         13      1255051 1254968 Leu     CAA     0       0       60.46
    tp.assembly.567468735.1         14      251108  251037  Asp     GTC     0       0       59.48
    tp.assembly.567468735.1         15      250520  250449  Asp     GTC     0       0       59.48

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

# all output needs the gff header
print "##gff-version 3\n";

## globals
my $i=1;

## parse the file
foreach my $line (<$ifh>){
	my @cols = split /[\t]/, $line;
	chomp @cols;
	my $contig = $cols[0];

    if ($contig =~ /^(.+?)\s+$/) {
        $contig = $1;
    }

    ## skip the header lines
    next if $contig eq 'Sequence' || $contig eq 'Name' || $contig eq '--------';
    
	my $start = $cols[2];
	my $stop = $cols[3];
	my $trna = $cols[4];
	my $score = $cols[8];
	my $anticodon = $cols[5];
	
	## Check if tRNAScan-SE predicted a possible intron (value in column 6 > 0)
	if ($cols[6] > 0) {
   		if ($start < $stop){
			my $intron_start = $cols[6]-1; ## offset intron boundaries
			my $intron_end = $cols[7]+1;
			print "$contig\ttRNAScan-SE\tgene\t$start\t$stop\t$score\t+\t.\tID=tRNA-$trna$i\_gene\n";
			print "$contig\ttRNAScan-SE\ttRNA\t$start\t$stop\t$score\t+\t.\tID=tRNA-$trna$i\_tRNA;Name=tRNA-$trna;anticodon=$anticodon \n";
			print "$contig\ttRNAScan-SE\texon\t$start\t$intron_start\t$score\t+\t.\tID=tRNA-$trna$i\_exon;Note=contains predicted Intron\n";
			print "$contig\ttRNAScan-SE\texon\t$intron_end\t$stop\t$score\t+\t.\tID=tRNA-$trna$i\_exon;Parent=tRNA-$trna$i\_exon\n";
			$i++;
    	}else{
			my $intron_start = $cols[6]+1; ## offset intron boundaries
			my $intron_end = $cols[7]-1;
			print "$contig\ttRNAScan-SE\tgene\t$stop\t$start\t$score\t-\t.\tID=tRNA-$trna$i\_gene\n";
			print "$contig\ttRNAScan-SE\ttRNA\t$stop\t$start\t$score\t-\t.\tID=tRNA-$trna$i\_tRNA;Name=tRNA-$trna;anticodon=$anticodon \n";
			print "$contig\ttRNAScan-SE\texon\t$stop\t$intron_end\t$score\t-\t.\tID=tRNA-$trna$i\_exon;Note=contains predicted Intron\n";
			print "$contig\ttRNAScan-SE\texon\t$intron_start\t$start\t$score\t-\t.\tID=tRNA-$trna$i\_exon;Parent=tRNA-$trna$i\_exon\n";
			$i++;
		}
	} else{
		if ($start < $stop){
			print "$contig\ttRNAScan-SE\tgene\t$start\t$stop\t$score\t+\t.\tID=tRNA-$trna$i\_gene\n";
			print "$contig\ttRNAScan-SE\ttRNA\t$start\t$stop\t$score\t+\t.\tID=tRNA-$trna$i\_tRNA;Name=tRNA-$trna;anticodon=$anticodon \n";
			print "$contig\ttRNAScan-SE\texon\t$start\t$stop\t$score\t+\t.\tID=tRNA-$trna$i\_exon\n";
			$i++;
		}else{
			print "$contig\ttRNAScan-SE\tgene\t$stop\t$start\t$score\t-\t.\tID=tRNA-$trna$i\_gene\n";
			print "$contig\ttRNAScan-SE\ttRNA\t$stop\t$start\t$score\t-\t.\tID=tRNA-$trna$i\_tRNA;Name=tRNA-$trna;anticodon=$anticodon \n";
			print "$contig\ttRNAScan-SE\texon\t$stop\t$start\t$score\t-\t.\tID=tRNA-$trna$i\_exon\n";
			$i++;
		}
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

