#!/usr/bin/env perl

=head1 NAME

reorient_fasta.pl - Adjusts the sequence in a FASTA file to a new origin.

=head1 SYNOPSIS

USAGE: reorient_fasta.pl
            --input_file=/path/to/some_file.fasta
            --output_file=/path/to/some_file.fasta
          [ --new_origin=1500
          ]

=head1 OPTIONS

B<--input_file,-i>
    The path to a single FASTA flat file which should contain a single sequence.

B<--output_file,-o>
    Path of output file to be created.

B<--new_origin,-n>
    Optional.  Use this to adjust the origin of the molecule.  It should be a positive number between
    1 and the end of the molecule, and it will be recircularized at that point.  Use 0 or omit this
    option to leave origin as is.  Use of this option assumes a circular molecule.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Adjusts the sequence in a FASTA file to a new origin.

=head1  INPUT

One FASTA file enters.

=head1  OUTPUT

One FASTA file leaves.

=head1  CONTACT

    Jonathan Crabtree
    jcrabtree@som.umaryland.edu

=cut

use strict;
use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_file|i=s',
                          'output_file|o=s',
                          'new_origin|n=i',
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

my $seqin = Bio::SeqIO->new( -file => $options{input_file} );
my $seqout = Bio::SeqIO->new( -file => ">$options{output_file}", -format => 'FASTA' );

## entry is a Bio::Seq::RichSeq
while ( my $entry = $seqin->next_seq ) {

    my $molecule_length = length($entry->seq);

    ## fix the sequence
    my $new_residues = $entry->subseq( $options{new_origin}, $molecule_length ) . 
                       $entry->subseq( 1, $options{new_origin} - 1 );

    $entry->seq( $new_residues );

    ## now write the entry back out
    $seqout->write_seq($entry);
}


exit(0);


sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}

sub _logdie {
    my $msg = shift;

    print STDERR "\n$msg\n\n";
    _log($msg);
    exit(1);
}

sub fatal {
    my $msg = shift;
    print STDERR "\n$msg\n\n";
    exit(1);
}

sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( input_file output_file );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            fatal("--$option is a required option");
        }
    }
    
    ## input_file or input_list is required
    if ( ! exists $$options{input_file} && ! exists $$options{input_list} ) {
        fatal("must define input with either --input_file or --input_list");
    }
    
    ## handle some defaults
    $options{new_origin} = 0  unless (defined $options{new_origin});
}
