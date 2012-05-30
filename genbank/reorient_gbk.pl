#!/usr/bin/env perl

=head1 NAME

reorient_gbk.pl - Performs operations on a GBK file such as adjusting to a new origin or 
    reverse/inverting the sequence.

=head1 SYNOPSIS

USAGE: reorient_gbk.pl 
            --input_file=/path/to/some_file.gbk
            --output_file=/path/to/some_file.gbk
          [ --new_origin=1500
          ]

=head1 OPTIONS

B<--input_file,-i>
    The path to a single GenBank flat file that can contain one or more entries. 

B<--output_file,-o>
    Path of output file to be created.

B<--new_origin,-n>
    Optional.  Use this to adjust the origin of the molecule.  It should be a positive number between
    1 and the end of the molecule, and it will be recircularized at that point.  Use 0 or omit this
    option to leave origin as is.  Use of this option assumes a circular molecule.

B<--reverse_molecule,-r>
    Optional.  NOT YET IMPLEMENTED
    
B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script can be used to change the origin of features within any GBK file, and assumes a circular
molecule.

=head1  INPUT

Any Genbank flatfile containing one or more records.  Only the coordinates and sequence of 
each entry are manipulated.

=head1  OUTPUT

The resulting GBK file should contain all of the header, features and sequence information as the
original, but with the coordinates of the features adjusted as well as the sequence reordered to
reflect these new coordinates.

=head1  LIMITATIONS

The script doesn't allow adjustment of the origin to a location WITHIN an annotated feature.  It 
will check for this and die a horrible death should it occur.

Also, currently the BioPerl methods used here allow for the recoordination of the annotated
features but doesn't resort them to their proper place in the document.  For most software this
shouldn't matter, but this should be fixed when possible.

Doesn't currently work with multi-part / fragmented features (such as spliced genes).

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

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
my $seqout = Bio::SeqIO->new( -file => ">$options{output_file}", -format => 'Genbank' );

## entry is a Bio::Seq::RichSeq
while ( my $entry = $seqin->next_seq ) {

    ## is the user completely confused?  maybe just too much free time on their hands?
    #   if '1' just skip any manipulations and spit the unmodified record back out
    if ( $options{new_origin} != 1 ) {

        my $molecule_accession = $entry->accession;
        my $molecule_length = length($entry->seq);

        _log("INFO: parsing annotation for accession: $molecule_accession");
        _log("INFO: molecule length ($molecule_length)");

        ## look through the features (Bio::SeqFeature::Generic)
        for my $feature ( $entry->get_SeqFeatures ) {

            ## the source feature will throw this off
            next if $feature->primary_tag eq 'source';

            ## features have start, end, strand
            if ( $options{new_origin} ) {

                ## origin change can't be within this feature
                if ( overlap( $feature, $options{new_origin}) ) {
                    _logdie( "ERROR: new origin defined falls within an existing feature." );
                }

                if ( $feature->start >= $options{new_origin} ) {
                    $feature->start( $feature->start - $options{new_origin} + 1 );
                    $feature->end( $feature->end - $options{new_origin} + 1 );   

                } elsif ( $feature->start < $options{new_origin} ) {
                    $feature->start( $molecule_length -  $options{new_origin} + $feature->start + 1 );
                    $feature->end( $molecule_length -  $options{new_origin} + $feature->end + 1 );
                }
            }
        }

        ## fix the sequence
        my $new_residues = $entry->subseq( $options{new_origin}, $molecule_length ) . 
                           $entry->subseq( 1, $options{new_origin} - 1 );

        $entry->seq( $new_residues );
    }

    ## now write the entry back out
    $seqout->write_seq($entry);
}


exit(0);

sub overlap {
    my ($feature, $point) = @_;
    
    if ( $feature->strand == 1 ) {
    
        if ( $point >= $feature->start && $point <= $feature->end  ) {
            return 1;
        }
    
    } elsif ( $feature->strand == -1 ) {
        if ( $point >= $feature->end && $point <= $feature->start  ) {
            return 1;
        }    
    }
    
    return 0;
}


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
    
    ## warn them of their strange mistake
    if ( $options{new_origin} == 1 ) {
        print STDERR "WARNING: setting --new_origin=1 makes the script not actually do anything.\n";
    }
    
    
    if ( $options{reverse_molecule} ) {
        print STDERR "\n\nERROR: --reverse_molecule option not yet implemented\n\n";
        exit(1);
    }
}
