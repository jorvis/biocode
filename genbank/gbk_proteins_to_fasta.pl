#!/usr/bin/env perl

=head1 NAME

gbk_proteins_to_fasta.pl - Extracts a set of proteins from a GBK record.  


=head1 SYNOPSIS

USAGE: gbk_proteins_to_fasta.pl
            --input_file=/path/to/some_file.gbk
            --output_file=/path/to/outfile.fsa
          [ --fabricate_ids=1 ]

=head1 OPTIONS

B<--input_file,-i>
    Input GBK file.  Can be a single or multi-entry file.

B<--output_file,-o>
    Output multi-FASTA file to be created.

B<--fabricate_ids,-f>
    optional.  This script first considers 'locus_tag' values for IDs, then 'protein_id'.  If 
    they're not present it will report failure unless this option is passed, causing IDs to be 
    made up from the feature coordinates like: CDS_656_1204.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Several GBK -> FASTA converters will skip the feature table and just output the source (genomic)
sequence at the end of the file.  This script reads all the CDS entries in a GBK file and
writes a multi-FASTA file of their sequences.

=head1  INPUT

The input Genbank file should have sections like this:

     CDS             1306..1674
                     /locus_tag="RF_p02"
                     /translation="MSPNPDFITIVKIANYFNCAVDQVVGRRKFLPSINLIVSFNNPD
                     LNDINSNLCNFLKAKLSQDNISPYLLSKNIGFSKKIIHCFLKANSPYKMLSTNVIIAL
                     ADYFNVSVDDMIERYPTTKQ"

It can have more attributes - I've just removed those that are ignored.  If there is no
locus_tag or protein_id attribute the script will fail unless the --fabricate_ids option is passed with
a value of 1.  If so, IDs will be derived from the coordinates.  The sequence above, for
example, would get an ID of CDS_1306_1674.

It should be noted that the sequence is extracted from the /translation attribute, rather
than piecing the CDS together from the coordinates.

=head1  OUTPUT

The output is a single FASTA file containing the polypeptide sequences.

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Bio::SeqIO;
use FindBin;
use lib "$FindBin::Bin/../lib";
use bioUtils;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_file|i=s',
                          'output_file|o=s',
                          'fabricate_ids|f=s',
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

open(my $ofh, ">$options{output_file}") || die "failed to open output file $options{output_file}: $!";

my $seqio = Bio::SeqIO->new( -file => $options{input_file} );

## entry is a Bio::Seq::RichSeq
while ( my $entry = $seqio->next_seq ) {
    
    ## look through the features (Bio::SeqFeature::Generic)
    for my $feature ( $entry->get_SeqFeatures ) {
        
        ## skip it unless it's a CDS
        next if $feature->primary_tag ne 'CDS';

        my $feature_id;
        my $product = '';
        
        if ( $feature->has_tag('locus_tag') ) {
            $feature_id = ( $feature->get_tag_values('locus_tag') )[0];
        } elsif ( $feature->has_tag('protein_id') ) {
            $feature_id = ( $feature->get_tag_values('protein_id') )[0];
        } elsif ( $options{fabricate_ids} ) {
            ## no locus_tag  was found - we'll need to derive an ID.
            $feature_id = 'CDS_' . $feature->start . '_' . $feature->end;
        } else {
            die "found a CDS with no locus_tag at position " . $feature->start . ".  read about the --fabricate_ids option\n";
        }
        
        if ( $feature->has_tag('product') ) {
            $product = ( $feature->get_tag_values('product') )[0];
        }
        
        my $translation = '';
        if ( $feature->has_tag('translation') ) {
            $translation = ($feature->get_tag_values('translation'))[0];

        ## no translation tag provided - let bioperl try        
        } else {
            $translation = $feature->seq->translate->seq();
        }
        
        if ( length $translation < 1 ) {
            die "failed to get a translation for CDS $feature_id\n";
        }
        
        print $ofh ">$feature_id $product\n" . characters_per_line($translation, 60) . "\n";
    }
}


exit(0);


sub _log {
    my $msg = shift;

    print STDERR "$msg\n";
    print $logfh "$msg\n" if $logfh;
}



sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( input_file output_file );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## handle some defaults
    $options{fabricate_ids} = 0 unless defined $options{fabricate_ids};
    if ( $options{fabricate_ids} != 0 && $options{fabricate_ids} != 1 ) {
        die "--fabricate_ids option value must be either 0 or 1\n";
    }
}
