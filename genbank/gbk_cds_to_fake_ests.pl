#!/usr/bin/env perl

=head1 NAME

gbk_cds_to_fake_ests.pl - extracts the coding portion of 'complete CDS' Genbank records and 
exports them (as if they were EST sequences) in FASTA format.  Useful for PASA runs.

=head1 SYNOPSIS

USAGE: gbk_cds_to_fake_ests.pl 
            --input_file=/path/to/some_file.gbk
            --output_file=/path/to/outfile.fsa
          [ --fabricate_ids=1 ]

=head1 OPTIONS

B<--input_file,-i>
    Input GBK file.  Can be a single or multi-entry file.

B<--output_file,-o>
    Output multi-FASTA file to be created.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

If we want to use gbk models as if they were ESTs we need to extract the coding portions of
the sequence out of the annotation and export them in FASTA format.

=head1  INPUT

The input Genbank file should have sections like this:

     CDS             join(324..390,444..1717)
                     /gene="SpoC1-C1C"
                     /protein_id="AAA33324.1"

It can have more attributes, they'll just be ignored.  The input file can also have multiple
genbank entries, each with multiple CDS.

=head1  OUTPUT

The output is a single FASTA file with multiple sequences within it.  The IDS of each entry in
the output file will be composed of the LOCUS of each entry with an incrementing id appended 
to the end.  For example, of this entry:

    LOCUS       ENU35731                1141 bp    DNA     linear   PLN 30-JUL-1996
     
     ....
     
         CDS             276..959
                     /gene="orlA"

The exported ID would be ENU35731_est_1

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Bio::SeqIO;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_file|i=s',
                          'output_file|o=s',
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

my $seqio = Bio::SeqIO->new( -file => $options{input_file}, -format => 'genbank' );

## entry is a Bio::Seq::RichSeq
while ( my $entry = $seqio->next_seq ) {
    
    if ( ! $entry->accession ) {
        die "got an entry without an accession - quitting";
    }
    
    my $parent_seq = Bio::PrimarySeq->new( -seq => $entry->seq );
    
    my $post_id_num = 1;
    
    ## look through the features (Bio::SeqFeature::Generic)
    for my $feature ( $entry->get_SeqFeatures ) {
        
        ## skip it unless it's a CDS
        next if $feature->primary_tag ne 'CDS';
       
        my $loc = $feature->location;
        my $seq = $parent_seq->subseq($loc);
        
        print $ofh ">", $entry->accession, "_est_$post_id_num\n$seq\n";
        
        $post_id_num++;
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
    #$options{fabricate_ids} = 0 unless defined $options{fabricate_ids};

}
