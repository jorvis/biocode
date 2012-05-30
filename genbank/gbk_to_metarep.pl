#!/usr/bin/env perl

=head1 NAME

gbk_to_metarep.pl - Converts GBK files into a format for input into JCVI's METAREP tool

=head1 SYNOPSIS

USAGE: gbk_to_metarep.pl 
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

This script converts GBK files into a format for input into JCVI's METAREP tool here:

    http://jcvi.org/metarep/

This is useful if you've done your own assemblies and annotations from metagenomic datasets
and want to compare annotation profiles.

=head1  INPUT

Standard GBK format with annotations.

=head1  OUTPUT

    Column	Field ID	Description	Range/Type	Example	Multi-Valued
    1	peptide_id	Unique peptide ID	text	contig_6_1025_1_1	 no
    2	library_id	Library ID	text	Mock-stg-Illumina-PE-200904-SOAP_contig	no
    3	com_name	Functional description	text	sugar ABC transporter, periplasmic sugar-binding protein	yes
    4	com_name_src	Source of functional description assignment	text	RF|YP_167057.1|56696696|NC_003911 for hmm_id	yes
    5	go_id	Gene Ontology ID	text	GO:0009265	yes
    6	go_src	Source of Gene Ontology assignment	text	PF02511	yes
    7	ec_id	Enzyme Commission ID	text	2.1.1.148	yes
    8	ec_src	Source of Enzyme Commission ID	text	PRIAM	yes
    9	hmm_id	Hidden Markov Model hits (e.g. PFAM/TIGRFAM)	text	PF02511	yes
    10	blast_taxon	NCBI taxonomy ID	integer	246194	yes
    11	blast_evalue	BLAST E-Value	double	1.78E-20	no
    12	blast_pid	BLAST percent Identity	float	0.93	no
    13	blast_cov	BLAST sequence coverage of shortest sequence	float	0.82	no
    14	filter	Any filter tag (categorical variable)	text	repeat||duplicate	yes


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

print STDERR "INFO: processing file: $options{input_file}\n";

my $seqio = Bio::SeqIO->new( -file => $options{input_file} );

## entry is a Bio::Seq::RichSeq
while ( my $entry = $seqio->next_seq ) {
    
    my $seq_id = '';
    
    if ( $entry->accession_number ) {
        $seq_id = $entry->accession_number;
    } elsif ( $entry->id ) {
        $seq_id = $entry->id;
    } elsif ( $entry->display_id ) {
        $seq_id = $entry->display_id;
    } else {
        print STDERR "WARNING: no sequence identifier found\n";
    }
    
    ## look through the features (Bio::SeqFeature::Generic)
    for my $feature ( $entry->get_SeqFeatures ) {
        
        ## skip it unless it's a CDS
        next if $feature->primary_tag ne 'CDS';

        my $feature_id;
        my $product = '';
        my $ec_num = '';
        
        my @output = ('', $seq_id, '', '', '', '', '', '', '', '', '', '', '', '');
        
        
        if ( $feature->has_tag('protein_id') ) {
            $feature_id = ( $feature->get_tag_values('protein_id') )[0];
            $output[0] = $feature_id;
            
        } elsif ( $feature->has_tag('locus_tag') ) {
            $feature_id = ( $feature->get_tag_values('locus_tag') )[0];
            $output[0] = $feature_id;
            
        } else {
            print STDERR "WARNING: found a CDS with no locus_tag or protein_id at position " . $feature->start . ".\n";
        }
        
        if ( $feature->has_tag('product') ) {
            $product = ( $feature->get_tag_values('product') )[0];
            $output[2] = $product;
        }
        
        if ( $feature->has_tag('EC_number') ) {
            $ec_num = ( $feature->get_tag_values('EC_number') )[0];
            $output[4] = $ec_num;
        }
        
        print $ofh join("\t", @output), "\n";
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
}
