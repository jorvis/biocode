#!/usr/bin/env perl

=head1 NAME

filter_fasta_by_ids.pl - extract a subset of sequences from a multi-FASTA 
file based on an input list of IDs or a comma-separated ID or IDs (no spaces!)

=head1 SYNOPSIS

USAGE:  filter_fasta_by_ids.pl
              --id_list=/path/to/some/file.list OR --id=some_id1,some_id2
	          --fasta_file=/path/to/input/fasta.file
	          --output_file=/path/to/output.file
	      
This script prints a subset of FASTA sequence based on input IDs.  The inputs are a 
multi-FASTA file and one of the following:

	- a text file containing a list of IDs
	- a single ID or comma separated list of IDs on the command line
	
=head1 OPTIONS

B<--id_list, -l>
   A text file containing a list of IDs of interest
   
B<--fasta_file, -f>
   A multi-fasta file
   
B<--output_file, -o>
   An output file name for the multi-FASTA file
   
B<--id, -i>
   An ID or comma separated list of IDs (no spaces!)
   
B<--help, -h>
   This help message
   
=head1 EXAMPLES

    ./filter_fasta_by_ids.pl --id=1232,1435 --fasta_file=/path/to/file.fasta 
    --output_file=subset.fasta
 
=head1  CONTACT

    Anu Ganapathy
    alenas@gmail.com   
   
=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'id_list|l=s',
                          'id|i=s',
                          'fasta_file|f=s',
                          'output_file|o=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
my ($seq_id_ref, $fasta_file, $output_file) = &check_parameters(%options);

## store IDs
my %sequence_id = %$seq_id_ref;

## read through fasta
open(FILE, "$fasta_file");
open(OUTFILE, ">$output_file") || die "Could not open file $output_file: $!";

my $current_seq_id = "";

while(<FILE>){

    if($_ =~ /^>/){ ## defline
		my $id = $_;
		$id =~ s/^>//;
		($id) = split(/\s+/, $id);
		
		## if this ID is one we want ...
		if(exists($sequence_id{$id})){

	    	## reset loop variables for this ID
		    $current_seq_id = $id;
		    print OUTFILE $_;
		    next;
		}else{
		    $current_seq_id = "";
		}

    }elsif($_ =~ /^\s+$/){ ## blank line
		next;
    }else{ ## sequence line

        ## if we're currently watching an ID ...
		if($current_seq_id){ 
	    	## append its sequence
	  	  print OUTFILE $_;
		}
    }
}

close FILE;
close OUTFILE;

sub check_parameters {
	my (%options) = @_;
	
	## ids
	my %ids = ();
	
	## check if ids came from file
	if(defined $options{id_list}){
		## make sure file exists
		pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if (!-e $options{id_list});
	
		## store id list
		if(-e $options{id_list}){
			open(FILE, $options{id_list});
			while(my $line = <FILE>){
				chomp $line;
				$ids{$line} = "";
			}
			close FILE;
		}
	## if ids came from command line set
	}elsif(defined $options{id}){
		## id set could be csv
		my @id_set = split(',', $options{id});
		
		## initialize hash
		foreach (@id_set){
			$ids{$_} = "";
		}
	## no id list selected - die w/ usage info
	}else{
		die "Please enter an ID list via --id_list /path/to/file or --id someid1 \n";
	}
	
	## check fasta file
	die "Please enter a multi-FASTA file via --fasta_file " if(!defined $options{fasta_file} || !-e $options{fasta_file});
	
		
	## check output file
	die "Please enter an output file via --output_file " if(!defined $options{output_file});
	
	return (\%ids, $options{fasta_file}, $options{output_file});
	
}
