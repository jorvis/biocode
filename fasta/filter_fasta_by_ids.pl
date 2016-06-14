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

B<--mode, -m>
   Optional.  Specifies whether you want to include or exclude the IDs
   submitted from the output file.  Default=include.
   
B<--id, -i>
   An ID or comma separated list of IDs (no spaces!)
   
B<--help, -h>
   This help message
   
=head1 EXAMPLES

    ./filter_fasta_by_ids.pl --id=1232,1435 --fasta_file=/path/to/file.fasta 
    --output_file=subset.fasta
 
=head1  CONTACT

    Anu Ganapathy
    alenas at gmail

    highly modified by:
    Joshua Orvis
    jorvis at gmail
   
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
                          'mode|m=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
my ($seq_id_ref, $fasta_file, $output_file, $mode) = &check_parameters(%options);

## store IDs
my %sequence_id = %$seq_id_ref;

## read through fasta
open(FILE, "$fasta_file");
open(OUTFILE, ">$output_file") || die "Could not open file $output_file: $!";

my $keep_line = 0;

while(<FILE>){
    if($_ =~ /^>/){ ## defline
		my $id = $_;
		$id =~ s/^>//;
		($id) = split(/\s+/, $id);

        if ($mode eq 'include') {
            if (exists $sequence_id{$id}) {
                $keep_line = 1;
            } else {
                $keep_line = 0;
            }
        } else {
            if (exists $sequence_id{$id}) {
                $keep_line = 0;
            } else {
                $keep_line = 1;
            }
        }
    }

    print OUTFILE $_ if $keep_line;
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
        if (! -e $options{id_list}) {
            die "Error, file defined by --id_list not found\n";
        }
	
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

    if (defined $options{mode}) {
        if ($options{mode} ne 'include' && $options{mode} ne 'exclude') {
            die "--mode option value must be either 'include' or 'exclude'";
        }
    } else {
        $options{mode} = 'include';
    }
    
	## check fasta file
	die "Please enter a multi-FASTA file via --fasta_file " if(!defined $options{fasta_file} || !-e $options{fasta_file});
		
	## check output file
	die "Please enter an output file via --output_file " if(!defined $options{output_file});
	
	return (\%ids, $options{fasta_file}, $options{output_file}, $options{mode});
	
}
