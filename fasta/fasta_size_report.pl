#!/usr/bin/env perl

=head1 NAME

fasta_size_report.pl - print a file containing fasta ID, sequence size, the defline, 
and source file (NOTE: Only verified on a single file input)

=head1 SYNOPSIS

USAGE: fasta_size_report.pl 
            --file=/path/to/some_file.fsa --file=/path/to/some/other/file.fsa
                       OR
            --list=/path/to/some_file.list
                       OR
            --dir=/path/to/somedir
          [ --summary=1 ]

=head1 OPTIONS

B<--file,-f>
    a single or multiple fasta files of interest

B<--list,-l>
    an input list containing a set of fasta files

B<--dir,-d>
    a directory containing .fsa, .fasta, .fa files 

B<--summary,-s>
    optional.  If passed, also exports summary lines (each begins with a # symbol)

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script takes an input fasta file, a file containing a list of fasta files, or a 
directory of fasta files and prints a synopsis file containing the ID, the size of the 
sequence, the original defline, and the originating file.

=head1  INPUT

There are several potential inputs:

- A single fasta file

- A list of fasta files:
    
    fasta.list
      /path/to/A.fasta
      /path/to/B.fasta

- A directory of fasta files with the extension .fsa, .fa, or .fasta

=head1  OUTPUT

The output is a tabular list of:

ID    SIZE    DEFLINE    FILE
1     100     >1 blah    /path/to/A.fasta

Optionally, if the -s option is set to 1, also exports summary lines, each beginning
with a comment.

=head1  CONTACT

    Anu Ganapathy
    alenas@gmail.com

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;




my %options = ();
my $results = GetOptions (\%options, 
                          'file|f=s@',
                          'list|l=s',
                          'dir|d=s',
                          'summary|s',
                          'log|l=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

## make sure everything passed was peachy and retrieve files of interest
my @files = &check_parameters(\%options);

my $summary_record_count = 0;
my $summary_size = 0;

## print mapping for each file
foreach my $file (@files){
    &print_size($file);
}

if ( $options{summary} ) {
    print "\nSummary:\n" .
          "\tRecord count   : $summary_record_count\n" .
          "\tSum of residues: $summary_size\n";
}

exit(0);


sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @input_types = qw( file dir list );
    my @inputs = ();

    ## file
    if  ( defined $options->{file} ) {
        @inputs = @{$options->{file}};
    }

    ## list
    if ( defined $options->{list} ) {
	open FILE, $options->{list} || die "Could not open input list $options->{list}: $!";
	while (my $file = <FILE>){
	    chomp;
	    push @inputs, $file;
	}
	close FILE;
    }

    ## dir (accept all 3 ways specified in documentation)
    if( defined $options->{dir} ) {
        die "Invalid path give for directory input: $options->{dir} " if (! -d $options->{dir} );
        my $dir = $options->{dir};
        my @files = <$dir/*.fsa>;
        @inputs = (@inputs, @files);
        @files = <$dir/*.fa>;
        @inputs = (@inputs, @files);
        @files = <$dir/*.fasta>;
        @inputs = (@inputs, @files);
    }

    ## check to make sure we have files 
    die "Script requires at least one fasta file input via the -f or -l options" if (scalar @inputs == 0);

    ## check to make sure they all exist
    foreach my $file (@inputs){
        die "Error with file $file.  Check if it exists and is readable" if (!-e $file || !-r $file);
    }

    return @inputs;
}

sub print_size {
    my ($file) = @_;

    open (my $FIN, $file) || die "Unable to read fasta file ($file): $!";

    ## initialize fields of interest
    my $len     = 0; 
    my $defline = "";
    my $id      = "";
    
    ## go through file
    while (my $line = <$FIN>) {
	chomp($line);
	
	## header line
	if ($line =~ /^>/) {

	    ## if we processed a line already
	    if ($len > 0) {
		print "$id\t$len\t$defline\t$file\n";
        $summary_record_count++;
        $summary_size += $len;
        $len = 0;
	    }

	    $defline = $line;

	    $line =~ s/^>//;
	    ($line) = split(/\s+/, $line);
	    $id = $line;
	}
	## sequence line
	else {
	    $line =~ s/\s+//g;  ## remove spaces
	    $len += length($line);
	}
    }
    print "$id\t$len\t$defline\t$file\n";
    $summary_record_count++;
    $summary_size += $len;


    close($FIN);

}
