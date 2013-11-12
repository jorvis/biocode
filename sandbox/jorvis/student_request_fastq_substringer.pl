#!/usr/bin/env perl

=head1 DESCRIPTION

Problem definition from C. Adler:

It is a single FASTQ file.  Its read length is 100bp, and it was a paired-end RNA-seq 
run. I think I might have explained poorly, though.

I'd like to create a script to get 4 files from this 1 file - 1 where the 100bp is 
chopped to the first 33bp, 1 where the 100bp is chopped to the first 50bp, 1 where 
the 100bp is chopped to the first 75bp, and the 1 for the 100bp original read length 
(which I think could just be the file as it stands). 

=cut

use warnings;
use strict;

my $infile_path = shift || die "Pass me an input FASTQ file!";

open(my $ifh, $infile_path) || die "failed to open input file: $!";

my @entry;
my @desired_sizes = (33, 50, 75, 100);
my $file_basename;

if ( $infile_path =~ /(.+)\.fastq$/i || $infile_path =~ /(.+)\.fq$/i) {
    $file_basename = $1;
} else {
    die "Expected the input file to end in either .fastq or .fq";
}

my %outfiles = ();
for my $size ( @desired_sizes ) {
    my $output_path = $file_basename . ".$size.fastq";
    open($outfiles{$size}, ">$output_path") || die "Failed to create file: $output_path because: $!";
}

while (my $line = <$ifh>) {
    chomp $line;
    
    if ( scalar(@entry) == 4 ) {
        process_entry(\@entry);
        @entry = ();
    }

    push @entry, $line;
}

process_entry(\@entry);


sub process_entry {
    my $lines = shift;

    if ($$lines[0] !~ /^\@/) {
        die "ERROR: attempted to process an entry which did not start with this symbol \@";
    }

    for my $size (@desired_sizes) {
        print {$outfiles{$size}} "$$lines[0]\n";

        if ( length($$lines[1]) < $size ) {
            die "ERROR: Sequence size requested ($size) is greater than actual sequence length: $$lines[1]";
        }

        print {$outfiles{$size}} substr($$lines[1], 0, $size) . "\n";
        print {$outfiles{$size}} "$$lines[2]\n";
        print {$outfiles{$size}} substr($$lines[3], 0, $size) . "\n";
    }
}
