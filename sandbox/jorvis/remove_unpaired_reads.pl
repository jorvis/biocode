#!/usr/bin/env perl

=head1 DESCRIPTION

This was written quickly for a needed analysis, but could later be abstracted to be more
useful.  It should be used if you have two FASTA files of paired-end reads but where you're
not sure all reads are present.  This can happen with techniques like digital normalization.

Given any two FASTA files with naming convention ending with -1 or -2 to indicate direction
this will write new versions of both where each member is guaranteed to have a mate in the
other file.  There are no expectations of input order.

Memory usage is the amount of RAM required to store a hash of the IDs in the first file only.
This is at the expense of a bit of runtime, as the first file is parsed twice and the second
once during execution.

=cut

use strict;

my $file1 = shift;
my $file2 = shift;

my %file1_ids;

# first index all those in file one
open(my $ifh1, $file1) || die "Failed to open file1: $!";

while (my $line = <$ifh1>) {
    if ($line =~ /^\>(\S+)\-1/) {
        $file1_ids{$1} = 1;
    }
}
close($ifh1);

## read through file 2, only printing those where an entry was found in file1
open(my $ifh2, $file2) || die "Failed to open file2: $!";
open(my $ifh2_out, ">$file2.filtered") || die "Failed to create output file for file2: $!";

my $keep_seq = 0;

while (my $line = <$ifh2>) {
    if ($line =~ /^\>(\S+)\-2/) {
        my $seq_id = $1;

        if (exists $file1_ids{$seq_id}) {
            $keep_seq = 1;
            $file1_ids{$seq_id}++;
        } else {
            $keep_seq = 0;
        }
    }

    print $ifh2_out $line if $keep_seq;
}


## Finally, loop back through the first and export any pairs within it that were found
open($ifh1, $file1) || die "Failed to open file1, second loop: $!";
open(my $ifh1_out, ">$file1.filtered") || die "Failed to create output file for file1: $!";

while (my $line = <$ifh1>) {
    if ($line =~ /^\>(\S+)\-1/) {
        my $seq_id = $1;

        if ($file1_ids{$seq_id} == 2) {
            $keep_seq = 1;
        } elsif ($file1_ids{$seq_id} > 2 ) {
            print STDERR "ERROR: Found more than one instance of basename $seq_id in file2\n";
            $keep_seq = 0;
        } else {
            $keep_seq = 0;
        }
    }

    print $ifh1_out $line if $keep_seq;
}
close($ifh1);
