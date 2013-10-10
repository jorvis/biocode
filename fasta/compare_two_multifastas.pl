#!/usr/bin/env perl

=head1 DESCRIPTION

    This script reads two multi-fasta files and generates a report on their content, also 
    checking to make sure all sequences found in one are also found in the other.  Sequence
    comparisons are case-sensitive.  ID differences are not considered, but are reported in
    the case of errors.
    
    RAM is cheap, and this quick and dirty script can use lots of it.  :)
    
    This script also writes a mapping of the matching IDs to /tmp/id.map

=head1 USAGE

    ./compare_two_multifastas.pl file1.fsa file2.fsa

=head1 CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use Data::Dumper;

if ( @ARGV != 2 ) {
    print "\nUSAGE: $0 file1.fsa file2.fsa\n\n";
    exit(1);
}

my ($file1, $file2) = @ARGV;

my $differences_found = 0;
my $multiple_match_count = 0;

my $DEBUG_MODE = 0;
open(my $debug_fh, ">/tmp/compare.debug.out") || die "couldn't write debugging file: $!";


## structure like:
#   $h->{id} = {
#        seq => 'ATG.....',
#        ## number of times this identical sequence is found in the other file
#        match_count => 2,
#    }
my $seqs1 = read_seqs_from_file( $file1 );
print $debug_fh Dumper($seqs1) if $DEBUG_MODE;

my $seqs2 = read_seqs_from_file( $file2 );
print $debug_fh Dumper($seqs2) if $DEBUG_MODE;

print "\nsequence counts:\n";
print scalar(keys %$seqs1), "\t$file1\n";
print scalar(keys %$seqs2), "\t$file2\n";

if ( scalar(keys %$seqs1) != scalar(keys %$seqs2) ) {
    $differences_found++;
}

print "\nchecking for unique sequences in $file1\n";
my @uniq_in_f1 = ();

open(my $idmapfh, ">/tmp/id.map") || die "couldn't create id.map output file: $!";

for my $s1_id ( keys %$seqs1 ) {
    
    my @matches = ();
    my $s1_len = length $$seqs1{$s1_id}{seq};
    
    
    for my $s2_id ( keys %$seqs2 ) {
    
        ## first just check and make sure they're the same length
        if ( length $$seqs2{$s2_id}{seq} == $s1_len ) {
        
            ## now check for identical residues
            if ( $$seqs1{$s1_id}{seq} eq $$seqs2{$s2_id}{seq} ) {
                push @matches, $s2_id;
                print $idmapfh "$s1_id\t$s2_id\n";
            } else {
                ## can print this for debugging
                #print "\twarning: $s1_id matches $s2_id in length but not sequence.  looking for others\n";
            }
        }
    }

    if ( scalar @matches > 1 ) {
        $multiple_match_count++;
    }

    if (! scalar @matches ) {
        #print "\tentry $s1_id unique to file1\n";
        push @uniq_in_f1, $s1_id;
        $differences_found++;
    }
}

print "checking for unique sequences in $file2\n";
my @uniq_in_f2 = ();

for my $s2_id ( keys %$seqs2 ) {
    
    my @matches = ();
    my $s2_len = length $$seqs2{$s2_id}{seq};
    
    for my $s1_id ( keys %$seqs1 ) {
        ## first just check and make sure they're the same length
        if ( length $$seqs1{$s1_id}{seq} == $s2_len ) {
        
            ## now check for identical residues
            if ( $$seqs2{$s2_id}{seq} eq $$seqs1{$s1_id}{seq} ) {
                push @matches, $s1_id;
            } else {
                ## can print this for debugging
                #print "\twarning: $s2_id matches $s1_id in length but not sequence.  looking for others\n";
            }
        }
    }
    
    if (! scalar @matches ) {
        #print "\tentry $s2_id unique to file2\n";
        push @uniq_in_f2, $s2_id;
        $differences_found++;
    }
}

print "\n$multiple_match_count entries in file 1 had multiple matches in file 2\n";

if ( $differences_found ) {
    print "\nthere were $differences_found sequence differences found between the files\n";
    exit(1);

} else {
    print "\nthere were no sequence differences found between the files\n\n";
    exit(0);
}

if ( scalar @uniq_in_f1 ) {
    print "\nThe following entries were unique to $file1:\n";

    for (@uniq_in_f1) {
        print "\t$_\n";
    }
}

if ( scalar @uniq_in_f2 ) {
    print "\nThe following entries were unique to $file2:\n";

    for (@uniq_in_f2) {
        print "\t$_\n";
    }
}

sub read_seqs_from_file {
    my ($file) = @_;
    
    my $data = {};
    
    open(my $ifh, "<$file") || die "can't read input file $file: $!";
    
    my $current_seq_id = '';
    my @current_seq_parts = ();
    
    while ( my $line = <$ifh> ) {
        chomp $line;
        
        if ( $line =~ /\>(\S+)/ ) {
            
            my $new_seq_id = $1;
            
            if ( $current_seq_id ) {
                $$data{$current_seq_id} = new_seq( join('', @current_seq_parts) );
            }
            
            $current_seq_id = $new_seq_id;
            @current_seq_parts = ();
            
        } else {
            push @current_seq_parts, uc($line);
        }
    }
    
    $$data{$current_seq_id} = new_seq( join('', @current_seq_parts) );
    
    return $data;
}

sub new_seq {
    my ($residues) = @_;
    
    $residues =~ s|\s||g;
    $residues =~ s|\*||g;
    
    return {
        seq => $residues,
        match_count => 0,
    };
    
}
