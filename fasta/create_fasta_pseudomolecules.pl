#!/usr/bin/env perl

=head1 NAME

create_fasta_pseudomolecules.pl - reads a collection of FASTA files and input options to create
one or more pseudomolecules.

=head1 SYNOPSIS

USAGE: create_fasta_pseudomolecules.pl 
            --input_fasta_list=/path/to/some.list
            --map_file=/path/to/some.map
            --output_file=/path/to/out.fsa
          [ --input_fasta_file=/path/to/some.fsa
            --unmapped_output=/path/to/unmapped.fsa
            --linker_sequence=NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN
          ]

=head1 OPTIONS

B<--input_fasta_list>
    A list file of FASTA sequences.  Each individual file can contain one or more sequences
    within it.  IDs for each sequence are taken to be the first string of characters after
    the header symbol > and up to the first whitespace.  These must be unique across your
    entire collection of input files.

B<--map_file> 
    This tab-delimited file species the order and orientation in which to piece together the
    pseudomolecule from the individual sequences.

B<--output_file>
    FASTA file that will be created by this script.  It will contain only the mapped sequences,
    with the unmapped optionally being written using the --unmapped_output option.

B<--input_fasta_file>
    Optional.  You may pass a single input file instead of a list of input files using
    this option.

B<--unmapped_output> 
    Optional.  If passed, the sequences present in the input files but not part of the mapping
    will be written to this file.

B<--linker_sequence>
    Optional.  This sequence will be inserted after each sequence stitched together into
    a pseudomolecule.  The default sequence (NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN) will be
    inserted if this isn't specified and contains 6-frame translational stop codons.  This
    can be an empty string (but that's a little on the crazy side)

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script creates output pseudomolecule FASTA files from your input sequences and a mapping
file.

=head1  INPUT

The FASTA files provided by either the --input_fast_list or --input_fasta_file options (or both)
can each contain multiple sequences.  IDs are accepted to be the first string after each > symbol
until the first whitespace.  These IDs must be unique across the entire sequence set (the script
will report a failure if they aren't).

The map files provide the order and orientation of the sequences within the pseudomolecules to be
created.  The tab-delimited format should have these columns:

    new_pseudomolecule_id       seq_id      direction

The first column allows you to create multiple pseudomolecules with a single run of this script.
Example rows would be:

    cdo.pseudomolecule.1.1      cdo.assembly.4.1        +
    cdo.pseudomolecule.1.1      cdo.assembly.2.1        -
    cdo.pseudomolecule.1.1      cdo.assembly.16.1       +
    cdo.pseudomolecule.2.1      cdo.assembly.18.1       +
    cdo.pseudomolecule.2.1      cdo.assembly.17.1       +

This will create an output file with two pseudomolecules.  The first (cdo.pseudomolecule.1.1) would
contain three assemblies.  First would be assembly 4 in forward orientation, then assembly 2 reverse
complemented and finally assembly 16 forward.

=head1  OUTPUT

The output is the result of the mapping in FASTA format.

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_fasta_list=s',
                          'map_file=s',
                          'output_file=s',
                          'input_fasta_file=s',
                          'unmapped_output=s',
                          'linker_sequence=s',
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

## need to read in the fasta files
my @fasta_infiles = parse_inputs();
_log("INFO: found " . scalar(@fasta_infiles) . " input FASTA files");

my $seqs = {};

## read in the sequence files
for my $file ( @fasta_infiles ) {
    _log("INFO: reading sequences in file $file");
    open(my $seq_fh, "<$file") || die "can't open sequence file $file: $!";
    
    my $current_id;
    ## using an array here is faster and uses less memory
    my @current_seq = ();
    
    while (my $line = <$seq_fh>) {
        next if $line =~ /^\s*$/;
        
        if ( $line =~ /^\>(\S+)/ ) {
            my $new_id = $1;
        
            ## current_id will be undefined if this is the first sequence in the file
            if ( defined $current_id ) {
                ## Bad Things if this sequence already exists
                if ( exists $$seqs{$current_id} ) {
                    die "duplicate input sequence ($current_id) found.  quitting.";
                }
            
                ## save the current sequence
                $$seqs{$current_id} = join('', @current_seq);
            
                ## remove whitespace from the sequence
                $$seqs{$current_id} =~ s|\s+||g;
                
                ## reset the array
                @current_seq = ();
            }
            
            ## set the new ID
            $current_id = $new_id;
            
        } else {
            push @current_seq, $line;
        }
    }
    
    ## remove whitespace from the last sequence
    $$seqs{$current_id} = join('', @current_seq);
    $$seqs{$current_id} =~ s|\s+||g;
}

open( my $outfh, ">$options{output_file}" ) || die "failed to create output file: $!";

open( my $mapfh, "<$options{map_file}" ) || die "failed to read map file: $!";

my $mapped_assembly_ids = {};
my $mapped_pseudomolecule_ids = {};

while ( my $line = <$mapfh> ) {
    chomp $line;
    next if $line =~ /^\s*$/;
    
    my @cols = split(/\t/, $line);
    if ( scalar @cols != 3 ) {
        _log("WARN: skipping the current map line");
        next;
    }
    
    ## if this pseudomolecule ID isn't mapped write the header line before continuing
    if ( ! exists $$mapped_pseudomolecule_ids{$cols[0]} ) {
        print $outfh ">$cols[0]\n";
        $$mapped_pseudomolecule_ids{$cols[0]}++;
    }
    
    ## make sure we know about the sequence in the map file
    if ( ! exists $$seqs{$cols[1]} ) {
        die("sequence ($cols[1]) defined in map file was not found among input FASTA seqs.");
    }
    
    ## add the current assembly with linker
    $$mapped_assembly_ids{$cols[1]}++;
    
    if ( $cols[2] eq '+' ) {
        _log("INFO: writing assembly $cols[1] to $cols[0] in forward orientation");
        print $outfh "$$seqs{$cols[1]}$options{linker_sequence}";
        
    } elsif ( $cols[2] eq '-' ) {
        my $seq = reverse $$seqs{$cols[1]};
           $seq =~ tr/ATGCatgc/TACGtacg/;
        
        _log("INFO: writing assembly $cols[1] to $cols[0] in reverse orientation");
        print $outfh "$seq$options{linker_sequence}";
    
    } else {
        die("unrecognized direction ($cols[2]) defined for sequence $cols[1]");
    }
}

## if requested dump any unmapped sequences
if ( defined $options{unmapped_output} ) {
    open(my $unmapped_fh, ">$options{unmapped_output}") || die "failed to create unmapped output file: $!";
    
    for my $seq_id ( keys %$seqs ) {
        if ( ! exists $$mapped_assembly_ids{$seq_id} ) {
            print $unmapped_fh ">$seq_id\n$$seqs{$seq_id}\n";
        }
    }
}


exit(0);

sub parse_inputs {
    my @files = ();

    if ( defined $options{input_fasta_list} ) {
        open(my $ilistfh, "<$options{input_fasta_list}") || die "failed to read input list: $!";

        while ( <$ilistfh> ) {
            chomp;
            next if /^\s*$/;

            push @files, $_;
        }
    }

    if ( defined $options{input_fasta_file} ) {
        push @files, $options{input_fasta_file};
    }
        return @files;
    }


sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( map_file output_file );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ## either input_fasta_list or input_fasta_file is required
    if ( ! defined $$options{input_fasta_list} && ! defined $$options{input_fasta_file} ) {
        die "you must pass either --input_fasta_list or --input_fasta_file";
    }
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## handle some defaults
    $$options{linker_sequence} = 'NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN' unless ( defined $$options{linker_sequence});
}
