#!/usr/bin/env perl

=head1 NAME

write_fasta_from_gff.pl - Extracts the protein or CDS sequences from a GFF3 file

=head1 SYNOPSIS

USAGE: write_fasta_from_gff.pl 
            --gff=/path/to/some_file.gff
            --output=/path/to/somefile.faa
          [ --type=protein|cds
            --fasta=/path/to/somefile.fna
          ]

=head1 OPTIONS

B<--gff,-g>
    Input file following the GFF3 specification.

B<--output,-o>
    Output FASTA file to be created.

B<--type,-t>
    Optional.  Should the output records be the CDS nucleotides or translated as
    proteins?  (values are 'cds' or 'protein', default=protein)

B<--fasta,-f>
    Optional.  If the FASTA entries for the underlying assemblies isn't 
    included in the GFF3 document you'll need to pass a FASTA file here.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This extracts the CDS entries from a GFF3 file, optionally translates them, and
exports a FASTA file of the CDS or resulting polypeptide sequences.

=head1  INPUT

The input file should follow the official GFF3 specification found here:

    http://www.sequenceontology.org/gff3.shtml

=head1  OUTPUT

A FASTA file with the CDS entries extracted and optionally translated to polypeptides.

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../lib";
use bioUtils;

my %options = ();
my $results = GetOptions (\%options, 
                          'gff|g=s',
                          'output|o=s',
                          'fasta|f=s',
                          'type|t=s',
                          'use_locus_tags!',
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

## structure like: 
#   $h{$mol_id}{$mRNA_id}, [ { min => N, max => N, strand => ?, phase => N } ];
my %coords = ();

my %molecules = ();

open(my $gffin_fh, $options{gff}) || die "failed to open input gff file: $options{gff}: $!";
open(my $fastaout_fh, ">$options{output}") || die "failed to create FASTA output file: $!";

my $in_fasta = 0;
my $current_seq_id;
my @current_seq_lines = ();

my %mRNA_locus_tags = ();

my %locus_lookup = ();

while ( my $line = <$gffin_fh> ) {
    chomp $line; 
    
    if ( $line =~ /^\#/ ) {
        next;
    }
    
    my @cols = split("\t", $line);
    
    if ( $line =~ /^\#\#FASTA/ ) {
        $in_fasta = 1;
        next;

    } elsif ( $in_fasta ) {
        if ( $line =~ /^>(\S+)/ ) {
            my $next_seq_id = $1;

            if ( $current_seq_id ) {            
                if ( exists $molecules{$current_seq_id} ) {
                    die "found duplicate molecule ($current_seq_id) in FASTA definitions";
                }            

                $molecules{$current_seq_id} = join( '', @current_seq_lines );
                $molecules{$current_seq_id} =~ s/\s//g;
                $molecules{$current_seq_id} = uc( $molecules{$current_seq_id} );
            }
            
            $current_seq_id = $next_seq_id;
            @current_seq_lines = ();
        } else {
            push @current_seq_lines, $line;
        }
    
    } else {
        next unless scalar @cols == 9;
    
        ## check for Parent definitions with others after 
        my $parent;
        my $locus_tag;
        my $feat_id;
        

        if ( $cols[2] eq 'mRNA') {
            if ($options{use_locus_tags}) {
                if ( $cols[8] =~ /locus_tag\=(.+?)\;/ ) {
                    $locus_tag = $1;
                } elsif ( $cols[8] =~ /locus_tag\=(.+)/ ) {
                    $locus_tag = $1;
                } else {
                    die "Failed to find locus tag on mRNA feature line: $line";
                }

                $cols[8] =~ /ID\=(.+?)\;/ || die "Failed to find ID attribute on the following mRNA line: $line";
                $feat_id = $1;

                $locus_lookup{$feat_id} = $locus_tag;
            }
        }

        if ( $cols[2] eq 'CDS' ) {

            if ( $cols[8] =~ /Parent\=(.+?)\;/ ) {
                $parent = $1;
            } elsif ( $cols[8] =~ /Parent\=(.+)/ ) {
                $parent = $1;
            } else {
                die "CDS features require Parent features to be defined";
            }

            push @{$coords{$cols[0]}{$parent}}, { min => $cols[3], max => $cols[4], 
                                        strand => $cols[6], phase => $cols[7] };

        }
    }
}

## do the last one and reset
if ($current_seq_id) {
    $molecules{$current_seq_id} = join( '', @current_seq_lines );
    $molecules{$current_seq_id} =~ s/\s//g;
    $molecules{$current_seq_id} = uc( $molecules{$current_seq_id} );
}

## sanity check that any sequence data was actually found
if (scalar keys %molecules == 0 && ! $options{fasta}) {
    die "ERROR: No sequence data found embedded in GFF3.  Perhaps you forgot to pass the path to a FASTA file via --fasta ? ";
}

$current_seq_id = '';
@current_seq_lines = ();

if ( $options{fasta} ) {
    open(my $fasta_fh, "<$options{fasta}") || die "failed to read FASTA file: $!";
    
    while (my $line = <$fasta_fh>) {
        if ( $line =~ /^>(\S+)/ ) {
            my $next_seq_id = $1;
            
            if ( $current_seq_id ) {
                if ( exists $molecules{$current_seq_id} ) {
                    die "found duplicate molecule ($current_seq_id) in FASTA definitions";
                }

                $molecules{$current_seq_id} = join( '', @current_seq_lines );
                $molecules{$current_seq_id} =~ s/\s//g;
                $molecules{$current_seq_id} = uc( $molecules{$current_seq_id} );
            }
            
            $current_seq_id = $next_seq_id;
            @current_seq_lines = ();
        } else {
            push @current_seq_lines, $line;
        }
    }
}

## do the last one and reset
$molecules{$current_seq_id} = join( '', @current_seq_lines );
$molecules{$current_seq_id} =~ s/\s//g;
$molecules{$current_seq_id} = uc( $molecules{$current_seq_id} );

for my $mol_id ( keys %coords ) {

    ## make sure this molecule exists in FASTA
    if ( ! exists $molecules{$mol_id} ) {
        die "ERROR: molecule ($mol_id) referenced in GFF but not found in FASTA\n";
    }

    for my $mrna_id ( keys %{$coords{$mol_id}} ) {
        my $complete_cds = '';
        my $strand;

        for my $cds ( sort {$$a{min} <=> $$b{min} } @{$coords{$mol_id}{$mrna_id}} ) {
            $complete_cds .= substr($molecules{$mol_id}, 
                                   $$cds{min} - 1,
                                   $$cds{max} - $$cds{min} + 1);
            
            $strand = $$cds{strand};
        }
        
        if ( $strand eq '-' ) {
            $complete_cds = reverse $complete_cds;
            $complete_cds =~ tr/ATGC/TACG/;
        }

        my $display_id = $mrna_id;

        if ($options{use_locus_tags}) {
            $display_id = $locus_lookup{$mrna_id};
        }
        
        if ($options{type} eq 'protein' ) {
            my $protein_seq = translate_sequence( $complete_cds, $mrna_id );
            print $fastaout_fh ">$display_id\n" . characters_per_line($protein_seq,60) . "\n";
        } else {
            print $fastaout_fh ">$display_id\n" . characters_per_line($complete_cds,60) . "\n";
        }
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
    my @required = qw( gff output );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    if ( defined $$options{type} ) {
        $$options{type} = lc $$options{type};
    
        if ( $$options{type} ne 'cds' && $$options{type} ne 'protein' ) {
            die "option --type must have a value of either cds or protein";
        }
    } else {
        $$options{type} = 'protein';
    }
}













