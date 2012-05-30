#!/usr/bin/env perl

=head1 NAME

write_protein_fasta_from_gff.pl - Extracts the protein sequences from a GFF3 file

=head1 SYNOPSIS

USAGE: write_protein_fasta_from_gff.pl 
            --gff=/path/to/some_file.gff
            --output=/path/to/somefile.faa
          [ --fasta=/path/to/somefile.fna
          ]

=head1 OPTIONS

B<--gff,-g>
    Input file following the GFF3 specification.

B<--output,-o>
    Output FASTA file to be created.

B<--fasta,-f>
    Optional.  If the FASTA file for the underlying assemblies isn't 
    included in the document you'll need to pass the file here.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This extracts the CDS entries from a GFF3 file, translates them, and
exports a FASTA file of the resulting polypeptide sequences.

=head1  INPUT

The input file should follow the official GFF3 specification found here:

    http://www.sequenceontology.org/gff3.shtml

=head1  OUTPUT

A FASTA file with the CDS entries extracted and translated to polypeptides.

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'gff|g=s',
                          'output|o=s',
                          'fasta|f=s',
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

open(my $gffin_fh, $options{gff}) || die "failed to get gff input option";
open(my $fastaout_fh, ">$options{output}") || die "failed to create FASTA output file: $!";

my $in_fasta = 0;
my $current_seq_id;
my @current_seq_lines = ();

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
$molecules{$current_seq_id} = join( '', @current_seq_lines );
$molecules{$current_seq_id} =~ s/\s//g;
$molecules{$current_seq_id} = uc( $molecules{$current_seq_id} );

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
            $complete_cds = substr($molecules{$mol_id}, 
                                   $$cds{min} - 1,
                                   $$cds{max} - $$cds{min} + 1);
            
            $strand = $$cds{strand};
        }
        
        if ( $strand eq '-' ) {
            $complete_cds = reverse $complete_cds;
            $complete_cds =~ tr/ATGC/TACG/;
        }
        
        my $protein_seq = translate_sequence( $complete_cds, $mrna_id );
        
        print $fastaout_fh ">$mrna_id\n$protein_seq\n";
    }
}




exit(0);

sub translate_sequence {
    my ($dna, $label) = @_;
    
    my $polypeptide_seq = '';
    
    for(my $i=0; $i<(length($dna) - 2); $i+=3) {
        $polypeptide_seq .= codon2aa( substr($dna,$i,3) );
    }
    
    ## report if this had an internal stop
    if ( $polypeptide_seq =~ /\*.*?.$/ ) {
        _log("ERROR: internal stops detected in feature with feature_id: $label: ($polypeptide_seq)");
    }
    
    return $polypeptide_seq;
}

sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if (exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{
        return 'X';
    }
}

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
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## handle some defaults
    #$options{optional_argument2}   = 'foo'  unless ($options{optional_argument2});
}
