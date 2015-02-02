#!/usr/bin/env perl

=head1 NAME

annotation_statistics.pl - produces general annotation statistics from a variety 
of annotation sources, though initially written for files in GFF3 format describing 
eukaryotic annotation.

=head1 SYNOPSIS

USAGE: annotation_statistics.pl 
            --input_file=/path/to/some_file.gff3
            --format=gff3
          [ --input_list=/path/to/some.list
            --source=maker
          ]

=head1 OPTIONS

B<--input_file,-i>
    If you have a single input file, use this option.

B<--input_list,-l>
    Optional.  If you have multiple input files, create a list file where each line is
    a path to one of them.

B<--format,-f>
    The only currently-supported format is gff3.  Feel free to pass other options though
    to see a nice little error message.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Reads a file and prints selected statistics on the annotation features found within
it.  Currently only GFF3 is supported.  See the sections below for details.

=head1  INPUT

Because so many of the statistics rely on reading the underlying sequence, it must be stored
for access while this script runs.  Depending on how much is referenced in your input files,
this script could use considerable memory.

The GFF3 specification should be followed and this script very much assumes (improperly 
perhaps) that the 'gene' features in the GFF file are defined before any subfeatures of 
each gene.  

=head1  OUTPUT

Exported statistics are:

    Molecule count
	Gene count
	Single-exon gene count
	Mean gene length
	Mean CDS length
	Mean introns per gene
	Mean intron length
	Mean intergenic region length
	%GC of protein-coding regions

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_file|i=s',
                          'input_list|l=s',
                          'format|f=s',
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

my $files = get_files_from_opts( \%options );

=head1 ANNOTATION DATA STRUCTURE

Regardless of the format, input files are parsed and returned as a data structure like
this one.  All coordinate are 0-based, interbase.

$annot = 
{
    $molecule_id = {
        seq => $nucleotide_sequence,
        genes => [
            {
                id      => ?,
                fmin    => ?,  ## these coordinates are for the overall gene range
                fmax    => ?,
                strand  => ?,
                exons  => [   ## this is an array for just the coding portions of the gene
                    { fmin => ?, fmax => ? }    ## strand would be redundant
                ],
                CDS => [
                    { fmin => ?, fmax => ? }    ## strand would be redundant
                ]
            }, ...
        ]
        
    }, ...
}

=cut

my $annot = parse_annotation_files( $files, $options{format} );

#use Data::Dumper;
#print Dumper( $annot );

my $stats = {
    molecule_count => 0,
    gene_count => 0,
    gene_sum_length => 0,  ## used to calculate mean length
    mean_gene_length => 0,
    joined_exon_count => 0,
    joined_exon_sum_length => 0,
    joined_exon_gc_count => 0,
    joined_exon_gc_perc => 0,
    mean_joined_exon_length => 0,
    single_exon_gene_count => 0,
    intron_count => 0,
    mean_intron_count_per_gene => 0,
    jointed_intron_sum_length => 0,
    intergenic_sum_length => 0,
};

for my $mol_id ( keys %$annot ) {
    $$stats{molecule_count}++;

    ## do we have the sequence for this molecule so that we can export statistics
    #   based on sequence content?
    my $have_mol_seq = 0;
    if ( exists $$annot{$mol_id}{seq} && length($$annot{$mol_id}{seq}) > 0 ) {
        $have_mol_seq = 1;
    } else {
        print STDERR "WARN: we don't appear to have a sequence for molecule ($mol_id)\n";
    }

    ## these are used for intergenic sequence tracking
    my $last_gene_pos = -1;

    if ( $$annot{$mol_id}{genes} == undef ) {
        $$annot{$mol_id}{genes} = [];
    }

    ## iterate through the genes, in order
    for my $gene ( sort { $$a{fmin}<=>$$b{fmin} } @{$$annot{$mol_id}{genes}} ) {
        $$stats{gene_count}++;
        $$stats{gene_sum_length} += abs( $$gene{fmax} - $$gene{fmin} );
        
        my $this_gene_joined_exon_length = 0;
        
        $$stats{single_exon_gene_count}++ if scalar @{$$gene{exons}} == 1;
        $$stats{intron_count} += scalar(@{$$gene{exons}}) - 1 if scalar(@{$$gene{exons}}) > 0;
        
        for my $exon ( @{$$gene{exons}} ) {
            $this_gene_joined_exon_length += abs( $$exon{fmax} - $$exon{fmin} );
            my $exon_seq = substr($$annot{$mol_id}{seq}, $$exon{fmin}, $$exon{fmax} - $$exon{fmin} );
            $$stats{joined_exon_gc_count} += ( $exon_seq =~ tr/GCgc// );
        }
        
        $$stats{joined_exon_sum_length} += $this_gene_joined_exon_length;
        $$stats{joined_exon_count}++ if $this_gene_joined_exon_length;
        $$stats{jointed_intron_sum_length} += calculate_intron_length( $gene );
        
        if ( $last_gene_pos != -1 ) {
            $$stats{intergenic_sum_length} += $$gene{fmin} - $last_gene_pos;
        }
        
        $last_gene_pos = $$gene{fmax};
    }
}

$$stats{mean_gene_length} = sprintf("%.1f", $$stats{gene_sum_length}/$$stats{gene_count} );
$$stats{mean_joined_exon_length} = sprintf("%.1f", $$stats{joined_exon_sum_length}/$$stats{gene_count} );
$$stats{mean_intron_count_per_gene} = sprintf("%.1f", $$stats{intron_count}/$$stats{gene_count} );
$$stats{mean_intron_length} = sprintf("%.1f", $$stats{jointed_intron_sum_length}/$$stats{intron_count} );
$$stats{mean_intergenic_length} = sprintf("%.1f", $$stats{intergenic_sum_length}/($$stats{gene_count} - 1) );
$$stats{joined_exon_gc_perc} = sprintf("%.1f", ($$stats{joined_exon_gc_count}/$$stats{joined_exon_sum_length}) * 100 );

print "\nMolecule count: $$stats{molecule_count}\n";
print "Gene count: $$stats{gene_count}\n";
print "Single-exon gene count: $$stats{single_exon_gene_count}\n";
print "Mean gene length: $$stats{mean_gene_length}\n";
print "Mean CDS length : $$stats{mean_joined_exon_length}\n";
print "Mean introns per gene: $$stats{mean_intron_count_per_gene}\n";
print "Mean intron length: $$stats{mean_intron_length}\n";
print "Mean intergenic region length: $$stats{mean_intergenic_length}\n";
print "\%GC of protein-coding regions: $$stats{joined_exon_gc_perc}\n";

exit(0);

sub calculate_intron_length {
    my $gene = shift;
    
    return 0 unless scalar(@{$$gene{exons}}) > 1;
    
    my $base_count = 0;
    my $last_pos = -1;
    
    for my $cds ( sort { $$a{fmin}<=>$$b{fmin} } @{$$gene{CDS}} ) {
        if ( $last_pos != -1 ) {
            $base_count += $$cds{fmin} - $last_pos;
        }
        
        $last_pos = $$cds{fmax};
    }
    
    return $base_count;
}


sub parse_gff3_file {
    my ($infile, $ann) = @_;
    
    open(my $ifh, $infile) || die "can't read input GFF3 file ($infile): $!";
    
    my $in_fasta_section = 0;
    my $current_fasta_seq_id;
    my @current_fasta_seq = ();
    
    my $current_gene;     ## one of the elements of 'genes' in $annotation as the file is parsed
    my $current_gene_mol;
    
    while (my $line = <$ifh>) {
        chomp $line;
        
        if ( $line =~ /^\#\#FASTA/ ) {
            $in_fasta_section = 1;
            next;
        }
        
        ## skip other comments
        next if $line =~ /^\#/;
        
        if ( $in_fasta_section ) {
            if ( $line =~ /\>(\S+)/ ) {
                my $new_id = $1;
            
                ## we've hit a new sequence.  store the last one
                if ( $current_fasta_seq_id ) {
                    $$ann{$current_fasta_seq_id}{seq} = join('', @current_fasta_seq);
                    $$ann{$current_fasta_seq_id}{seq} =~ s|\s||g;
                }
                
                ## then reset
                $current_fasta_seq_id = $new_id;
                @current_fasta_seq = ();
            } else {
                push @current_fasta_seq, $line;
            }
        
        ## within here we must be in feature definitions
        } else {
            my @cols = split("\t", $line);
            next unless scalar @cols == 9;
            
            my $feat_id;
            
            if ( $cols[8] =~ /ID\=(.+?)\;/ || $cols[8] =~ /ID\=(.+?)/ ) {
                $feat_id = $1;
            } else {
                die "failed to find an ID on the following GFF line: $line\n";
            }
            
            my ($fmin, $fmax, $strand) = get_coordinates_from_gff(\@cols);
            
            if ( $cols[2] eq 'gene' ) {
                ## store the existing gene
                if ( defined $current_gene ) {
                    push @{$$ann{$current_gene_mol}{genes}}, $current_gene;
                }
                
                ## then reset
                $current_gene_mol = $cols[0];
                $current_gene = {
                    id      => $feat_id,
                    fmin    => $fmin,
                    fmax    => $fmax,
                    strand  => $strand,
                    exons  => [],
                };
            } elsif ( $cols[2] eq 'exon' ) {
                push @{$$current_gene{exons}}, { fmin => $fmin, fmax => $fmax };
                
            } elsif ( $cols[2] eq 'CDS' ) {
                push @{$$current_gene{CDS}}, { fmin => $fmin, fmax => $fmax };
            } 
        }
    }
    
    ## don't forget to check the last FASTA entry, if we're processing it
    if ( $current_fasta_seq_id ) {
        print STDERR "INFO: adding sequence for ID ($current_fasta_seq_id) length: ";
        $$ann{$current_fasta_seq_id}{seq} = join('', @current_fasta_seq);
        $$ann{$current_fasta_seq_id}{seq} =~ s|\s||g;
        print STDERR length($$ann{$current_fasta_seq_id}{seq}) . "\n";
        
    }
    
    ## store the last gene
    if ( $current_gene ) {
        push @{$$ann{$current_gene_mol}{genes}}, $current_gene;
    }
}

## trivial, but saves a bit of coding since it's needed so many times
sub get_coordinates_from_gff {
    my $cols = shift;
    
    my ($fmin, $fmax) = ($$cols[3] - 1, $$cols[4]);
    my $strand = $$cols[6] eq '+' ? 1 : -1;
    
    return ($fmin, $fmax, $strand);
}


sub parse_annotation_files {
    my ($files, $format) = @_;
    
    ## see ANNOTATION DATA STRUCTURE perldoc section for these
    my $annotation = {};
    
    for my $file ( @$files ) {
        _log("INFO: processing file $file");
        
        if ($format eq 'gff3') {
            parse_gff3_file($file, $annotation);
        } else {
            die "Unknown file format.  You shouldn't see this, since it should have been checked while parsing parameters.";
        }
    }
    
    return $annotation;
}


sub get_files_from_opts {
    my $opts = shift;
    
    my @files = ();
    
    if ( $$opts{input_file} ) {
        push @files, $$opts{input_file};
    }
    
    if ( $$opts{input_list} ) {
        open(my $ifh, $$opts{input_list}) || die "failed to open input list: $!";
        
        while (my $line = <$ifh>) {
            chomp $line;
            
            if ( $line =~ /\S/ ) {
                push @files, $line;
            }
        }
    }
    
    _log("INFO: found " . scalar(@files) . " file(s) to process");
    
    return \@files;
}


sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}


sub check_parameters {
    my $options = shift;
    
    ## input must have been specified
    unless ( $$options{input_file} || $$options{input_list} ) {
        die "You must specify input with either --input_file or --input_list";
    }
    
    ## make sure required arguments were passed
    #my @required = qw( format );
    #for my $option ( @required ) {
    #    unless  ( defined $$options{$option} ) {
    #        die "--$option is a required option";
    #    }
    #}
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    $options{format} = 'gff3' unless defined $options{format};
    
    my %accepted_input_types = (gff3 => 1);
    $$options{format} = lc $$options{format};
    
    if (! exists $accepted_input_types{$$options{format}} ) {
        die "$$options{format} is not a supported input format. Please specify one of the following: " .
            join( ' ', keys %accepted_input_types);
    }
    
    ## handle some defaults
    #$options{optional_argument2}   = 'foo'  unless ($options{optional_argument2});
}





















