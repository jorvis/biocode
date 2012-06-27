#!/usr/bin/env perl

=head1 NAME

report_gc_content_by_feature_type.pl - Reads a GFF file and reports sum length and GC content
aggregated by feature type

=head1 SYNOPSIS

USAGE: report_gc_content_by_feature_type.pl 
            --gff_file=/path/to/some_file.gff 
            --fasta_file=/path/to/some_file.fa
            --molecule_id=aa.assembly.81787
           

=head1 OPTIONS

B<--gff_file,-g>
    An input GFF file.  Should conform to the GFF3 spec, untested otherwise.

B<--fasta_file,-f>
    A single or multi-FASTA file that must contain an entry for --molecule_id

B<--molecule_id,-m>
    The ID of the molecule you're interested in.  This is the first column in the GFF
    file and must correspond to an identifier in the FASTA file passed.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Useful for analyzing sequence bias by different feature types, this script will
export a tab delimited report of each feature type found in a GFF file and some
statistics on each.  (see OUTPUT section below).  A few other types are also
calculated that aren't defined in the file, such as intron and telomere.

=head1  INPUT

A multi-FASTA file and a GFF file.

=head1  OUTPUT

Example output:

    ## Data by feature
    ## feature_type	summed_length	GC_count	GC_percent
    chromosome	2540030	867789	34.16
    gene	2103644	752352	35.76
    mRNA	2098384	749699	35.73
    polypeptide	2051944	738571	35.99
    intron	742678	192378	25.90
    telomere	6226	1882	30.23
    rRNA	4146	1989	47.97
    tRNA	1114	664	59.61

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'gff_file|g=s',
                          'fasta_file|f=s',
                          'molecule_id|m=s',
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

my $molecule_seq = get_molecule_sequence_from_fasta( $options{fasta_file}, $options{molecule_id} );

## make sure we got something for the molecule
if ( ! $molecule_seq ) {
    die "Failed to return a sequence for molecule id $options{molecule_id} in file $options{fasta_file}\n";
} else {
    print STDERR "INFO: loaded molecule $options{molecule_id} (" . length($molecule_seq) . " bp)\n";
}


open(my $gff_fh, $options{gff_file}) || die "failed to open GFF3 file: $!";

## each type here will have 3 values, feat_count, length, gc_count
my %features = ();
my %last_feat_fmax_by_type = ();

while (my $line = <$gff_fh>) {
    chomp $line;
    my @cols = split("\t", $line);
    
    ## gff data rows have 9 columns
    next if $line =~ /^\#/;
    next unless scalar(@cols) == 9;
    
    my $mol_id  = $cols[0];
    my $type    = $cols[2];
    my $fmin    = $cols[3] - 1;
    my $fmax    = $cols[4];
    
    $cols[8] =~ /ID\=(.+?)\;/;
    my $feat_id = $1 || $cols[8];
    
    ## quick sanity check
    if ($fmin > $fmax) {
        die "Error: on the following line, the stop coordinate is less than the start coordinate:\n$line";
    }

    next unless $mol_id eq $options{molecule_id};
    
    ## when we hit an exon, we need to calculate the previous intron and, if this
    #   is the first exon, the telomeric region
    if ( $type eq 'exon' ) {
        if ( defined $last_feat_fmax_by_type{exon} ) {
            ## it's possible to have overlapping exons on opposite strands
            ## if so, there's no intron to calculate
            if ( $fmin < $last_feat_fmax_by_type{exon} ) {
                print STDERR "Warning: overlapping exon at coordinates (${fmin}-$fmax)\n";
            } else {
                my $intron_seq = substr( $molecule_seq, $last_feat_fmax_by_type{exon}, $fmin - $last_feat_fmax_by_type{exon} );
                $features{intron}{length} += length $intron_seq;
                $features{intron}{gc_count} += $intron_seq =~ tr/gcGC//;
                $features{intron}{feat_count} += 1;
            }
            
        } else {
            my $telomere_seq = substr( $molecule_seq, 0, $fmin );
            $features{telomere}{length} += length $telomere_seq;
            $features{telomere}{gc_count} += $telomere_seq =~ tr/gcGC//;
            $features{telomere}{feat_count} += 1;
        }
    }
    
    my $feat_seq = '';
    
    if ( ! exists $last_feat_fmax_by_type{$type} || $fmax > $last_feat_fmax_by_type{$type} ) {
        ## don't double count if these is a feature of the same type overlapping this one
        if ( exists $last_feat_fmax_by_type{$type} && $fmin < $last_feat_fmax_by_type{$type} ) {
            $feat_seq = substr( $molecule_seq, $last_feat_fmax_by_type{$type}, $fmax - $last_feat_fmax_by_type{$type} );
        } else {
            $feat_seq = substr( $molecule_seq, $fmin, $fmax - $fmin );
        }

        $features{$type}{length} += length($feat_seq);
        $features{$type}{gc_count} += $feat_seq =~ tr/gcGC//;
        $features{$type}{feat_count} += 1;
        $last_feat_fmax_by_type{$type} = $fmax;
    } else {
        print STDERR "Warning: feature ($feat_id) at (${fmin}-$fmax) seems to be contained within previous feature of the same type.\n";
    }
    
}

## hack for kelly, who wants mRNAs to only be the exon portions and doesn't have files which
#   have CDS defined to deal with.
$features{mRNA}{gc_count} = $features{exon}{gc_count};
$features{mRNA}{length} = $features{exon}{length};

## calculate the telomeric region after the last exon
my $last_telomere_seq = substr($molecule_seq, $last_feat_fmax_by_type{exon}, length($molecule_seq) - $last_feat_fmax_by_type{exon} );
$features{telomere}{length} += length $last_telomere_seq;
$features{telomere}{gc_count} += $last_telomere_seq =~ tr/gcGC//;
$features{telomere}{feat_count} += 1;

print "## Data by feature\n";
print "## feature_type\tfeature_count\tsummed_length\tGC_count\tGC_percent\n";

## print them out, sorted by most prevalent
for my $type (sort {$features{$b}{length} <=> $features{$a}{length}} keys %features) {
    print "$type\t$features{$type}{feat_count}\t$features{$type}{length}\t$features{$type}{gc_count}\t" .
          sprintf("%.2f", ( ($features{$type}{gc_count}/$features{$type}{length}) * 100)) . "\n";
}


exit(0);


sub get_molecule_sequence_from_fasta {
    my ($file, $mol_id) = @_;
    
    open(my $ifh, $file) || die "failed to read FASTA file: $!";
    
    my @seq_parts = ();
    my $in_seq_region = 0;
    
    while (my $line = <$ifh> ) {
        chomp $line;
        
        if ( $line =~ /^>(\S+)/ ) {
            my $this_mol_id = $1;
            
            if ( $this_mol_id eq $mol_id ) {
                $in_seq_region = 1;
            } else {
                $in_seq_region = 0;
            }
        } else {
            push @seq_parts, $line if $in_seq_region;
        }
    }
    
    return join('', @seq_parts);
}


sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( gff_file fasta_file molecule_id );
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
