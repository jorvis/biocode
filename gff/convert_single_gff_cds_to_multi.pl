#!/usr/bin/env perl

=head1 NAME

convert_single_gff_cds_to_multi.pl - Converts a GFF file between coding conventions related to the
CDS feature in GFF.

=head1 SYNOPSIS

USAGE: convert_single_gff_cds_to_multi.pl 
            --input=/path/to/some_file.gff 
            --output=/path/to/new_file.gff

=head1 OPTIONS

B<--input,-i>
    This input GFF expects the convention of one CDS per gene.

B<--output,-o>
    Output file to be created, where each CDS is the coding portion on an exon

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

The GFF format is great as long as the usage conventions described in the specification are 
followed, which is rarely the case.  This script corrects a specific issue I often find where
each gene only has a single CDS that spans multiple exons.  This script corrects this to 
follow the GFF3 canonical specification where each CDS represents the coding portion of each exon:

    http://www.sequenceontology.org/gff3.shtml

It is not magic and will not auto-correct any other errors your GFF has.  Check the INPUT
section below for expected input conventions.

=head1  INPUT

Example of the input type where one CDS spans multiple exons:

In terms of the order of features in the file, this script only assumes that a feature is defined
in a line before it is referenced as a Parent.

It also assumes the following relationships

    mRNA parent of exon
    mRNA parent of CDS
    gene parent of mRNA
    
    1 mRNA per gene
    1 (input) CDS per gene

=head1  OUTPUT

Comment lines are retained, but will all be at the top of the document.  This is conventional for 
much GFF header information but will completely foobar any file where comments embedded between 
feature lines is important.

Any other FASTA or non-columnar lines will be placed at the end of the output file.

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'input|i=s',
                          'output|o=s',
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

open(my $ifh, "<$options{input}") || die "can't read input file: $!";
open(my $ofh, ">$options{output}") || die "can't read output file: $!";

my %genes = ();
my @comment_lines = ();
my @other_feats = ();
my @other_lines = ();

## for ID lookups
my %mRNA2gene = ();


while (my $line = <$ifh>) {
    chomp $line;
    
    if ( $line =~ /^\#/ ) {
        push @comment_lines, $line;
        next;
    }
    
    my @cols = split("\t", $line);
    
    if ( scalar @cols < 9 ) {
        push @other_lines, $line;
        next;
    }
    
    my $feat_type = $cols[2];
    my ($id) = $cols[8] =~ /^ID\=(.+?)\;/;
    my $gene_id;
    
    if ( $feat_type eq 'mRNA' ) {
        ($gene_id) = $cols[8] =~ /Parent\=(.+?)\;/;
        $mRNA2gene{$id} = $gene_id;
        
    } elsif ( $feat_type eq 'CDS' || $feat_type eq 'exon' ) {
        my ($parent) = $cols[8] =~ /Parent\=(.+?)\;/;
        
        if ( ! exists $mRNA2gene{ $parent } ) {
            die ("failed to lookup gene ID for $parent");
        }
        
        $gene_id = $mRNA2gene{ $parent };
    }
    
    if ( $feat_type eq 'gene' ) {
        $genes{ $id }{line} = $line;
    
    } elsif ( $feat_type eq 'mRNA' ) {
        $genes{$gene_id}{mRNA} = \@cols;
    
    } elsif ( $feat_type eq 'CDS' ) {
        $genes{$gene_id}{CDS_id} = $id;
        $genes{$gene_id}{CDS_min} = $cols[3];
        $genes{$gene_id}{CDS_max} = $cols[4];
        
    } elsif ( $feat_type eq 'exon' ) {
        push @{ $genes{$gene_id}{exons} }, \@cols;
    
    } elsif ( $feat_type eq 'contig' ) {
        push @other_feats, $line;
        
    } else {
        die("panic!  don't know what to do with feat type: $feat_type");
    }
}

## first write the comments to the output file
for ( @comment_lines ) {
    print $ofh "$_\n";
}

## write features not handled in this script
for ( @other_feats ) {
    print $ofh "$_\n";
}

## write out the annotations
for my $g_id ( keys %genes ) {
    print $ofh "$genes{$g_id}{line}\n";
    
    ## print the mRNA
    if ( defined $genes{$g_id}{mRNA} ) {
        print $ofh join("\t", @{$genes{$g_id}{mRNA}}), "\n";
    } else {
        print "no mRNA found for $g_id\n";
    }
    
    my $cds_id = $genes{$g_id}{CDS_id};
    
    ## print the CDS and exons
    for my $colref ( @{ $genes{$g_id}{exons} } ) {
        
        my @cds_row = @$colref;
        $cds_row[2] = 'CDS';
        $cds_row[8] =~ s|ID=.+?\;|ID=$cds_id\;|;
        
        ## contained cds fragment
        if ( $$colref[3] >= $genes{$g_id}{CDS_min} && $$colref[4] <= $genes{$g_id}{CDS_max} ) {
            #print "CONTAINED: $$colref[3] >= $genes{$g_id}{CDS_min} && $$colref[4] <= $genes{$g_id}{CDS_max} ) {\n";
            ## just print it
            print $ofh join("\t", @cds_row), "\n";
       
        ## two cases here
        } elsif ( $$colref[3] < $genes{$g_id}{CDS_min} ) {
            $cds_row[3] = $genes{$g_id}{CDS_min};
            
            if ( $$colref[4] > $genes{$g_id}{CDS_max} ) {
                $cds_row[4] = $genes{$g_id}{CDS_max};
            }
            
            print $ofh join("\t", @cds_row), "\n";
            
        ## left-side overlap
        } elsif ( $$colref[4] > $genes{$g_id}{CDS_max} && $$colref[3] <= $genes{$g_id}{CDS_max} ) {
            #print "LEFT: $$colref[4] > $genes{$g_id}{CDS_max} && $$colref[3] <= $genes{$g_id}{CDS_max} ) {\n";
            $cds_row[4] = $genes{$g_id}{CDS_max};
            print $ofh join("\t", @cds_row), "\n";
        
        }

        ## print the actual exons
        print $ofh join("\t", @$colref), "\n";
    }
}

## finally write 'other' lines out (usually FASTA)
for ( @other_lines ) {
    print $ofh "$_\n";
}

exit(0);


sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( input output );
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
