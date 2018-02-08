#!/usr/bin/env perl

=head1 NAME

correct_gff_feature_order.pl - Reorders rows of a GFF file to conform with sane, expected norms.

=head1 SYNOPSIS

USAGE: correct_gff_feature_order.pl 
            --input=/path/to/some_file.gff 
            --output=/path/to/new_file.gff

=head1 OPTIONS

B<--input,-i>
    This input GFF expects the convention of one CDS per gene.

B<--output,-o>
    Output file to be created.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

While the GFF3 specification doesn't necessarily impose a restriction on the order of features 
within the document, a more defined ordering can help downstream software because of practical 
considerations.  GFF entries are often processed in groups where each 'group' contains the features
that define a gene and its children.  If these are scattered all over the document parsing
scripts have to first load the whole thing in memory.

=head1  INPUT and OUTPUT

This script makes very few assumptions on the overall contents of the input file.  The following
columns are used for the final ordering (and applied in this order:)

    1 - molecule
    3 - feature type
    4 - min location (of gene feature in each group)

Parental relationships assumed/supported:

    mRNA parent of exon
    mRNA parent of CDS
    gene parent of mRNA

Comment lines are retained, but will all be at the top of the document.  This is conventional for 
much GFF header information but will completely foobar any file where comments embedded between 
feature lines is important.

WARNING: this script does not currently handle multiple features on different lines with the same
ID value.  They'll simply be collapsed, which is doubtfully the desired behavior.

Remember that if you have a file which is nothing but match features with no parental
relationships you can easily do a sort on the command line like this:

    sort -t $'\t' -k 1,1 -k 4n,4 -k 5n,5 foo.gff3 > foo.sorted.gff3

Any other FASTA or non-columnar lines will be placed at the end of the output file.

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

## structure like:
my @comment_lines = ();
my @other_lines = ();
my @other_feats = ();

## for all gene graph features [ gene, mRNA, CDS, exon ]
## data structure like:
#    $h{$id} => {
#        type => 'gene',
#        parent => $id || undef,
#        cols => [columns],
#        
#        ## only gene type have this
#        #   hashref of ids of child elements of this gene
#        children => { cds => [], exons => [], mRNAs => [] },
#    }
my %features = ();

## keeps all the row information for segmented features.
## data structure like:
#   $h{$id} = [ [columns], [columns], ...  ]
my %segmented_features = ();

## used to keep references of genes on molecules for sorting purposes
my %molecules = ();


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
    my $id = gff3_get_column_9_value( $cols[8], 'ID' );
    #my ($parent) = $cols[8] =~ /Parent\=(.+?)\;/;
    my $parent = gff3_get_column_9_value( $cols[8], 'Parent' );
    
    if ( exists $features{$id} ) {
        push @{$segmented_features{$id}}, \@cols;
    } else {
        $features{$id} = {
            type => $feat_type,
            parent => $parent,
            cols => \@cols,
        };
    }
    
    if ( $feat_type eq 'gene' ) {
        ## make sure the molecule doesn't have a gene txhere already
        if ( exists $molecules{$cols[0]} && $molecules{$cols[0]}{$cols[3]} ) {
            die "found more than one gene at position $cols[3] on molecule $cols[0]\n";
        } else {
            $molecules{$cols[0]}{$cols[3]} = $id;
        }
    
    } elsif ( $feat_type eq 'mRNA' ||
              $feat_type eq 'CDS' || 
              $feat_type eq 'exon' || 
              $feat_type eq 'five_prime_UTR' ||
              $feat_type eq 'polypeptide' ||
              $feat_type eq 'three_prime_UTR' ||
              $feat_type eq 'tRNA' ||
              $feat_type eq 'rRNA' ||
              $feat_type eq 'non_canonical_five_prime_splice_site' ||
              $feat_type eq 'non_canonical_three_prime_splice_site'
        ) {
    
        ## nothing to see here, carry on
        
    
    } elsif ( $feat_type eq 'contig' ||
              $feat_type eq 'region' ) {
        
        push @other_feats, $line;
    
    } else {
        die("panic!  don't know what to do with feat type: $feat_type");
    }
}

## first write the comments to the output file
for ( @comment_lines ) {
    print $ofh "$_\n";
}

## pair genes with their child features
for my $feat_id ( keys %features ) {

    my $parent_id;
    if ( $features{$feat_id}{type} eq 'gene' ) {
        $parent_id = undef;
    } else {
        $parent_id = ultimate_parent_of($features{$feat_id}{parent});
    }
    
    my $type = $features{$feat_id}{type};

    ## if it has a parent, make sure we know about it
    if (defined $parent_id && ! exists $features{ $parent_id } ) {
        die "failed to find parent ID ($parent_id) for feature $feat_id";
    }
    
    if ( $features{$feat_id}{type} ne 'gene' ) {
        #print "adding: features{$parent_id}{children}{$type} }, $feat_id;\n";
        push @{ $features{$parent_id}{children}{$type} }, $feat_id;
    }
}

for my $molecule_id ( sort keys %molecules ) {

    print "DEBUG: Found " . scalar(keys %{$molecules{$molecule_id}}) . " genes to export on molecule $molecule_id\n";
    
    for my $start ( sort { $a<=>$b } keys %{$molecules{$molecule_id}} ) {
        my $feat_id = $molecules{$molecule_id}{$start};
        
        ## print the gene first
        print $ofh join("\t", @{$features{$feat_id}{cols}}), "\n";

        ## 5' UTRs
        for my $utr_id ( @{$features{$feat_id}{children}{five_prime_UTR}} ) {
            print $ofh join("\t", @{$features{$utr_id}{cols}}), "\n";
        }
        
        ## then the mRNAs
        for my $mRNA_id ( @{$features{$feat_id}{children}{mRNA}} ) {
            print $ofh join("\t", @{$features{$mRNA_id}{cols}}), "\n";
        }

        ## then the rRNAs
        for my $rRNA_id ( @{$features{$feat_id}{children}{rRNA}} ) {
            print $ofh join("\t", @{$features{$rRNA_id}{cols}}), "\n";
        }

        ## then the tRNAs
        for my $tRNA_id ( @{$features{$feat_id}{children}{tRNA}} ) {
            print $ofh join("\t", @{$features{$tRNA_id}{cols}}), "\n";
        }

        for my $child_type ( qw(exon CDS polypeptide ) ) {
            if ( exists $features{$feat_id}{children}{$child_type} ) {
                for my $child_id ( sort { $features{$a}{cols}[3] <=> $features{$b}{cols}[3] } @{$features{$feat_id}{children}{$child_type}} ) {
                    print $ofh join("\t", @{$features{$child_id}{cols}}), "\n";

                    ## if this is a segmented feature, print those too
                    if (exists $segmented_features{$child_id}) {
                        for my $child ( $segmented_features{$child_id} ) {
                            for my $elem ( sort {$$a[3] <=> $$b[3]} @$child ) {
                                print $ofh join("\t", @$elem), "\n";
                            }
                        }
                    }
                }
            }
        }

        ## 3' UTRs
        for my $utr_id ( @{$features{$feat_id}{children}{three_prime_UTR}} ) {
            print $ofh join("\t", @{$features{$utr_id}{cols}}), "\n";
        }
    }
}




## write features not handled in this script
for ( @other_feats ) {
    print $ofh "$_\n";
}



## finally write 'other' lines out (usually FASTA)
for ( @other_lines ) {
    print $ofh "$_\n";
}

exit(0);

sub ultimate_parent_of {
    my $child_id = shift;
    
    if ( $features{$child_id}{parent} ) {
        return ultimate_parent_of( $features{$child_id}{parent} );
    } else {
        return $child_id;
    }
}


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
