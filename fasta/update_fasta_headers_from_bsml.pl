#!/usr/bin/env perl

=head1 NAME

update_fasta_headers_from_bsml.pl - Add annotation attributes to FASTA headers

=head1 SYNOPSIS

USAGE: fasta_size_filter.pl 
            --fasta_input_list=/path/to/some_file.list
            --bsml_input_list=/path/to/some_file.list 
            --output_file=/path/to/some_file.fsa 

=head1 OPTIONS

B<--fasta_input_list,-f>
    A text file that contains the paths to all the FASTA input files, one per line.

B<--bsml_input_list,-b>
    A text file that contains the paths to all the BSML input files, one per line.

B<--output_file,-o>
    File in which to write output.

B<--format,-r>
    Format in which you'd like the output FASTA headers.  See the OUTPUT section below
    for details.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

The point of this script is to give the user complete control over the headers in a
FASTA file by allowing the specification of the format on the command line.  This
is done by using a set of tokens that will be replaced by the actual values from
each annotated gene.  See the OUTPUT section for a list of these tokens and
examples.

=head1  INPUT

This script assumes that the FASTA headers are like:

    >gbb118a.polypeptide.1002.0 

Where the values here correspond to the id attributes in the polypeptide features of the BSML file.  Any content
after the ID portion of the source headers will be ignored.

For the BSML, an example entry with full annotation is:

    <Feature class="transcript" id="gbb118a.transcript.750.0">
      <Attribute name="assignby" content="ysebasti"></Attribute>
      <Attribute name="auto_annotate_toggle" content="0"></Attribute>
      <Attribute name="auto_comment" content="IDENT INFO:\ncom_name = phosphate ABC transporter, permease protein from RF|NP_212350.1|15594561|NC_001318\ngene_sym =  from RF|NP_212350.1|15594561|NC_001318\nec       =  from RF|NP_212350.1|15594561|NC_001318\nrole_id = 143 from guess_role\ngo_term =  from \nautoannotate = 0\nBEST BER INFO:\nacc: RF|NP_212350.1|15594561|NC_001318\ncom_name: phosphate ABC transporter, permease protein\nspecies of closest match: Borrelia burgdorferi B31\n\n"></Attribute>
      <Attribute name="date" content="Mar  9 2009  2:52PM"></Attribute>
      <Attribute name="gene_product_name" content="phosphate ABC transporter, permease protein"></Attribute>
      <Attribute name="gene_symbol" content="flgD"></Attribute>
      <Attribute name="start_site_editor" content="auto"></Attribute>
      <Interval-loc complement="0" startpos="0" endpos="441"></Interval-loc>
      <Cross-reference database="TIGR_Gbb118a" identifier="ORF00214" id="_10" identifier-type="feat_name"></Cross-reference>
      <Cross-reference database="TIGR_Gbb118a" identifier="BBU118A_0221" id="_11" identifier-type="locus"></Cross-reference>
      <Cross-reference database="TIGR_Gbb118a" identifier="BBU118A_0221" id="_12" identifier-type="display_locus"></Cross-reference>
      <Cross-reference database="Genbank" identifier="EEG98688" id="_13"></Cross-reference>
      <Cross-reference database="NCBI_gi" identifier="225369235" id="_14"></Cross-reference>
      <Link rel="sequence" href="#gbb118a.transcript.750.0_seq"></Link>
      <Attribute-list>
        <Attribute name="TIGR_role" content="143"></Attribute>
        <Attribute name="assignby" content="autoAnno"></Attribute>
        <Attribute name="date" content="Mar  4 2009  3:38PM"></Attribute>
      </Attribute-list>
    </Feature>

=head1  OUTPUT

The --format option allows the user to specify the exact format of the output FASTA
headers through the use of a specific set of tokens that are replaced by the
annotated values for each polypeptide.

If you start with a header like this:

    >gbb118a.polypeptide.999.0 Some annotation here
    
You could use the --format value "-=LOCUS=- -=XREF_GENBANK=- -=PRODUCT=- (-=GENESYM=-)" to produce:

    >BBU118A_0221 EEG98688 phosphate ABC transporter, permease protein (flgD)

Notice that you can enter any characters in the string you want (like the parentheses above)
and they will be maintained.   Only the tokens in the format -=LIKETHIS=- will be
replaced.

This script can only pull what is in the associated BSML files.  So if you specify
any tokens in the output string that aren't actually present for any given gene that
portion will just be empty.

The full list of currently supported tokens:

    ORIG_ID - the original ID in the header without reformatting
    ORIG_DESC - the original content of the header after ID and first space
    PRODUCT - the gene product name
    GENESYM - the gene symbol
    FEAT_NAME - identifier used in the TIGR legacydb for this feature, if available.
    LOCUS - locus specifier
    XREF_? - The ? is replaced by any database (Genbank, NCBI_gi, etc.) like above


=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use XML::Twig;

my %options = ();
my $results = GetOptions (\%options, 
                          'fasta_input_list|f=s',
                          'bsml_input_list|b=s',
                          'output_file|o=s',
                          'format|r=s',
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

my @input_fasta_seqs = list_file_members( $options{fasta_input_list} );
my @input_bsml_seqs = list_file_members( $options{bsml_input_list} );

## the output file handle
my $ofh;
if ( defined $options{output_file} ) {
    $ofh = get_conditional_write_fh( $options{output_file} );
}

## keyed on transcript ID, contains subkeys:
#       product
#       symbol
#
#   loading will error if there are duplicates
my %annotations;

## id lookup
my %polypeptide2transcript;

## this doesn't include those starting with XREF_, which are variable
my %known_keys = (
    'ORIG_ID' => 1,
    'ORIG_DESC' => 1,
    'PRODUCT' => 1,
    'GENESYM' => 1,
    'FEAT_NAME' => 1,
    'LOCUS' => 1,
);



for my $bsml_file ( @input_bsml_seqs ) {
    my $ifh = get_conditional_read_fh( $bsml_file );
    
    my $twig = XML::Twig->new(
        twig_handlers => {
            'Feature[@class="transcript"]' => \&process_transcript,
            'Feature-group' => \&process_feature_group,
        },
    );
    
    print STDERR "processing file $bsml_file\n";
    $twig->parse( $ifh );
}


for my $fasta_file ( @input_fasta_seqs ) {
    my $ifh = get_conditional_read_fh( $fasta_file );
    
    while ( my $line = <$ifh> ) {
        chomp $line;
    
        if ( $line =~ /\>(\S+)\s*(.*)/ ) {
            
            my $orig_id = $1;
            my $orig_desc = $2;
            my $orig_line = $line;
            $line = ">$options{format}";
            
            my $t_id; # transcript ID
            
            ## bad things if we don't know what this is
            if ( exists $polypeptide2transcript{$orig_id} ) {
                $t_id = $polypeptide2transcript{$orig_id};
            } else {
                die "ID found in FASTA ($orig_id) not found in BSML polypeptide -> transcript lookup.";
            }
            
            if ( ! exists $annotations{ $t_id } ) {
                die "Found mapping, but no annotation entry for transcript ($t_id)";
            }
            
            ## process the keys
            ###################
            if ( $options{format} =~ /-=ORIG_ID=-/ ) {
                $line =~ s/-=ORIG_ID=-/$orig_id/g;
            }
            
            if ( $options{format} =~ /-=ORIG_DESC=-/ ) {
                $line =~ s/-=ORIG_DESC=-/$orig_desc/g;
            }
            
            if ( $options{format} =~ /-=PRODUCT=-/ ) {
                $line =~ s/-=PRODUCT=-/$annotations{$t_id}{product}/g;
            }
            
            if ( $options{format} =~ /-=GENESYM=-/ ) {
                $line =~ s/-=GENESYM=-/$annotations{$t_id}{symbol}/g;
            }
            
            if ( $options{format} =~ /-=FEAT_NAME=-/ ) {
                $line =~ s/-=FEAT_NAME=-/$annotations{$t_id}{feat_name}/g;
            }
            
            if ( $options{format} =~ /-=LOCUS=-/ ) {
                $line =~ s/-=LOCUS=-/$annotations{$t_id}{locus}/g;
            }
            
            while ( $options{format} =~ /-=(XREF_([^=]+))=-/g ) {
                my $orig_db = $1;
                my $db = lc($2);
                
                my $xref = $annotations{$t_id}{"xref_$db"};
                $line =~ s/-=$orig_db=-/$xref/g;
            }
        }

        ## strip the line of consecutive whitespace (could happen with token gaps)
        $line =~ s/\s+/ /g;
        
        ## strip the end of the line
        $line =~ s/\s+$//;
        
        print $ofh "$line\n";
    }
}


exit(0);

sub process_transcript {
    my ($twig, $transcript) = @_;
    
    my $transcript_id = $transcript->{att}->{id};
    
    ## make sure we haven't seen this already
    if ( exists $annotations{$transcript_id} ) {
        die "ERROR: found duplicate instance of transcript ($transcript_id) - this shouldn't happen.";
    }
    
    $annotations{$transcript_id} = {};
    
    for my $attribute ( $transcript->children('Attribute') ) {
        my $name = $attribute->{att}->{name};
        
        if ( $name eq 'gene_product_name' ) {
            $annotations{$transcript_id}{product} = $attribute->{att}->{content};
            
        } elsif ( $name eq 'gene_symbol' ) {
            $annotations{$transcript_id}{genesym} = $attribute->{att}->{content};
        } 
    }
    
    for my $xref ( $transcript->children('Cross-reference') ) {
        if ( defined $xref->{att}->{'identifier-type'} ) {
            $annotations{$transcript_id}{ $xref->{att}->{'identifier-type'} } = $xref->{att}->{identifier};
        } else {
            #print "adding an xref_annotation: nnotations{$transcript_id}{ 'xref_' . lc($xref->{att}->{'database'}) } = $xref->{att}->{identifier};\n";
            $annotations{$transcript_id}{ ('xref_' . lc($xref->{att}->{'database'})) } = $xref->{att}->{identifier};
        }
    }
}


sub process_feature_group {
    my ($twig, $fg) = @_;
    
    my ($transcript_id, $polypeptide_id);
    
    ## only do this if it's gene grouping
    next unless $fg->{att}->{'group-set'} =~ /\.gene\./;
    
    for my $fgm ( $fg->children('Feature-group-member') ) {
        if ( $fgm->{att}->{'feature-type'} eq 'transcript' ) {
            $transcript_id = $fgm->{att}->{'featref'};
            
        } elsif ( $fgm->{att}->{'feature-type'} eq 'polypeptide' ) {
            $polypeptide_id = $fgm->{att}->{'featref'};
        }
    }
    
    if ( defined $transcript_id && defined $polypeptide_id ) {
        $polypeptide2transcript{$polypeptide_id} = $transcript_id;
    } else {
        die "failed to make transcript => polypeptide ID mapping. got '$transcript_id' -> '$polypeptide_id'";
    }
}

sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}

sub get_conditional_write_fh {
    my $path = shift;
    my $fh;
    
    if ( $path =~ /\.gz$/ ) {
        open( $fh, ">:gzip", "$path" ) || die "failed to create output file: $!";
    } else {
        open( $fh, ">$path" ) || die "failed to create output file: $!";
    }
    
    return $fh;
}

sub get_conditional_read_fh {
    my $path = shift;
    
    ## auto-handle gzipped files
    if (! -e $path && -e "$path.gz") {
        $path .= '.gz';
    }

    my $fh;
    if ( $path =~ /\.gz$/ ) {
        open($fh, "<:gzip", $path) || die "can't read file $path: $!";
    } else {
        open($fh, "<$path") || die "can't read file $path: $!";
    }
    
    return $fh;
}

sub list_file_members {
    my $list_file = shift;
    my @seqs;
    
    if ( -e $list_file ) {

        open(my $list_fh, "<$list_file") || die "can't open input list: $!";

        while ( <$list_fh> ) {
            next if ( /^\s*$/ );

            chomp;

            if ( -e $_ || -e "$_.gz") {
               push @seqs, $_;
            } else {
                die "the following input file doesn't exist: $_";
            }
        }

    } else {
        die "the following input list doesn't exist: $list_file";
    }
    
    return @seqs;
}

sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( fasta_input_list bsml_input_list output_file format );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
}









