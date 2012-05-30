#!/usr/bin/env perl

=head1 NAME

extract_sequences_by_type.pl - extracts sequences of features in a chado database by their type.

=head1 SYNOPSIS

USAGE: extract_sequences_by_type.pl
            --database=eha3
            --user=someuser
            --password=somepass
            --feature_type=polypeptide
            --output_directory=/path/to/somedir
          [ --database_type=mysql
            --server=SYBTIGR 
            --organism_id=4
            --uniquename_filter="eha3%"
            --use_feature_names
            --output_list=/path/to/some.list
            --log=/path/to/some.log ]


=head1 OPTIONS

B<--database,-d>
    Database name to connect to.

B<--user,-u>
    User account with select privileges on the specified database.

B<--password,-p>
    Password for user account specified.

B<--feature_type,-f>
    Type of feature for which you want to extract sequences, such as assembly or
    polypeptide.  Corresponds to cvterm.name.

B<--output_directory,-o>
    Directory where the output files should be written.  This direcetory must 
    already exist.

B<--database_type>
    Optional.  Database type (vendor.)  Currently supports 'mysql' or 'postgresql' (default = mysql).

B<--server,-s>
    Optional.  Database server to connect to (default = localhost).

B<--organism_id,-r>
    Optional.  If you want to limit your extraction to a single organism ID (only
    applies to a database with multiple organisms loaded, such as the comparative
    databases.)

B<--length_cutoff,-e>
    Optional.  If you want to filter the sequences extracted by length cutoff you
    can specify the cutoff length, in base pairs, here (default = 1).

B<--output_list,-t>
    Optional.  Writes a list file containing the full paths to each of the
    sequence files written.

B<--uniquename_filter,-n>
    Optional.  Allows you to filter the results returned by the value in 
    feature.uniquename.  Pass a SQL-like string, like "eha3%"

B<--use_feature_names,-u>
    Optional, will use the feature.name as the primary identifier in the 
    output fasta file, instead of the uniquename.

B<--log,-l> 
    Optional.  Full path to a log file to create.

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script can be used to extract the sequences of any feature in a chado schema which
has the feature.residues field populated, writing an individual FASTA-formatted file
for each feature extracted.  If you attempt to extract a sequence type which has no
residues, the script will report an error and quit.

This script will not currently extract features whose sequences must be derived from
other features.  The script is also currently limited to features where feature.is_analysis
is set to 0.

=head1  INPUT

Required input includes the information necessary to connect to the database.  The 
feature_type option corresponds to the cvterm.name table, with values such as 'assembly'
or 'polypeptide' (corresponding to SO terms.)

=head1  OUTPUT

This script will note produce any terminal output but will write a log if the --log 
option is used.  A single file for each feature will be written to the output_directory
specified in FASTA format.  As an example, if you run the script to extract polypeptide
features from the EHA3 database the output file created may be called this:

    eha3.polypeptide.10037.1.fsa
    
And the contents of the file would look like this:

    >eha3.polypeptide.10037.1
    MIELLLLFFLFYCYSEKIDLKKNRLDVVKDKLNNVIDSTTAYLLLNNDKMIGSYGRKKLI
    NGKVIHLPTEEEIQIWKKSNDYNGILVVPDHLLDSTHLNLFEQLKVQGVVFYDSNSTHYV
    SQDYRSSFNPNASEVAWKKYTYPIMFTPTSIVYFTQNNGYLQIEDEMLADGDSEICLRRG
    MCQPIGGQSLVGRLVENHLSESILLSAQLDGVSIFRDATPAHEQIVSGEATLLSVLETLR
    PLMKGAKKGIDFGFFEGETFGNMGSKTFSKSNDSYSHIIHFGPLGLTEGSLFVYSKDQVP
    QDIHKRIVNKTNEYPHCALSSFKNGIKMYLSDHEELYKGNIGTHQDYLPYITELCDNIES
    ITQIILQIIFNRSDIETLQKELEIQYDCEKYHRLYHCLAYDLSCDYFNETIGIIQGNSNS
    YSSIYRETKITPRSKAIHDYLNRIVTMGVCKEECDTSIKNGTVTYLGSYGSDIDISKNKV
    VGSSGDAWTESNWKQVKLTFYQPRTFGDVTLITLGVGLTLMTIIAFFTVIFKLNLL*

Each line contains 60 characters.

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
$|++;

my %options = ();
my $results = GetOptions (\%options, 
                          'database|d=s',
                          'database_type|a=s',
                          'server|s=s',
                          'user|u=s',
                          'password|p=s',
                          'feature_type|f=s',
                          'organism_id|r=i',
                          'length_cutoff|e=i',
                          'output_directory|o=s',
                          'output_list|t=s',
                          'uniquename_filter|n=s',
                          'use_feature_names|u',
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

## create the list file if requested
my $listfh;
if (defined $options{output_list}) {
    open($listfh, ">$options{output_list}") || die "can't create output list file: $!";
}

_log("connecting to the database.");

my $dsn = '';

if ( $options{database_type} eq 'mysql' ) {
    $dsn = "dbi:mysql:database=$options{database};host=$options{server}";

} elsif ( $options{database_type} eq 'postgresql' ) {
    $dsn = "DBI:Pg:dbname=$options{database};host=$options{server}";
}

_log("attempting to create database connection");
my $dbh = DBI->connect($dsn, $options{user}, $options{password}, {PrintError=>1, RaiseError=>1} );

## get the maximum length for this feature type
my $max_feature_length = 0;
my $qry = "SELECT max(f.seqlen), count(f.feature_id) " .
          "FROM feature f, cvterm c " .
          "WHERE f.type_id=c.cvterm_id " .
          "  AND c.name = ? " .
          "  AND f.is_analysis = ? " .
          "  AND f.is_obsolete = ? ";
my $feature_length_selector = $dbh->prepare( $qry );
   $feature_length_selector->execute( $options{feature_type}, 0, 0 );
   
$max_feature_length = ( $feature_length_selector->fetchrow_array )[0] || 0;
_log("setting maximum feature length to $max_feature_length based on the selected features.");

$feature_length_selector->finish();

## quick length check
if (! $max_feature_length ) {
    _log("quitting.  no features of type '$options{feature_type}' were found with residue length > 0");
    die("no features of type '$options{feature_type}' were found with residue length > 0");
}

## set the textsize to support this max length
## this is Sybase-specific, which I don't care about supporting anymore
#$dbh->do("set textsize $max_feature_length");

$qry = "SELECT f.uniquename, f.residues, f.name " .
       "FROM feature f, cvterm c " .
       "WHERE f.type_id=c.cvterm_id " .
       "  AND c.name = ? " .
       "  AND f.is_analysis = ? " .
       "  AND f.is_obsolete = ? ";
       
if ( $options{organism_id} ) {
    $qry .= "  AND f.organism_id = $options{organism_id} ";
}

if ( defined $options{uniquename_filter} ) {
    $qry .= "  AND f.uniquename LIKE '$options{uniquename_filter}' ";
}

my $feature_selector = $dbh->prepare($qry);
   $feature_selector->execute( $options{feature_type}, 0, 0 );

while ( my $row = $feature_selector->fetchrow_hashref ) {
    my $featid = ( $options{'use_feature_names'} ) ? $row->{'name'} : $row->{'uniquename'};

    ## only do polypeptides that actually have a sequence
    if ( length $row->{residues} > 0 ) {

        if ( length $row->{residues} >= $options{length_cutoff} ) {
            my $of_path = "$options{output_directory}/$featid.fsa";
            _log("writing $featid to $of_path");
            open(my $ofh, ">$of_path") || die "can't create $of_path: $!";
            print $ofh ">$featid\n";

            ## format in lines of length 60
            while ($row->{residues} =~ /(.{1,60})/g) {
                print $ofh "$1\n";
            }

            ## add it to the list file if necessary
            print $listfh "$of_path\n" if ( $listfh );
        } else {
            print "STDERR skipping $featid\n";
            _log("skipping $featid because it is shorter than the length cutoff");        
        }
        
    } else {
        print STDERR "skipping $featid because it has no sequence\n";
        _log("skipping $featid because it has no sequence");
    }
}

$feature_selector->finish();

$dbh->disconnect();

exit(0);



sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( database user password feature_type output_directory );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ## database should be lower-case
    ## KG: why should the database be lower-case?
    #$$options{database} = lc $$options{database};
    
    ## handle some defaults
    $$options{server} = 'localhost' unless defined $$options{server};
    $$options{organism_id} = 0 unless defined $$options{organism_id};
    $$options{length_cutoff} = 1 unless defined $$options{length_cutoff};
    $$options{database_type} = 'mysql' unless defined $$options{database_type};
}
