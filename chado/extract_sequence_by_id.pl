#!/usr/bin/env perl

=head1 NAME

extract_sequence_by_id.pl - extract the sequence of a feature in Chado by its ID.  Currently this
is limited to feature.uniquename, but will be expanded as needed.

=head1 SYNOPSIS

USAGE: extract_sequences_by_id.pl
            --database=eha3
            --user=someuser
            --password=somepass
          [ --uniquename=eha3.assembly.3425.1
            --database_type=mysql
            --server=yourserverhost
            --output_file=/path/to/somefile.fsa
            --log=/path/to/some.log ]


=head1 OPTIONS

B<--database,-d>
    Database name to connect to.

B<--user,-u>
    User account with select privileges on the specified database.

B<--password,-p>
    Password for user account specified.

B<--uniquename,-n>
    Optional.  One method of identifying sequences to extract.  Corresponds to feature.uniquename

B<--database_type>
    Optional.  Database type (vendor.)  Currently supports 'mysql' or 'postgresql' (default = mysql).

B<--server,-s>
    Optional.  Database server to connect to (default = localhost).

B<--output_file,-o>
    Optional.  Can specify the output file (default = STDOUT).

B<--log,-l> 
    Optional.  Full path to a log file to create.

B<--help,-h>
    This help message

=head1  DESCRIPTION

The goal of this script is to allow the extraction of the feature of any sequence in a Chado
database.  Currently rather limited, it allows selection by feature.uniquename value only
from a MySQL or Pg database.

=head1  INPUT

Required input includes the information necessary to connect to the database as well as
feature identifying info.  Current methods include:

    uniquename - the value in the feature.uniquename field.

This also assumes that feature.residues is also populated.  The methods for access will 
expand as needed.

=head1  OUTPUT

Features will be exported in FASTA format either on STDOUT or to a file (if --output_file is
passed)

Each FASTA sequence line contains 60 characters.

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
                          'user|u=s',
                          'password|p=s',
                          'uniquename|n=s',
                          'database_type|a=s',
                          'server|s=s',
                          'output_file|o=s',
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

_log("connecting to the database.");

my $dsn = '';

if ( $options{database_type} eq 'mysql' ) {
    $dsn = "dbi:mysql:database=$options{database};host=$options{server}";

} elsif ( $options{database_type} eq 'postgresql' ) {
    $dsn = "DBI:Pg:dbname=$options{database};host=$options{server}";
}

_log("attempting to create database connection");
my $dbh = DBI->connect($dsn, $options{user}, $options{password}, {PrintError=>1, RaiseError=>1} );

## manage the output file
my $ofh;

if ( $options{output_file} ) {
    
    open($ofh, ">$options{output_file}") || die "failed to create output file: $!";
    
} else {
    $ofh = *STDOUT;
}

## whatever method is used to pull sequences needs to result in the following columns (@results):
#   0. identifier   1. descriptive name (optional)  2. sequence
my $qry;
my @qry_opts = ();
my @results;

if ( $options{uniquename} ) {
    $qry = qq{
        SELECT uniquename, '', residues
          FROM feature
         WHERE uniquename = ?
    };
    
    push @qry_opts, $options{uniquename};
}

my $dsh = $dbh->prepare($qry);
   $dsh->execute(@qry_opts);

while ( my $row = $dsh->fetchrow_arrayref ) {
    @results = @$row;
    last;
}

$dsh->finish();

if ( scalar(@results) == 3 ) {
    
    print $ofh ">$results[0] $results[1]\n";
    
    ## format in lines of length 60
    while ($results[2] =~ /(.{1,60})/g) {
        print $ofh "$1\n";
    }
    
    print $ofh "\n";
    
} else {
    die("failed to find a sequence with the passed criteria");
}

$dbh->disconnect();

exit(0);



sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( database user password );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ## one of the following must be passed
    unless ( defined $$options{uniquename} ) {
        die "you must define IDs to extract using the --uniquename option."
    }
    
    ## handle some defaults
    $$options{server} = 'localhost' unless defined $$options{server};
    $$options{database_type} = 'mysql' unless defined $$options{database_type};
}













