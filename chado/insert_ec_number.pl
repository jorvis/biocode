#!/usr/bin/env perl

=head1 NAME

insert_ec_number.pl - Utility script to insert a single EC number into a Chado database

=head1 SYNOPSIS

USAGE: insert_ec_number.pl 
            --ec_number=5.3.3.13
            --name=Isomerases||Intramolecular oxidoreductases||Transposing C==C bonds||Polyenoic fatty acid isomerase
            --database=abs2
            --user=someuser
            --password=somepass
          [ --definition=(5Z,8Z,11Z,14Z,17Z)-eicosapentaenoate = (5Z,7E,8E,14Z,17Z)-(5Z,8Z,11Z,14Z,17Z)-eicosapentaenoate = (5Z,7E,8E,14Z,17Z)-eicosapentaenoate
            --server=localhost
            --database_type=postgresql
          ]

=head1 OPTIONS

B<--ec_number>
    Can be fully-qualified like 4.1.1.85 or partial like 4.1.-.-

B<--name>
    This is the display name for the EC number and is formatted as the path the names for each
    step in the path separated by || symbols.  See the example above.

B<--database>
    Database name to connect to.

B<--user>
    User account with select, insert, and update privileges on the specified database.

B<--password>
    Password for user account specified.

B<--definition>
    Optional.  Text string that represents the reaction catalyzed by this enzyme.

B<--server>
    Optional.  Server to connect to (default = localhost).

B<--database_type>
    Optional.  Database type (vendor.)  Currently supports 'mysql' or 'postgresql' (default = mysql).

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Some enterprising young annotator has manually curated an EC number that hasn't made it into
our update cycle yet, which will cause a lookup failure.  This script can be used to insert
it for maximum annotator happiness.

=head1  INPUT

Command-line options only.  Normally our EC numbers are loaded from an ontology, which provides
identifiers like 'EC:00004613'.  These aren't going to be available here, so the loader looks at
the highest ID in the database and increments from that.  This, of course, assumes that updates
either won't be done or will be done based on the actual EC number instead, and that these
internally generated identifiers won't be exported.

=head1  OUTPUT

For each EC number, does the following INSERTS into the database

    - one row in dbxref for the EC number (5.3.3.13)
    - one row in dbxref for the EC ID (EC:00004613)
    - one row in cvterm for the new term's name and definition
    - one row in cvterm_dbxref to link the ID with its annotation
    - one row in cvterm_relationship linking 5.3.3.13 -> 5.3.3.-

Example representation of this EC entry:

    [Term]
    id: EC:00004613
    name: Isomerases||Intramolecular oxidoreductases||Transposing C==C bonds||Polyenoic fatty acid isomerase
    def: "(5Z,8Z,11Z,14Z,17Z)-eicosapentaenoate = (5Z,7E,8E,14Z,17Z)-(5Z,8Z,11Z,14Z,17Z)-eicosapentaenoate = (5Z,7E,8E,14Z,17Z)-eicosapentaenoate"
    comment: "-!- The enzyme from the red alga Ptilota filicina catalyzes the-!- The enzyme from the red alga Ptilota filicina catalyzes theisomerization of skip dienes (methylene-interrupted double bonds) in-!- The enzyme from the red alga Ptilota filicina catalyzes theisomerization of skip dienes (methylene-interrupted double bonds) ina broad range of fatty acids and fatty-acid analogues, such as-!- The enzyme from the red alga Ptilota filicina catalyzes theisomerization of skip dienes (methylene-interrupted double bonds) ina broad range of fatty acids and fatty-acid analogues, such asarachidonate and gamma-linolenate, to yield a conjugated triene."
    synonym: "PFI.Eicosapentaenoate cis-delta(5,8,11,14,17)-eicosapentaenoate cis-delta(5)-trans-delta(7,9)-cis-delta(14,17) isomerase." []
    xref_analog: EC:5.3.3.13
    is_a: EC:00000323

    +-----------+-------+-------------+---------+-------------+
    | dbxref_id | db_id | accession   | version | description |
    +-----------+-------+-------------+---------+-------------+
    |      9980 |     8 | EC:00004613 | 1.0     | NULL        | 
    |      9981 |     8 | 5.3.3.13    | current | NULL        | 
    +-----------+-------+-------------+---------+-------------+

    +------------------+-----------+-----------+-------------------+
    | cvterm_dbxref_id | cvterm_id | dbxref_id | is_for_definition |
    +------------------+-----------+-----------+-------------------+
    |             4859 |      5122 |      9981 |                 0 | 
    +------------------+-----------+-----------+-------------------+

+-----------+-------+----------------------------------------------------------------------------------------------------+-----------+
| cvterm_id | cv_id | name                                                                                               | dbxref_id |
+-----------+-------+----------------------------------------------------------------------------------------------------+-----------+
|      5122 |     7 | Isomerases||Intramolecular oxidoreductases||Transposing C==C bonds||Polyenoic fatty acid isomerase |      9980 | 
+-----------+-------+----------------------------------------------------------------------------------------------------+-----------+


=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut

use strict;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'ec_number=s',
                          'name=s',
                          'definition=s',
                          'database=s',
                          'user=s',
                          'password=s',
                          'server=s',
                          'database_type=s',
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

my $dsn = '';

if ( $options{database_type} eq 'mysql' ) {
    $dsn = "dbi:mysql:database=$options{database};host=$options{server}";

} elsif ( $options{database_type} eq 'postgresql' ) {
    $dsn = "DBI:Pg:dbname=$options{database};host=$options{server}";
}

_log("attempting to create database connection");
my $dbh = DBI->connect($dsn, $options{user}, $options{password}, {PrintError=>1, RaiseError=>1} );

## populated in the dbxref.version field for EC:N ontology-ish IDs
my $ec_version = '1.0';

## $h{$name} = $id
my %db;
$db{EC} = get_db_id('EC');

my %cv;
$cv{EC} = get_cv_id('EC');

my %next_ids;

foreach ( qw( dbxref cvterm cvterm_dbxref cvterm_relationship) ) {
    $next_ids{$_} = get_next_canonical_id( $_ );
}

## number portion of new IDs like EC:00004613
my $next_ec_id_num = get_next_ec_id_num();

## reused
my $qry;

$qry = qq{
    INSERT INTO dbxref ( dbxref_id, db_id, accession, version )
    VALUES ( ?, $db{EC}, ?, ? )
};
my $dbxref_inserter = $dbh->prepare($qry);

$qry = qq{
    INSERT INTO cvterm ( cvterm_id, cv_id, name, definition, dbxref_id, is_obsolete, is_relationshiptype )
    VALUES (?, $cv{EC}, ?, ?, ?, 0, 0)
};
my $cvterm_inserter = $dbh->prepare($qry);

$qry = qq{
    INSERT INTO cvterm_dbxref (cvterm_dbxref_id, cvterm_id, dbxref_id, is_for_definition) 
    VALUES( ?, ?, ?, 0 )
};
my $cvterm_dbxref_inserter = $dbh->prepare($qry);

$qry = qq{
    INSERT INTO cvterm_relationship ( cvterm_relationship_id, type_id, subject_id, object_id )
    VALUES ( ?, ?, ?, ?);
};
my $cvterm_relationship_inserter = $dbh->prepare($qry);

my $base_cvterm_id = get_base_cvterm_id( $options{ec_number} );

## insert the EC id
my $ec_row_id = $next_ids{dbxref}++; ## dbxref.dbxref_id
my $ec_id = 'EC:' . sprintf("%08d", $next_ec_id_num++);
_log("INFO: inserting row into dbxref with dbxref_id: $ec_row_id");
$dbxref_inserter->execute( $ec_row_id, $ec_id, $ec_version );

my $cvterm_row_id = $next_ids{cvterm}++;
_log("INFO: inserting row into cvterm with cvterm_id: $cvterm_row_id");
$cvterm_inserter->execute( $cvterm_row_id, $options{name}, $options{definition}, $ec_row_id );

my $ecnum_row_id = $next_ids{dbxref}++;
_log("INFO: inserting row into dbxref with dbxref_id: $ecnum_row_id");
$dbxref_inserter->execute( $ecnum_row_id, $options{ec_number}, 'current' );

my $cvdbx_id = $next_ids{cvterm_dbxref}++;
_log("INFO: inserting row into cvterm_dbxref with cvterm_dbxref_id: $cvdbx_id");
$cvterm_dbxref_inserter->execute( $cvdbx_id, $cvterm_row_id, $ecnum_row_id );

my $cvtr_id = $next_ids{cvterm_relationship}++;
_log("INFO: inserting row into cvterm_relationship with cvterm_relationship_id: $cvtr_id");
$cvterm_relationship_inserter->execute( $cvtr_id, get_cvterm_id('is_a'), $cvterm_row_id, $base_cvterm_id );

$cvterm_relationship_inserter->finish();
$cvterm_inserter->finish();
$dbxref_inserter->finish();
$dbh->disconnect();

exit(0);

## for an EC number like 5.3.3.13 this returns the cvterm corresponding to the entry for
#  5.3.3.-.   If 5.3.3.- was passed, would return that for 5.3.-.- 
sub get_base_cvterm_id {
    my $child_ec = shift;
    my $base_ec;
    
    my @parts = split('\.', $child_ec);
    my @new = ();
    
    ## this must be 4 parts
    if ( scalar @parts != 4 ) {
        die "encountered an EC number that doesn't have 4 parts: $child_ec, @parts";
    }
    
    foreach my $i ( reverse(0..3) ) {
        $new[$i] = '-';
    
        if ( $parts[$i] != '-' ) {
            ## slurp the rest in
            map $new[$_]=$parts[$_], ( reverse 0 .. ($i - 1) );
            last;
        }
    }
    
    $base_ec = join('.', @new);
    
    my $qry = qq{
        SELECT cvterm.cvterm_id
          FROM cvterm
          JOIN cvterm_dbxref cvdbx ON cvterm.cvterm_id = cvdbx.cvterm_id
          JOIN dbxref ON cvdbx.dbxref_id = dbxref.dbxref_id
         WHERE cvterm.cv_id = $cv{EC}
           AND dbxref.accession = ?
    };
    my $dsh = $dbh->prepare($qry);
       $dsh->execute( $base_ec );

    my $base_cvterm_id;

    while ( my $row = $dsh->fetchrow_hashref ) {
        $base_cvterm_id = $$row{cvterm_id};
        last;
    }
    
    if (defined $base_cvterm_id ) {
        return $base_cvterm_id;
    } else {
        die( "Unable to find cvterm_id corresponding to base accession $base_ec.  Check the base term?" );
    }
}

sub get_next_ec_id_num {
    my $next_num = 1;
    
    my $qry = qq{
        SELECT max(accession) AS accession
          FROM dbxref
         WHERE db_id = $db{EC}
           AND accession LIKE "EC:\%"
    };
    my $accession_selector = $dbh->prepare($qry);
    $accession_selector->execute();
    
    while (my $row = $accession_selector->fetchrow_hashref ) {
        if ( $$row{accession} =~ /^EC\:(\d+)/ ) {
            $next_num = $1 + 1;
        } else {
            die "expected accession to look like EC:00005555, instead got: $$row{accession}\n";
        }
    }
    
    $accession_selector->finish();
    
    return $next_num;
}

sub _log {
    my $msg = shift;

    print "$msg\n";
    print $logfh "$msg\n" if $logfh;
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( database user password ec_number name );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## handle some defaults
    $$options{server} = 'localhost' unless defined $$options{server};
    $$options{database_type} = 'mysql' unless defined $$options{database_type};
    $$options{definition} = '' unless defined $$options{definition};

    ## database type must be either mysql or postgresql
    if ( $$options{database_type} ne 'mysql' && $$options{database_type} ne 'postgresql' ) {
        die "value for option --database_type must be either 'mysql' or 'postgresql'";
    }
}

sub get_cvterm_id {
    my $name = shift;
    my $qry = qq{
        SELECT cvterm_id
        FROM cvterm
        WHERE name = ?
    };
    
    my $cvterm_selector = $dbh->prepare($qry);
       $cvterm_selector->execute($name);
    
    my $cvterm_id = ( $cvterm_selector->fetchrow_array )[0];
    $cvterm_selector->finish();

    if ( $cvterm_id ) {
        _log("INFO: got cvterm_id $cvterm_id for name $name");
        return $cvterm_id;
    } else {
        _log("ERROR: failed to retrieve cvterm_id for name $name");
        die "ERROR: failed to retrieve cvterm_id for name $name\n";
    }
}

sub get_db_id {
    my $name = shift;
    my $qry = qq{
        SELECT db_id
        FROM db
        WHERE name = ?
    };
    
    my $db_selector = $dbh->prepare($qry);
       $db_selector->execute($name);
    
    my $db_id = ( $db_selector->fetchrow_array )[0];
    $db_selector->finish();

    if ( $db_id ) {
        _log("INFO: got db_id $db_id for name $name");
        return $db_id;
    } else {
        _log("ERROR: failed to retrieve db_id for name $name");
        die "ERROR: failed to retrieve db_id for name $name\n";
    }
}


sub get_cv_id {
    my $name = shift;
    my $qry = qq{
        SELECT cv_id
        FROM cv
        WHERE name = ?
    };
    
    my $cv_selector = $dbh->prepare($qry);
       $cv_selector->execute($name);
    
    my $cv_id = ( $cv_selector->fetchrow_array )[0];
    $cv_selector->finish();

    if ( $cv_id ) {
        _log("INFO: got cv_id $cv_id for name $name");
        return $cv_id;
    } else {
        _log("ERROR: failed to retrieve cv_id for name $name");
        die "ERROR: failed to retrieve cv_id for name $name\n";
    }
}

## gets the next ID in a chado table with canonical names.  table called $foo and id called $foo_id
sub get_next_canonical_id {
    my $table_name = shift;
    
    my $last_id;
    
    my $qry = qq{
        SELECT max(${table_name}_id)
          FROM $table_name
    };
    
    my $max_selector = $dbh->prepare($qry);
       $max_selector->execute();
    
    $last_id = ( $max_selector->fetchrow_array )[0];
           
    $max_selector->finish();
    
    ## the table must exist else the query will fail.  no need to check for that.
    ## if ID was set there was at least one row in the table.
    if ( $last_id ) {
        ## rows in the table - add 1
        return $last_id + 1;
    } else {
        ## no rows in the table, return 1
        return 1;
    }
}
