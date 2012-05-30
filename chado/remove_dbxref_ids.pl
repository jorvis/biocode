#!/usr/bin/env perl

=head1 NAME

remove_dbxref_ids.pl - removes dbxref IDs (such as locus accessions) within a chado database

=head1 SYNOPSIS

USAGE: assign_dbxref_ids.pl
            --database=abs2
            --user=someuser
            --password=somepass
            --db_name=IGS
            --server=yourdbhost
          [ --organism_id
            --version
            --log=/path/to/some.log
            --help ]


=head1 OPTIONS

B<--database,-d>
    Database name to connect to.

B<--user,-u>
    User account with select, insert, and update privileges on the specified database.

B<--password,-p>
    Password for user account specified.

B<--db_name,-b>
    Name of database (corresponds to db.name)

B<--server,-s>
    Name of server where database is located.

B<--organism_id,-o>
    Optional.  Only process a specific organism_id (default to all organisms).

B<--version,-v>
    Optional.  Enter a value  to represent the version of the ID set.  This will be inserted
    into the dbxref.version field for all entries. (default all versions with the db name)

B<--log,-l> 
    Optional.  Full path to a log file to create.

B<--help,-h>
    This help message

=head1  DESCRIPTION

    Removes dbxrefs based on a database name.  Will remove rows from the feature_dbxref, 
    dbxref and db table based on passed in parameters.

    Rows are removed from the feature_dbxref and dbxref tables based on the organism_id
    and version information passed in.  It is assumed that if the entry for feature_dbxref
    table is removed, the corresponding dbxref should also be removed (and this is how the 
    script operates).  

    The db is only removed if there are no dbxref rows that reference that db.
    
=head1  INPUT
    
    The only input is the name of the database, db.name and database connection information.

=head1  OUTPUT

    This script will produce no output to the terminal unless an error occurs.  Check the log file
    for script run information.  If not log file is specified no output is produced.

=head1  CONTACT

    Kevin Galens
    kevingalens@gmail.com

=cut

use strict;
use warnings;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;
$|++;

## GLOBALS ########
my %options = ();
my $logfh;
###################

my $results = GetOptions (\%options, 
                          'database|d=s',
                          'user|u=s',
                          'password|p=s',
                          'db_name|b=s',
                          'server|s=s',
                          'organism_id|o=s',
                          'version|v=s',
                          'log|l=s',
                          'help|h') || pod2usage();


&check_options( \%options );
my $dbh = &connect( $options{'database'}, $options{'server'}, 
                    $options{'user'}, $options{'password'} );
die("Could not connect to database $options{'database'}")
    unless( $dbh );
_log("Connected to database $options{'database'}");

my $db_id = &get_db_id( $options{'db_name'} );
die("There was no db.name specified in database $options{'database'} with name $options{'db_name'}")
    unless( defined( $db_id ) );
_log("Database $options{'db_name'} has an id of $db_id");

my $count_remove_feat_dx = &remove_rows_from_fd_dx( $db_id, $options{'organism_id'},
                                                    $options{'version'} );
_log("Removed $count_remove_feat_dx rows from feature_dbxref");

my $count_db_removed = &remove_db( $db_id );

$dbh->disconnect;

sub get_db_id {
    my ($db_name) = @_;
    
    my $query = "SELECT db_id FROM db ".
        "WHERE name = '$db_name'";
    my $sth = $dbh->prepare( $query );
    $sth->execute();

    my $results = $sth->fetchall_arrayref();
    
    my $retval;
    
    if( @{$results} == 1 and @{$results->[0]} == 1 ) {
        $retval = $results->[0]->[0];
    }
    
    return $retval;
    
}

sub remove_rows_from_fd_dx {
    my ( $db_id, $org_id, $vers ) = @_;

    my $query = "DELETE fd, dx ".
        "FROM feature_dbxref fd, dbxref dx";

    if( $org_id ) {
        $query .= ", feature f";
    }

    $query .= " WHERE dx.dbxref_id = fd.dbxref_id ".
        "AND dx.db_id = $db_id ";

    if( $org_id ) {
        $query .= "AND f.feature_id = fd.feature_id ".
            "AND f.organism_id = $org_id ";
    }

    if( $vers ) {
        $query .= "AND dx.version = $vers";
    }

    my $sth = $dbh->prepare( $query );

    eval {
        $sth->execute();
    };
    if( $@ ) {
        die("Could not execute query: $DBI::errstr");
    }
    
    return scalar($sth->rows);
    
        
}

sub remove_db {
    my ($db_id) = @_;
    
    my $query = "SELECT count(*) FROM dbxref WHERE db_id = $db_id";
    my $sth = $dbh->prepare( $query );
    $sth->execute();
    my $results = $sth->fetchall_arrayref();

    my $count = $results->[0]->[0];
    
    if( $count == 0 ) {
        my $delete = "DELETE FROM db where db_id = $db_id";
        my $dsh = $dbh->prepare( $delete );
        eval {
            $dsh->execute();
        };
        if( $@ ) {
            die("Could not remove database $db_id ($options{'db_name'}) from database: $DBI::errstr");
        }
        _log("Removed db $db_id from database $options{'database'}");
    } else {
        _log("Did not remove db $db_id from database because there are still $count dbxrefs referencing the id");
    }

    return ($count) ? 0 : 1;
    
}

sub connect {
    my ($database, $server, $user, $password) = @_;

    my $dbh = DBI->connect("dbi:mysql:database=$database;host=$server", 
                           $user, $password, 
                           {PrintError=>1, RaiseError=>1});
    return $dbh;

}

sub check_options {
    my ($opts) = @_;
    
    my @reqs = qw(database user password db_name server);
    foreach my $req ( @reqs ) {
        die("Option --$req is required") unless( defined( $opts->{$req} ) );
    }
}

sub _log {
    my ($msg) = @_;
    print $logfh "$msg\n" if( $logfh );
}

