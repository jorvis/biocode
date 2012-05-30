#!/usr/bin/env perl

=head1 NAME

assign_dbxref_ids.pl - creates dbxref IDs (such as locus accessions) within a chado database

=head1 SYNOPSIS

USAGE: assign_dbxref_ids.pl
            --database=abs2
            --user=someuser
            --password=somepass
            --db_name=IGS
          [ --id_prefix=ABS1_
            --id_prefix_map=/path/to/mapfile.txt
            --feat_type=gene
            --server=SYBTIGR
            --organism_id=all
            --assign_missing=1
            --version=locus
            --zero_pad=6
            --log=/path/to/some.log ]


=head1 OPTIONS

B<--database,-d>
    Sybase database name to connect to.

B<--user,-u>
    User account with select, insert, and update privileges on the specified database.

B<--password,-p>
    Password for user account specified.

B<--db_name,-b>
    This corresponds to the db.name field for ID reference source.

B<--id_prefix,-i>
    Optional. Prefix for the IDs to be generated.  This will be followed by an incrementing number so that
    if you pass 'ABS1_' the IDs generated will be like 'ABS1_000001' and so on. Either --id_prefix or 
    --id_prefix_map needs to be used.

B<--id_prefix_map,-I>
    Optional. A map file containing molecule uniquename and id_prefixes.  Two column, tab-delimited file. If
    a molecule exists in the database and is not in the mapping file, script will use the value of the --id_prefix
    option for prefix.  If the --id_prefix option was not given, the script will fail.

B<--feat_type,-t>
    Optional.  Feature type on which to assign dbxref IDs (default = gene).

B<--server,-s>
    Optional.  Sybase server to connect to (default = localhost).

B<--organism_id,-o>
    Optional.  Only process a specific organism_id (default = all).

B<--assign_missing,-a>
    Optional. Will search for the db name, version, prefix given and only assign features
    which don't have an entry for that db.

B<--version,-v>
    Optional.  Enter a value to represent the version of the ID set.  This will be inserted
    into the dbxref.version field for all entries. (default = locus)

B<--zero_pad,-z>
    Optional.  Defines number of digits IDs will be zero-padded (default = 6).

B<--molecule_type,-m>
    Optional.  Accessions are assigned sequentially along molecules of this passed type.  Features
    assigned to multiple molecules of the same type aren't allowed (default = assembly). Can also be
    a comma separated list of molecule types.  Ex. assembly,plasmid,contig

B<--log,-l> 
    Optional.  Full path to a log file to create.

B<--help,-h>
    This help message

=head1  DESCRIPTION

There are occasions where an active project needs to have identifiers (such as locus names)
for software display purposes or tracking.  This script can be used to create actual locus
identifiers or internal identifers.

Within chado, both internal and locus ID references are stored in a common way:

    feature.feature_id -> feature_dbxref.feature_id -> 
    feature_dbxref.dbxref_id -> dbxref.dbxref_id -> dbxref.db_id -> db.db_id

Entering an id_prefix_map will give the user the ability to assign locus identifiers with 
different prefixes for different molecules in the same database.  Each prefix gets its own
set of numbers (each prefix represents its own namespace).
    
=head1  INPUT

The only input required is the information necessary to connect to the database and information
about the IDs to create.  If the db_name passed doesn't exist in the chado 'db' table it will
be created.  

The id_prefix_map should be in the following format:
    prefix1_    db.assembly_uniquename.1
    prefix2_    db.assembly_uniquename.2
    ...

Currently assumes an input database in MySQL of the chado schema.  Will first check (and fail)
if there are already entries in the database with the passed db_name, feat_type and version.

Considers only those features of the passed feature type (default transcript) where 
feature.is_obsolete = 0

Because this script assigns dbxrefs sequentially along each molecule (top strand first) it doesn't
support features located on multiple assemblies.

=head1  OUTPUT

This script will produce progress information on STDERR, and will write a log of all
statements performed if the --log option is used.

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
                          'db_name|b=s',
                          'id_prefix|i=s',
                          'id_prefix_map|I=s',
                          'feat_type|t=s',
                          'server|s=s',
                          'organism_id|o=i',
                          'version|v=s',
                          'zero_pad|z=i',
                          'assign_missing|a=s',
                          'molecule_type|m=s',
                          'log|l=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## toggling this will disable all inserts
my $debugging = 0;

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

_log("attempting to create database connection");
my $dbh = DBI->connect("dbi:mysql:database=$options{database};host=$options{server}", 
                        $options{user}, $options{password}, 
                        {PrintError=>1, RaiseError=>1});

#parse the molecule type
my @molecule_types = split(/[\s,]+/,$options{molecule_type});

#parse id_prefix_map file (if passed)
my $prefix_lookup = &parse_id_prefix_map_file( $options{id_prefix_map} ) if( $options{id_prefix_map} );

#if id_prefix was specified, use it as prefix default
if( $options{id_prefix} ) {
    $prefix_lookup->{default} = $options{id_prefix};
} else {
    #if an id_prefix was not passed, make sure all molecules in database are present
    #in the id_prefix_map file.
    &check_id_prefix_map( $prefix_lookup, \@molecule_types );
}

## will be reused throughout
my $qry;

## to get the ordering of accessions correct we need to verify that no features are
##  located on more than one molecule of the same type.
$qry = qq{
    SELECT f.uniquename, count(srcf.uniquename)
      FROM feature srcf, feature f, featureloc fl, cvterm srccv
     WHERE fl.srcfeature_id = srcf.feature_id
       AND fl.feature_id = f.feature_id
       AND srcf.type_id = srccv.cvterm_id
       AND f.is_obsolete = 0
       AND srccv.name = ?
  GROUP BY f.uniquename
    HAVING count(srcf.uniquename) > 1
};
my $multiloc_checker = $dbh->prepare($qry);

foreach my $molecule_type ( @molecule_types ) {
    $multiloc_checker->execute( $molecule_type );

    my $multiloc_count = 0;

    while ( my $feat = $multiloc_checker->fetchrow_hashref ) {
        $multiloc_count++;
        _log("ERROR: feature ($$feat{uniquename} found on multiple molecules of type $molecule_type");
    }

    if ( $multiloc_count ) {
        die "ERROR: $multiloc_count features found located on multiple molecules of the same type ".
            "[$molecule_type].  See log for details.";
    }
}
$multiloc_checker->finish();

## which feature type are we handling
my $feature_type_cvterm_id = get_cvterm_id( $options{feat_type} );

## now we need to check if the db of the passed name exists at all.
## if not defined in this call, it will be set when the database is created.
my $db_id = &get_db_id( $options{db_name} );


## CHANGES ##
if ( defined $db_id ) {
    _log("INFO: db with name $options{db_name} and id $db_id found");
    
    ## next check that no accessions already exist of the passed db name and version 
    #   on the feature type specified by the user
    $qry = qq{
        SELECT COUNT(dbx.dbxref_id)
          FROM feature_dbxref fdbx
                JOIN feature f ON fdbx.feature_id = f.feature_id
                JOIN dbxref dbx ON fdbx.dbxref_id = dbx.dbxref_id
         WHERE dbx.db_id = ?
           AND f.type_id = ?
           AND dbx.version = ?
    };
    my $dbxref_selector = $dbh->prepare($qry);
       $dbxref_selector->execute( $db_id, $feature_type_cvterm_id, $options{version} );

    my $dbxref_count = ( $dbxref_selector->fetchrow_array )[0];
    
    $dbxref_selector->finish();

        

    if ( defined $dbxref_count && $dbxref_count > 0 ) {
        unless( $options{'assign_missing'} ) {
            _log("ERROR: found existing $options{feat_type} dbxref entries for db_id $db_id and version $options{version}");
            die("ERROR: found existing $options{feat_type} dbxref entries for db_id $db_id and version $options{version}");
        }
    } else {
        if( $options{'assign_missing'} ) {
            my $msg = "ERROR: found no existing $options{feat_type} dbxref entries for db_id $db_id and version $options{version}. ".
                "There should exist some ids already when using the -assign_missing option";
            _log($msg);
            die($msg);
        } else {
            _log("INFO: no dbxref entries found for db_id $db_id and version $options{version}");
        }
    }

} else {
    _log("INFO: db with name $options{db_name} not found");
    
    if( defined $options{'assign_missing'} && $options{'assign_missing'} ) {
        my $msg = "ERROR: db name $options{'db_name'} does not exist in database. DB should exist when using ".
            "the --assign_missing option";
        _log($msg);
        die($msg);
    }

    ## the passed db doesn't exist.  create it.
    my $next_db_id = get_next_canonical_id( 'db' );
    
    $qry = qq{
        INSERT INTO db ( db_id, name )
             VALUES ( ?, ? );
    };
    
    _log("INFO: inserting new db with id $next_db_id and name $options{db_name}");
    
    my $db_inserter = $dbh->prepare($qry);
       $db_inserter->execute( $next_db_id, $options{db_name} ) unless $debugging;
    
    ## now set the db_id of the one we just entered.
    $db_id = $next_db_id;
}

#print "gathering features\n";
#######################
## feature gathering ##
#######################

## we know the db entry is OK to this point, so now we'll go through the features
my @feat_options = (\@molecule_types, $feature_type_cvterm_id, $options{'organism_id'});
if( $options{'assign_missing'} ) {
    push(@feat_options, ($options{'db_name'}, $options{'version'}) );
}
my @feature_ids_to_tag = &get_features( @feat_options );

if ( scalar @feature_ids_to_tag ) {
    _log("INFO: found " . scalar(@feature_ids_to_tag) . " features of type $options{feat_type} ");
} else {
    _log("WARN: no features found of type ");
}

########################
## feature processing ##
########################

my $next_dbxref_id = get_next_canonical_id( 'dbxref' );
my $next_feature_dbxref_id = get_next_canonical_id( 'feature_dbxref' );

$qry = qq{
    INSERT INTO dbxref ( dbxref_id, db_id, accession, version )
    VALUES ( ?, ?, ?, ? )
};
my $dbxref_inserter = $dbh->prepare($qry);

$qry = qq{
    INSERT INTO feature_dbxref ( feature_dbxref_id, feature_id, dbxref_id, is_current )
    VALUES ( ?, ?, ?, 1 )
};
my $feature_dbxref_inserter = $dbh->prepare($qry);

##initialize the id_num_part hash.  This contains a lookup of the next id
##for locus identifiers based on molecule name
my $id_num_part = {};
# map { $id_num_part->{$_} = 1 } keys %{$prefix_lookup};

for my $prefix ( keys %{$prefix_lookup} ) {
    $id_num_part->{$prefix} = get_next_accession_by_prefix( $db_id, $$prefix_lookup{$prefix} );
}

my $accession = '';
my $dbxref_id;
my $feature_dbxref_id;

for my $feat_info ( @feature_ids_to_tag ) {
    my ($feat_id, $mol_uniquename) = @{$feat_info};

    #search for the molecule in the prefix_lookup hash.  If it's not there then it wasn't
    #in the prefix_map file and default should be used.  If there is not default (--id_prefix
    # option) then die.  
    my $id_prefix;
    if( exists( $prefix_lookup->{$mol_uniquename} ) ) {
        $id_prefix = $prefix_lookup->{$mol_uniquename};
    } elsif( exists( $prefix_lookup->{'default'} ) ) {
        $id_prefix = $prefix_lookup->{'default'};
        $mol_uniquename = 'default';
    } else {
        die("$mol_uniquename was not in id_prefix_map (or no prefix map was passed) and there was ".
            "no default id_prefix name.");
    }

    $accession = $id_prefix . sprintf("%0${options{zero_pad}}d", $id_num_part->{$mol_uniquename}++);
    $dbxref_id = $next_dbxref_id++;
    $feature_dbxref_id = $next_feature_dbxref_id++;

    _log("INFO: inserting dbxref entry ($dbxref_id) for accession $accession, version $options{version}");
    $dbxref_inserter->execute( $dbxref_id, $db_id, $accession, $options{version} ) unless $debugging;

    _log("INFO: assigning dbxref_id $dbxref_id, accession $accession to feature_id $feat_id (feature_dbxref_id: $feature_dbxref_id)");
    $feature_dbxref_inserter->execute( $feature_dbxref_id, $feat_id, $dbxref_id ) unless $debugging;
    
}

$feature_dbxref_inserter->finish();
$dbxref_inserter->finish();

$dbh->disconnect();

exit(0);

sub get_features {
    my ($mol_types, $feat_type_id, $organism_id, $db_name, $version ) = @_;
    my @feature_ids_to_tag;
    
    _log("INFO: getting list of features to be tagged with dbxrefs");

    ## this is the source molecule(s), assembly by default
    my @molecule_cvterm_ids;
    foreach my $mol_type ( @{$mol_types} ) {
        push( @molecule_cvterm_ids, get_cvterm_id( $mol_type ) );
    }
    my $mol_cvterm_id_string = join( ",", @molecule_cvterm_ids );

    ## the ordering clause is added after the conditional below
    my $qry = qq{
        SELECT f.feature_id, f.uniquename, molecule.uniquename as mol, fl.fmin, fl.fmax
          FROM feature f
          JOIN featureloc fl ON fl.feature_id = f.feature_id
          JOIN feature molecule ON molecule.feature_id = fl.srcfeature_id
         WHERE f.type_id = ?
           AND molecule.type_id IN ( $mol_cvterm_id_string )
           AND f.is_obsolete = 0
        };
    my @feature_selector_args = ( $feat_type_id );

    ## if a db name and version is supplied, only return
    ## those features which DON'T have a dbxref for the
    ## specified db_name and version
    if( defined $db_name && defined $version ) {
        $qry = 
            "SELECT f.feature_id, f.uniquename, molecule.uniquename as mol, count(fd.feature_dbxref_id) as count ".
            "FROM feature f LEFT OUTER JOIN (feature_dbxref fd, dbxref d, db) ".
            "ON (f.feature_id = fd.feature_id AND fd.dbxref_id = d.dbxref_id AND ".
            "db.db_id = d.db_id AND db.name = ? AND d.version = ?), feature molecule, featureloc fl ".
            "WHERE f.type_id = ? AND f.is_obsolete = 0 AND f.feature_id = fl.feature_id ".
            "AND fl.srcfeature_id = molecule.feature_id AND molecule.type_id IN ($mol_cvterm_id_string)";
        @feature_selector_args = ( $db_name, $version, $feat_type_id );
    }

    ## are we also filtering based on organism_id?
    if ( $organism_id ) {
        _log("INFO: filtering feature list to include only organism ID: $organism_id");
        $qry .= " AND f.organism_id = ? ";
        push @feature_selector_args, $organism_id;
    }

    
    $qry .= " GROUP BY f.uniquename HAVING count = 0" if( defined $db_name && defined $version );
    $qry .= " ORDER BY molecule.uniquename, fl.fmin";

    my $feature_selector = $dbh->prepare($qry);
    $feature_selector->execute( @feature_selector_args );

    while ( my $feature = $feature_selector->fetchrow_hashref ) {
        push @feature_ids_to_tag, [$feature->{feature_id},$feature->{mol}];
    }

    $feature_selector->finish();
    @feature_ids_to_tag;
}

sub get_next_accession_by_prefix {
    my ($db_id, $prefix) = @_;
    
    my $qry = qq{
        SELECT accession
          FROM dbxref
         WHERE db_id = ?
           AND version = ?
    };
    
    my $selector = $dbh->prepare($qry);
       $selector->execute($db_id, $options{version});
       
    my $last_id = 0;
    
    while ( my $row = $selector->fetchrow_hashref ) {
        if ( $$row{accession} =~ /^${prefix}(\d+)/ ) {
            my $next_num = $1;
            
            ## make sure it's an integer
            if ( $next_num !~ /^\d+$/ ) {
                die("Got ($next_num) when expecting an integer as numeric portion of accession for prefix $prefix in db $db_id");
            }
            
            $last_id = $next_num if $next_num > $last_id;
        }
    }
    
    $selector->finish();
    
    return $last_id + 1;
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


## gets the numerical db id for a passed db.name.  returns undef if not found.
sub get_db_id {
    my $db_name = shift;

    my $db_id;

    ## name is a unique key, so we don't have worry about multiples
    my $qry = qq{
        SELECT db_id 
          FROM db
         WHERE name = ?
    };
    my $db_selector = $dbh->prepare($qry);
       $db_selector->execute( $db_name );
    
    while ( my $db = $db_selector->fetchrow_hashref ) {
        $db_id = $db->{db_id};
        last;
    }
    
    $db_selector->finish();
    
    return $db_id;
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

## parse an input prefix_id_map file
## format- 2 column tab delimited. 1st column prefix_id, 2nd column molecule name
## ex.
## ABC    abc.assembly.1
## ABCp   abcp.plasmid.1
sub parse_id_prefix_map_file {
    my ($map_file) = @_;
    my $retval;
    open(IN, "< $map_file" ) or die("Unable to open map file: $map_file [$!]");
    while( <IN> ) {
        chomp;
        my ($prefix,$molecule_name) = split(/\s+/);
        $retval->{$molecule_name} = $prefix;
    }
    close(IN);
    return $retval;
}

sub check_id_prefix_map {
    my ($prefix_lookup, $molecule_type ) = @_;
    my @db_molecule_names;
    
    #select all feature.uniquenames from database that have 
    #one of the specified molecule types
    my $query = 
        "SELECT uniquename from feature f ".
        "WHERE f.type_id = ? ".
        "AND is_obsolete = 0";

    my $molecule_query = $dbh->prepare( $query );
    
    foreach my $mt ( @{$molecule_type} ) {
        my $cvterm_id = &get_cvterm_id( $mt );
        $molecule_query->execute( $cvterm_id );
        while( my $results = $molecule_query->fetchrow_arrayref() ) {
            push(@db_molecule_names, $results->[0] );
        }
    }
    $molecule_query->finish();

    foreach my $molecule_name ( @db_molecule_names ) {
        unless( exists( $prefix_lookup->{$molecule_name} ) ) {
            _log("Molecule $molecule_name was found in database but not in ".
                 "id_prefix_map [$options{id_prefix_map}]");
            die("Molecule $molecule_name was found in database but not in ".
                "id_prefix_map [$options{id_prefix_map}]");
        } 
    }
    
    return 1;
    
}

sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( database user password db_name );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }

    # either id_prefix of id_prefix_map must be supplied
    my $flag = 0;
    if( $$options{id_prefix} ) {
        $$prefix_lookup{default} = $$options{id_prefix};
        $flag++;
    } elsif( $$options{id_prefix_map} ) {
        $flag++;
    }
    unless( $flag ) {
        die("Either --id_prefix or --id_prefix_map (or both) is required");
    }
    
    ## handle some defaults
    $$options{server} = 'localhost' unless defined $$options{server};
    $$options{version} = 'locus' unless defined $$options{version};
    $$options{feat_type} = 'gene' unless defined $$options{feat_type};
    $$options{zero_pad} = 6 unless defined $$options{zero_pad};
    $$options{molecule_type} = 'assembly' unless defined $$options{molecule_type};
    
    ## reset in case the user actually passes 'all'
    $$options{organism_id} = undef if $$options{organism_id} eq 'all';
}

