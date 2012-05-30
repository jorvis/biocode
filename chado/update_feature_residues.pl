#!/usr/bin/env perl

=head1 NAME

update_feature_residues.pl - For a given feature type, this script reads a Chado instance
and updates the feature.residues field for all instances of that type based on the molecule
it's located on.

=head1 SYNOPSIS

USAGE: update_feature_residues.pl 
            --database=eha3
            --user=someuser
            --password=somepass
            --feature_type=CDS
          [ --server=localhost
            --molecule_type=assembly
            --log=/path/to/some.log ]


=head1 OPTIONS

B<--database,-d>
    MySQL database name to connect to.

B<--user,-u>
    User account with select privileges on the specified database.

B<--password,-p>
    Password for user account specified.

B<--feature_type,-f>
    Feature type to be updated.  The only currently-supported (tested) types are CDS and polypeptide.

B<--molecule_type,-m>
    Optional.  Type of molecule to consider.  For most uses, the default 'assembly' is fine.

B<--server,-s>
    Optional.  MySQL database server to connect to (default = localhost).

B<--log,-l> 
    Optional.  Full path to a log file to create.

B<--help,-h>
    This help message


=head1  DESCRIPTION

This script creates the SQL necessary to update the feature.residues field for all features 
of a given type.  Currently limited to MySQL Chado instances.

=head1  INPUT

If the feature in question isn't located on anything it won't be updated.  Only those features
where feature.is_analysis = 0 will be included.

=head1  OUTPUT

This script doesn't actually modify the database.  Instead, it writes an output of
the needed SQL statements to perform the update.  If you want to create a file of the SQL statements
you should redirect to a file.  Logging information will still be printed to the screen.

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use Bio::DB::SwissProt;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'database|d=s',
                          'server|s=s',
                          'user|u=s',
                          'password|p=s',
                          'feature_type|f=s',
                          'molecule_type|m=s',
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

_log("attempting to create database connection");
my $dsn = "dbi:mysql:database=$options{database};host=$options{server}";
my $dbh = DBI->connect($dsn, $options{user}, $options{password}, {PrintError=>1, RaiseError=>1} );

my $feature_type_id  = get_cvterm_id( $options{feature_type} );
my $molecule_type_id = get_cvterm_id( $options{molecule_type} );

## structure like:
#   $h{$feature_id} = { on => mol_id, fmin => N, fmax => N, strand => 1|-1 }
my %features;

## structure like:
#   $h{molecule_feature_id} = TGATATA ....
my %molecules;

## get all features of this type located on anything.
my $qry = qq{
    SELECT feat.feature_id AS feat_id, feat.uniquename as feat_name, 
           ftl.fmin, ftl.fmax, ftl.strand, mol.feature_id AS mol_id, mol.uniquename AS mol_name
      FROM feature feat
        JOIN featureloc ftl ON feat.feature_id = ftl.feature_id
        JOIN feature mol ON ftl.srcfeature_id = mol.feature_id
     WHERE feat.type_id = $feature_type_id
       AND mol.type_id = $molecule_type_id
       AND feat.is_obsolete = 0
       AND feat.is_analysis = 0
};

my $feature_coord_selector = $dbh->prepare($qry);
   $feature_coord_selector->execute();

while ( my $row = $feature_coord_selector->fetchrow_hashref ) {
    _log( "INFO: $$row{feat_name} located on $$row{mol_name} at coordinates $$row{fmin} - $$row{fmax}, strand $$row{strand}" );
    $features{ $$row{feat_id} } = { on => $$row{mol_id}, fmin => $$row{fmin}, fmax => $$row{fmax}, strand => $$row{strand} };
    $molecules{ $$row{mol_id} } = '';
}



$feature_coord_selector->finish();

## now get the sequences for all the molecules.
for my $molecule ( keys %molecules ) {
    $molecules{$molecule} = get_molecule_seq( $molecule );
    _log( "INFO: Got sequence of length " . length($molecules{$molecule}) . " for molecule $molecule"  );
}

for my $feature ( keys %features ) {
    my $feat_seq = substr( $molecules{ $features{$feature}{on} }, $features{$feature}{fmin}, $features{$feature}{fmax} - $features{$feature}{fmin} );
    
    if ( $features{$feature}{strand} eq -1 ) {
        $feat_seq = reverse $feat_seq;
        $feat_seq =~ tr|atgcATGC|tacgTACG|;
    }
    
    ## we could now translate it here if that's appropriate for the feature type
    if ( $options{feature_type} eq 'polypeptide' ) {
        $feat_seq = translate_sequence( $feat_seq, $feature );
    }
    
    print "UPDATE feature SET residues = '$feat_seq' WHERE feature_id = $feature;\n";
}

$dbh->disconnect();

exit(0);

sub get_molecule_seq {
    my $feature_id = shift;
    
    my $seq;
    
    my $qry = qq{
        SELECT residues 
          FROM feature
         WHERE feature_id = ?
    };
    my $seq_selector = $dbh->prepare($qry);
       $seq_selector->execute($feature_id);
       
    while ( my $row = $seq_selector->fetchrow_hashref ) {
        $seq = $$row{residues};
        $seq =~ s/\s//g;
        last;
    }
    
    $seq_selector->finish;
    
    if ( $seq ) {
        return $seq;
    } else {
        die "failed to pull sequence residues for feature ID $feature_id\n";
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

sub _log {
    my $msg = shift;
    print STDERR "$msg\n";

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
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    if ( $$options{feature_type} ne 'CDS' &&
         $$options{feature_type} ne 'polypeptide' ) {
        die "the only currently-supported feature_types are CDS and polypeptide\n";
    }
    
    
    ## handle some defaults
    $$options{server} = 'localhost' unless defined $$options{server};
    $$options{molecule_type} = 'assembly' unless defined $$options{molecule_type};
}


sub translate_sequence {
    my ($dna, $label) = @_;
    
    my $polypeptide_seq = '';
    
    for(my $i=0; $i<(length($dna) - 2); $i+=3) {
        $polypeptide_seq .= codon2aa( substr($dna,$i,3) );
    }
    
    ## report if this had an internal stop
    if ( $polypeptide_seq =~ /\*.*?.$/ ) {
        _log("ERROR: internal stops detected in feature with feature_id: $label: ($polypeptide_seq)");
    }
    
    return $polypeptide_seq;
}

## taken from Tisdall's perl book, chapter 8.
#   modified to use * for stops instead of _
sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if (exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{
        return 'X';
    }
}
