#!/usr/bin/env perl

=head1 NAME

create_rotating_mysql_backup.pl - Backs up a MySQL database with N-count historical copies

=head1 SYNOPSIS

USAGE: create_rotating_mysql_backup.pl 
            --host=db_server 
            --database=yourdb
            --user=you 
            --password=whatever
            --backup_dir=/some/path
          [ --backup_count=10 ]

=head1 OPTIONS

B<--host,-s 
    Database server hostname.

B<--database,-d>
    MySQL database name.

B<--user,-u>
    MySQL user name with permissions to run the mysqldump utilities.

B<--password,-p>
    MySQL password for the user above.
    
B<--backup_dir,-b>
    Directory where the backup files are kept.    

B<--backup_count,-n>
    Optional.  How many rotating backups to keep (Default = 10)  

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Appropriate to put on a cron, this script runs the mysqldump command and,
given an output directory, maintains N-number rotational backups of a
given database.

=head1  INPUT

Any MySQL database and user credentials to run the mysqldump utility on it.

=head1  OUTPUT

A series of MySQL dump files named like:

    $host.$dbname.$datestamp.sql.dump

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
              'host|s=s',
              'database|d=s',
              'user|u=s',
              'password|p=s',
              'backup_dir|b=s',
              'backup_count|n=i',
              'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## play nicely, or not
#umask(0000);

## make sure all passed options are peachy
&check_parameters(\%options);


## quicker way than this to get YYYY-MM-DD?
my ($year, $mon, $day) = (localtime time)[5,4,3];
my $datestamp = sprintf("%d-%02d-%02d", $year + 1900, ++$mon, $day);

my $cmd = "mysqldump -h $options{host} -u $options{user} --password=$options{password} $options{database} > $options{backup_dir}/$options{host}.$options{database}.$datestamp.sql.dump";

print "$cmd\n";
`$cmd`;

## now list the databases in that directory and only keep the newest 10
my @backups;
for (glob "$options{backup_dir}/$options{host}.$options{database}.*.sql.dump") {
    push @backups, $_;
}

my $backup_count = scalar @backups;
if ( $backup_count > $options{backup_count} ) {
    for (my $i=0; $i<$options{backup_count}; $i++) {
        print "removing $backups[$i]\n";
        unlink($backups[$i]);
    }
}


sub check_parameters {
    unless ( defined $options{host} &&
             defined $options{database} && 
             defined $options{user} && 
             defined $options{password} &&
             defined $options{backup_dir}) {
             
        pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
    }
    
    $options{backup_count} = 10 if ! defined $options{backup_count};
}

