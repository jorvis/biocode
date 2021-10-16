#!/usr/bin/env perl

=head1 NAME

filter_columnar_file.pl - Filters a tab delimited text file on as many columnar values as you pass.

=head1 SYNOPSIS

    USAGE: filter_columnar_file.pl somefile.tab 0 foo 13 bar|dak

This would read the somefile.tab and return only those rows where column 1
matches the regex 'foo' and column 14 matches either 'bar' or 'dak'

=head1 OPTIONS

This script takes positional arguments.  First pass an tab-delimited input file path,
then any number of pairs of column number (0-based) and value.

=head1  DESCRIPTION

It's common to want to filter the rows of a tab-delimited file based on the values of 
one or more columns but still print all the columns of each row.  This script lets you 
filter based on any number of columns as well as give alternate values for each column.

=head1  INPUT

The input is simply a tab-delimited file and then any number of column number / value
pairs.  The column numbers start at 0.  The value arguments can be separated by pipes
to specify multiple options (see example above.)

=head1  OUTPUT

All columns for any matching rows are printed to STDOUT.

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use warnings;

if ( scalar(@ARGV) < 3 ) {
    print STDERR "\nERROR.  You must provide at least 3 arguments.  See the perldoc.\n\n";
    exit(1);
}

my $tab_file = shift @ARGV;
my %patterns = ();
my $max_column_number = 0;

## load any argument pairs
while ( scalar @ARGV >= 2 ) {
    my $col = shift @ARGV;
    my $value = shift @ARGV;
    
    ## col should look like a number
    if ( $col !~ /^\d+$/ ) {
        die "Option parsing error.  Expected ($col) to be an integer describing a column.";
    }
    
    $max_column_number = $col if $col > $max_column_number;
    my @regexes = split('\|', $value);
    $patterns{$col} = \@regexes;
}

open(my $ifh, $tab_file) || die "ERROR: failed to read input file: $!";

while (my $line = <$ifh>) {
    my @cols = split("\t", $line);
    
    my $valid_line = 1;
    
    if ( scalar(@cols) - 1 < $max_column_number ) {
        $valid_line = 0;
        next;
    }
    
    for my $col ( keys %patterns ) {
        
        my $matched_one_pattern = 0;
        
        for my $pattern ( @{$patterns{$col}} ) {
            if ( $cols[$col] =~ /$pattern/ ) {
                $matched_one_pattern++;
            }
        }
        
        if (! $matched_one_pattern) {
            $valid_line = 0;
            last;
        }
    }
    
    print $line if $valid_line;
}

