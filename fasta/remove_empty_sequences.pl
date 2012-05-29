#!/usr/bin/env perl

=head1 NAME

remove_empty_sequences.pl - Reads a multi-FASTA file and prints entries with empty 
sequences excluded.

=head1 SYNOPSIS

USAGE: remove_empty_sequences.pl 
            --input_file=/path/to/some_file.fsa

=head1 OPTIONS

B<--input_file,-i>
    A multi-FASTA file

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Many programs fail if passed multi-FASTA files that contain empty sequences.  This
script can be used to filter them out.

=head1  INPUT

A multi-FASTA file.

=head1  OUTPUT

The same sequences with entries with empty residues excluded.  Sequences are written
out to be 60bp per line.  Header lines of entries with no sequences are printed on
STDERR.

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
                          'input_file|i=s',
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

my $ifh = get_conditional_read_fh( $options{input_file} );

my $current_header;
my @current_seq_lines = ();

while ( my $line = <$ifh> ) {
    chomp $line;

    if ( $line =~ /^\>(.+)$/ ) {
        my $next_header = $1;
        
        if ( $current_header ) {
            process_seq( $current_header, \@current_seq_lines );
        }
        
        $current_header = $next_header;
        @current_seq_lines = ();
    } else {
        push @current_seq_lines, $line;
    }
}

## trailing one
process_seq( $current_header, \@current_seq_lines );

exit(0);

sub process_seq {
    my ($header, $lines) = @_;
    
    my $seq = join('', @$lines);
       $seq =~ s|\s||g;
    
    if ( length $seq ) {
        print ">$header\n";
        while ( $seq =~ /(.{1,60})/g ) {
            print "$1\n";
        }
    } else {
        print STDERR "EXCLUDING: $header\n";
    }
}


sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
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


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( input_file );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
}









