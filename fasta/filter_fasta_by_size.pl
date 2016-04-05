#!/usr/bin/env perl

=head1 NAME

filter_fasta_by_size.pl - filter the sequences within a FASTA file or list of FASTA files.

=head1 SYNOPSIS

USAGE: filter_fasta_by_size.pl 
            --input_list=/path/to/some_file.list 
            --min_size_cutoff=50
            --output_directory=/some/dir
          [ --input_file=/path/to/some_file.fsa
            --max_size_cutoff=100
            --output_file=/path/to/some_file.fsa 
            --limit_count=1000
          ]

=head1 OPTIONS

B<--input_list,-i>
    A text file that contains the paths to all the input files, one per line.

B<--input_file,-f>
    Optional. A single input file (can contain one or more sequences). Can be used 
    instead of the input_list option.

B<--min_size_cutoff,-s>
    Sequences smaller than this will not be reported to the output file(s).

B<--max_size_cutoff,-x>
    Optional. Sequences greater than this will not be reported to the output file(s).

B<--output_directory,-o>
    Directory used to write output files.  Files will contain the same name as the
    source files.

B<--output_file,-u>
    Optional.  If used instead of the output_directory option all output will be
    written to this single file.

B<--limit_count,-l>
    Optional.  Limit the result file to the first N sequences within the size cutoff.
    (Default = 0 [unlimited])

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script reads a single input file for a list and filters the FASTA sequences
within based on some size cutoff.  The format of the sequences are not changed
in any way (such as number of bases per line.)

=head1  INPUT

The input list is a plain-text file that contains to full paths to each FASTA input
file to be scanned.  Each FASTA file can have one more sequences.  Compressed input
files (via gzip) will be automatically read in their compressed state.

=head1  OUTPUT

The output will consist of the same FASTA sequences, unmodified but with those not 
meeting the size cutoff omitted.

=head1  CONTACT

    Joshua Orvis (and others)
    jorvis@gmail.com

=cut

use strict;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_list|i=s',
                          'input_file|f=s',
                          'min_size_cutoff|s=i',
                          'max_size_cutoff|x=i',
                          'output_directory|o=s',
                          'output_file|u=s',
                          'limit_count|m=i',
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

my @input_seqs = get_input_seqs();

## the output file handle
my $ofh;
if ( defined $options{output_file} ) {
    $ofh = get_conditional_write_fh( $options{output_file} );
}

my $c_seqs_exported = 0;

foreach my $input_file_path ( @input_seqs ) {
    _log("processing $input_file_path");
    
    ## handle the output file for this sequence
    if (! defined $ofh) {
        ## if we get here we must be writing one output file per input file
        ## conserve the file's basename
        my $basename = basename( $input_file_path );
        
        $ofh = get_conditional_write_fh("$options{output_directory}/$basename");
    }
    
    my $header = '';
    my $seq = '';
    my $seq_length = 0;
    my $ifh = get_conditional_read_fh($input_file_path);
    
    while ( <$ifh> ) {
        if ( /\>/ ) {
            ## if we hit a new header (not the first one) print the current sequence
            if ($header) {
                ## check that sequence meets both min and max criteria (as defined)
                my $has_valid_seq_length = 0;

                ## check min requirement
                if ( $seq_length >= $options{min_size_cutoff} ) {
                    $has_valid_seq_length = 1;
                } 

                 ## check max requirement
                if ( $has_valid_seq_length && defined($options{max_size_cutoff}) ) {
                    if ( $seq_length > $options{max_size_cutoff} ) {
                        $has_valid_seq_length = 0;
                    }
                }
        
                ## print output or log statement
                if ($has_valid_seq_length == 1) {
                    if ( $options{limit_count} && $c_seqs_exported >= $options{limit_count} ) {
                        _log("skipping the follow sequence because the limit_count has been reached:\n$header$seq");
                    } else {
                        print $ofh "$header$seq";
                        $c_seqs_exported++;
                    }
                } else {
                    _log("skipping the following sequence because its length is not within " .
                         "range (${seq_length}bp):\n$header$seq");
                }
            }
        
            ## reset for next round
            $header = $_;
            $seq_length = 0;
            $seq = '';
        } else {
            $seq .= $_;
            s/\s+//g;
            $seq_length += length $_;
        }
    }
    
    ## purge the last sequence of the file
    if ($header) {
        ## check that sequence meets both min and max criteria (as defined)
        my $has_valid_seq_length = 0;
    
        ## check min requirement
        if ( $seq_length >= $options{min_size_cutoff} ) {
            $has_valid_seq_length = 1;
        } 
    
        ## check max requirement
        if ( $has_valid_seq_length && defined($options{max_size_cutoff}) ){
            if ( $seq_length > $options{max_size_cutoff} ) {
                $has_valid_seq_length = 0;
            }
        }

        ## print output or log statement
        if ($has_valid_seq_length == 1) {
            if ( $options{limit_count} && $c_seqs_exported >= $options{limit_count} ) {
                _log("skipping the follow sequence because the limit_count has been reached:\n$header$seq");
            } else {
                print $ofh "$header$seq";
                $c_seqs_exported++;
            }
        } else {
            _log("skipping the following sequence because its length is not within " .
                 "range (${seq_length}bp):\n$header$seq");
        }
    }
    
    ## if the user specified an output directory we need to close
    ##  the current output file handle.  Another will be opened in
    ##  the next iteration.
    if ( defined $options{output_directory} ) {
        close $ofh;
        undef $ofh;
    }
}


exit(0);


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

sub get_input_seqs {
    my @seqs;
    
    ## check for a single file
    if ( defined $options{input_file} ) {
        ## make sure it exists
        if ( -e $options{input_file} ) {
            push @seqs, $options{input_file};
        } else {
            die "the following input file doesn't exist: $options{input_file}";
        }
    }
    
    ## check for a list
    if ( defined $options{input_list} ) {
        if ( -e $options{input_list} ) {
            
            open(my $list_fh, "<$options{input_list}") || die "can't open input list: $!";
            
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
            die "the following input list doesn't exist: $options{input_list}";
        }
    }
    
    return @seqs;
}

sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( min_size_cutoff );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## either input_list or input_file must be defined (both is OK)
    if (! defined $$options{input_list} && ! defined $$options{input_file} ) {
        die "either --input_list or --input_file must be defined";
    }
    
    ## either output_file or output_directory must be defined (but not both)
    if (! defined $$options{output_file} && ! defined $$options{output_directory} ) {
        die "either --output_file or --output_directory must be defined";
    }    
    
    if (defined $$options{output_file} && defined $$options{output_directory} ) {
        die "either --output_file or --output_directory must be defined, but not both";
    }  
    
    
    if (defined $$options{max_size_cutoff} && ($$options{max_size_cutoff} < $$options{max_size_cutoff}) ) {
        die "max_size_cutoff must be greater than min_size_cutoff";
    }  
    
    ## handle some defaults
    $options{limit_count} = 0  unless (defined $options{limit_count});
}









