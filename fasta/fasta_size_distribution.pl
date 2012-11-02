#!/usr/bin/env perl

=head1 NAME

fasta_size_distribution.pl - report the distribution sizes in a collection of fasta files

=head1 SYNOPSIS

USAGE: template.pl 
            --input_file=/path/to/file1.fsa,/path/to/file2.fsa
            --input_list=/path/to/some_file.list 
            --size_interval=50
          [ --complete_orfs_only=0 ]

=head1 OPTIONS

B<--input_file>
    A FASTA input file.  Can be comma-separated.

B<--input_list>
    A text file that contains the paths to all the input files, one per line.  (Can
    be used in place of or in addition to --input_file)

B<--size_interval>
    Will be the range of sizes reported by the script.  If you pass 50, for example, the
    script will show counts of seqs of length between 0-49, 50-99, 100-149, etc.

B<--gc_interval>
    Optional.  Will report size intervals broken down by GC content interval.  (default = none)

B<--complete_orfs_only>
    Optional.  Will only consider sequences with a start/stop, with alternative starts
    allowed.  This option can only be used with nucleotide sequences. (default = 0)

B<--debug> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script reads a collection of FASTA files and generates a simple report on their
size distribution given a user-specified size interval.

=head1  INPUT

The input list is a plain-text file that contains to full paths to each FASTA input
file to be scanned.  Each FASTA file can have one more sequences.  Compressed input
files (via gzip) will be automatically read in their compressed state.

=head1  OUTPUT

The output is a report of the number of sequences file at each interval specified
by the size_interval and gc_interval options.  All non-data lines begin with the # symbol,
allowing parsers to skip blank and commented lines.  Here's an example:

    command:
    
    fasta_size_distribution.pl --input_file=TS20.1_nodupNH.fasta --size_interval=50 --gc_interval=20

    # sequence size distribution

    # GC    size    count
    0-20    50-100  146
    0-20    100-150 51
    0-20    150-200 20
    0-20    200-250 17
    0-20    250-300 26
    0-20    300-350 4
    # 0-20 avg size: 122.0

    20-40   0-50    608
    20-40   50-100  16445
    20-40   100-150 14695
    20-40   150-200 14885
    20-40   200-250 26394
    20-40   250-300 44010
    20-40   300-350 1186
    20-40   350-400 3
    # 20-40 avg size: 202.3

    40-60   0-50    2962
    40-60   50-100  45963
    40-60   100-150 41615
    40-60   150-200 43832
    40-60   200-250 112630
    40-60   250-300 89943
    40-60   300-350 613
    40-60   350-400 1
    # 40-60 avg size: 196.5

    60-80   0-50    426
    60-80   50-100  5340
    60-80   100-150 3919
    60-80   150-200 4313
    60-80   200-250 14288
    60-80   250-300 10635
    60-80   300-350 53
    # 60-80 avg size: 199.6

    80-100  50-100  14
    80-100  100-150 1
    80-100  150-200 2
    # 80-100 avg size: 82.6


    # total seqs analyzed: 495040
    # total seqs skipped: 0
    # total seqs skipped: 98067424
    # smallest sequence length: 24
    # largest sequence length: 362
    # overall average length: 198.1
    # file count: 1

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Statistics::Descriptive;

my %options = ();
my $results = GetOptions (\%options,
                          'input_file=s', 
                          'input_list=s',
                          'size_interval=s',
                          'gc_interval=i',
                          'complete_orfs_only=i',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

my @input_files;

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

my $stat = Statistics::Descriptive::Full->new();

## structure like:
#   $h->{$min_gc_bound}{$min_size_bound} = $count;
my $sizes;

## structure like
#   $h->($min_gc_bound} = { count => N, sum => N }
my $avg_sizes;

my $file_count = 0;
my $seq_skipped_count = 0;
my $smallest = 100000000000000000000000000;
my $largest = 0;

foreach my $file ( @input_files ) {
    chomp $file;
    next if ( $file =~ /^\s*$/ );
    
    ##  catch compressed files
    if (! -e $file && -e "$file.gz") {
        $file .= '.gz';
    }
    
    my $ifh;
    if ( $file =~ /\.gz$/ ) {
        open( $ifh, "<:gzip", $file ) || die "can't read file $file: $!";
    } else {
        open ( $ifh, "<$file" ) || die "can't read file $file: $!";
    }
    
    $file_count++;
    
    ## get the sequences from the file
    my $seq = '';
    my $is_first = 1;
    while (<$ifh>) {
        chomp;
    
        if ( /^\>/ ) {
            if ( $is_first ) {
                $is_first = 0;
            } else {
                check_and_reset(\$seq);
            }
        } else {
            $seq .= $_;
        }
    }
    
    ## purge the last one
    check_and_reset(\$seq);
}

## reporting
print "\n# sequence size distribution\n";
print "\n# GC\tsize\tcount\n";

my $total_count = 0;
my $total_sum = 0;


for my $min_gc_bound ( sort {$a<=>$b} keys %$sizes ) {
    
    for my $min_size_bound ( sort {$a<=>$b} keys %{$$sizes{$min_gc_bound}} ) {

        print "$min_gc_bound-", ($min_gc_bound + $options{gc_interval}), "\t",
              "$min_size_bound-", ($min_size_bound + $options{size_interval}), "\t", 
              "$$sizes{$min_gc_bound}{$min_size_bound}\n";
 
        $total_count += $$sizes{$min_gc_bound}{$min_size_bound};         
    }
    
    print "# $min_gc_bound-",($min_gc_bound + $options{gc_interval})  ," avg size: " . sprintf( "%.1f", $$avg_sizes{$min_gc_bound}{sum} / $$avg_sizes{$min_gc_bound}{count} ) . "\n";
    $total_sum += $$avg_sizes{$min_gc_bound}{sum};
    
    print "\n";
}


print "\n# total seqs analyzed: $total_count\n";
print "# total seqs skipped: $seq_skipped_count\n";
print "# sum sequence length: $total_sum\n";
print "# smallest sequence length: $smallest\n";
print "# largest sequence length: $largest\n";
print "# overall average length: ", sprintf("%.1f", $total_sum / $total_count), "\n";
print "# overall median length: ", sprintf("%.1f", $stat->median()), "\n";
print "# file count: $file_count\n\n";

exit(0);


sub check_and_reset {
    my $seq = shift;
    $$seq =~ s/\s//g;

    ## are we checking for complete ORFs only?
    if ( $options{complete_orfs_only} ) {
        unless ( $$seq =~ /^(ATG|AUG|GTG|GUG|TTG|UUG)/i &&
                 $$seq =~ /(TAA|UAA|TAG|UAG|TGA|UGA)$/i ) {
            $$seq = '';
            $seq_skipped_count++;
            return;
        }
    }
    
    my $seq_len = length($$seq);
    
    # calculate GC
    my $gc_count = ( $$seq =~ tr/GC// );
    my $gc_perc = int( ($gc_count / $seq_len * 100 ) );
    my $min_gc_bound = int( $gc_perc/$options{gc_interval} ) * $options{gc_interval};
    
    ## need to find the size range
    $$sizes{$min_gc_bound}{ int( $seq_len / $options{size_interval} ) * $options{size_interval} }++;

    $$avg_sizes{$min_gc_bound}{count}++;
    $$avg_sizes{$min_gc_bound}{sum} += $seq_len;

    if ( $seq_len < $smallest ) {
        $smallest = $seq_len;
    }
    
    if ( $seq_len > $largest ) {
        $largest = $seq_len;
    }

    $stat->add_data($seq_len);

    ## reset for next sequence
    $$seq = '';
}

sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( size_interval );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ## handle some defaults
    $options{complete_orfs_only} = 0  unless (defined $options{complete_orfs_only});
    $options{gc_interval} = 100  unless (defined $options{gc_interval});
    
    ## parse the input params
    get_input_files();
    
    ## need to make sure files were passed
    if ( scalar @input_files < 1 ) {
        die "must specify input files with any combination of --input_file or --input_list";
    }
}


sub get_input_files {

    if ( $options{input_list} ) {
        open(my $ilh, $options{input_list}) || die "can't read input list: $!";
        
        while ( my $line = <$ilh> ) {
            chomp $line;
            next if $line =~ /^\s*$/;
            
            push @input_files, $line;
        }
    }
    
    if ( $options{input_file} ) {
        my @files = split(',', $options{input_file} );
        push @input_files, @files;
    }
}
