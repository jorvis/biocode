#!/usr/bin/perl

=head1 NAME

hmmer2htab.pl - creates an htab file from hmmscan or hmmpfam raw input depending on the version of hmmer passed using HmmTools.pm

=head1 SYNOPSIS

USAGE: hmmer2htab.pl
    --input_file=/path/to/hmmer.raw
    --output_htab=/path/to/hmmer.htab

=head1 OPTIONS

B<--input_file,-i>
    Raw output from hmmscan or hmmpfam run

B<--output_htab,-o>
    HTAB output file

B<--mldbm_file,-m>
    MLDBM perl data structure (tied hash) containing HMM information.  This was previously
    queried out of egad.hmm2.  See hmmlib_to_mldbm.pl for more information.

B<--prog_ver,-p>
    The version of hmmer [hmmer2 or hmmer3] used to produce the raw output. Default is hmmer3

B<--donotrun>
    skip pretty much any noticeble activity

B<--log,-l>
    Logfile.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

This script is used to convert the output from hmmscan or hmmpfam into htab using $ANNOT_DEVEL/hmm/bin/HmmTools.pm

=head1  INPUT

Output can be in multisequence format (generated from hmmscan or hmmpfam from multifasta input).

Can also be a list of file names.

example hmmscan raw output:

# hmmscan :: search sequence(s) against a profile database
# HMMER 3.0 (March 2010); http://hmmer.org/
# Copyright (C) 2010 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query sequence file:             /usr/local/projects/aengine/sagrawal/output_repository/translate_sequence/2722941276_translate/i1/g1/testApr6.polypeptide.2723277720.1.fsa
# target HMM database:             /export/test_sagrawal/combined_hmm.lib
# output directed to file:         /export/test_sagrawal/polytest1.raw
# prefer accessions over names:    yes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       testApr6.polypeptide.2723277720.1  [L=285]
Scores for complete sequence (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Model      Description
    ------- ------ -----    ------- ------ -----   ---- --  --------   -----------
    7.8e-76  254.9  35.0    8.9e-76  254.7  24.2    1.0  1  PF00950.11 ABC 3 transport family
    7.7e-42  143.7  31.1    8.7e-42  143.5  21.6    1.0  1  TIGR03770  anch_rpt_perm: anchored repeat-type ABC transport
    9.3e-07   28.0  30.3    9.3e-07   28.0  21.0    1.3  1  PF01032.12 FecCD transport family
  ------ inclusion threshold ------
       0.18   11.1   0.8       0.31   10.3   0.5    1.4  1  PF08690.4  GET complex subunit GET2


Domain annotation for each model (and alignments):
>> PF00950.11  ABC 3 transport family
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  254.7  24.2   2.2e-79   8.9e-76       1     257 [.      10     265 ..      10     266 .. 0.99 

Alignments for each domain:
== domain 1    score: 254.7 bits;  conditional E-value: 2.2e-79
                      PF00950.11   1 qyeflqrallasilvslvagilgsfvvlrrqsllgdalsHaalpGvalafllginlllgailvgliaavaiellkr 76  
                                     q+ f+ +all+s++vs+ +++l++f+vl+ ++l+gda+sHa++pG++la++lg++l++ga+++g+++ava ++lk 
testApr6.polypeptide.2723277720.1 10 QFTFMSHALLISLVVSIPCALLSVFLVLKSWALMGDAMSHAVFPGIVLAWILGLPLATGAFVAGVFCAVATGYLKD 85  
                                     689************************************************************************* PP

...
...
...

Internal pipeline statistics summary:
-------------------------------------
Query sequence(s):                         1  (285 residues)
Target model(s):                       16296  (3513130 nodes)
Passed MSV filter:                       920  (0.0564556); expected 325.9 (0.02)
Passed bias filter:                      260  (0.0159548); expected 325.9 (0.02)
Passed Vit filter:                        34  (0.0020864); expected 16.3 (0.001)
Passed Fwd filter:                        13  (0.000797742); expected 0.2 (1e-05)
Initial search space (Z):              16296  [actual number of targets]
Domain search space  (domZ):               4  [number of targets reported over threshold]
# CPU time: 0.52u 3.87s 00:00:04.39 Elapsed: 00:00:07.35
# Mc/sec: 136.22
//

example hmmscan raw output:

Logical Depth LDhmmpfam v1.5.4
Copyright (C) Logical Depth, Inc. All rights reserved.
TIGR and TIGR Affiliates 300 CPU-socket License

hmmpfam
HMMER 2.3-compatible (LDhmmpfam)

HMM file:      /usr/local/db/HMM_LIB/ALL_LIB_bin.HMM
Sequence file: /usr/local/annotation/MOORE/output_repository/translate_sequence/8422_translate_promoted/12/cya1.polypeptide.724320.1.fsa

Query sequence: cya1.polypeptide.724320.1
Accession:      [none]
Description:    [none]

Scores for sequence family classification (score includes all domains):
Model     Description                                   Score    E-value #D
--------  -----------                                   -----    ------- --
TIGR01730 RND_mfp: efflux transporter, RND family, MF   114.8    3.5e-31  1
TIGR00998 8a0101: efflux pump membrane protein          -97.8    1.2e-06  1
TIGR00999 8a0102: Membrane Fusion Protein cluster 2 p    -0.3    4.9e-06  1
TIGR01843 type_I_hlyD: type I secretion membrane fusi  -170.3    0.00012  1

...
...
...

Parsed for domains:
Model     Domain Seq-f Seq-t    HMM-f HMM-t      Score  E-value
--------  ------ ----- -----    ----- -----      -----  -------
TIGR03007  1/1      23   432 ..     1   510 []  -337.6     0.54
TIGR01000  1/1      29   442 ..     1   476 []  -289.8     0.25
TIGR01133  1/1      30   293 ..     1   368 []  -183.6      8.9
TIGR00998  1/1      32   434 ..     1   379 []   -97.8  1.2e-06

...
...
...

Alignments of top-scoring domains:
TIGR03007: domain 1 of 1, from 23 to 432: score -337.6, E = 0.54
                   *->eqllsYlkgiWrr.RwlfvavAwvVmivGwvvvyvlPdrYeAsarVY
                      e+ +   +   +++Rwl+ +v +  +i+ w                 
  cya1.polyp    23    EENRQNTTKNKQFpRWLIPIVILGGGITLWQ---------------- 53   

                   VDTQsvLrPLlkGlAvtPnvdqkirIlsrtLlS.....RpnLekVirmlD
                        + +PL+   + t n     + ++ +LlS+++++R +  +++ +++
  cya1.polyp    54 -----IFSPLVIPTTETNNQTPPPKPVETVLLSsgqgnRQV--RLLGQVE 96   

                   LDvgakspaqlEalitklqknIsIslagrdNLFtISYeDkdPelA.....
                   +  +a+   q  + ++k+  +   s++   + +    +D+d++ A  + +
  cya1.polyp    97 AGAKATLSSQVSGTVEKILVKEGDSITS--GMIVAILDDADGKIAlaeaq 144

=head1 OUTPUT

    Description of the output format (tab-delimited, one line per domain hit)

    col  perl-col   description
    1      [0]      HMM accession
    2      [1]      Date search was run (if available), otherwise date of htab parse
    3      [2]      Length of the HMM (not populated if -s is used)
    4      [3]      Search program
    5      [4]      Database file path
    6      [5]      Sequence accession
    7      [6]      Alignment start position on HMM match - hmm-f
    8      [7]      Alignment end position on HMM match - hmm-t
    9      [8]      Alignment start position on sequence - seq-f
    10     [9]      Alignment end position on sequence - seq-t
    11     [10]     frame (only populated if --frames search is run on nucleotide sequence)
    12     [11]     Domain score
    13     [12]     Total score
    14     [13]     Index of domain hit
    15     [14]
    16     [15]     HMM description (may be truncated by hmmsearch or hmmpfam if -s is used)
    17     [16]     Sequence description (may be truncated by hmmsearch or hmmpfam)
    18     [17]     Total score trusted cutoff (not populated if -s is used)
    19     [18]     Total score noise cutoff (not populated if -s is used)
    20     [19]     Expect value for total hit
    21     [20]     Expect value for domain hit
    22     [21]     Domain score trusted cutoff (egad..hmm2.trusted_cutoff2) (not populated if -s is used)
    23     [22]     Domain score noise cutoff (egad..hmm2.noise_cutoff2) (not populated if -s is used)
    24     [23]     Total score gathering threshold (not populated if -s is used)
    25     [24]     Domain score gathering threshold (not populated if -s is used)

=head1  CONTACT

    Sonia Agrawal
    sagrawal@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../lib";
use HmmTools;
use MLDBM 'DB_File';
use Fcntl qw( O_RDONLY );

####### GLOBALS AND CONSTANTS ###########
my @input_files;                   #Holds input files
my $output_htab;                   #Output htab file
my $output_alignment;              #Output alignment file
my $debug;                         #The debug variable
my $donotrun = 0;
my $alignment_format;              #Holds the format to print output formats
my %alignment_formats =            #Accepted output alignment formats
    ( 'mul' => 1,
      'mfs' => 1,
      'fasta' => 1 );
                         
########################################

my %options = ();
my $results = GetOptions (\%options, 
                          'donotrun=s',
                          'input_file|i=s',
                          'output_htab|o=s',
                          'mldbm_file|m=s',
			              'prog_ver|p=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;

if ($options{'donotrun'}) {
    #$logger->warn("Exiting because --donotrun says to.");
    exit(0);
}

# Check the options.
&check_parameters(\%options);

#Gather information about the hmms.
tie(my %hmm_info, 'MLDBM', $options{mldbm_file}, O_RDONLY ) or die("Unable to tie hash to $options{mldbm_file}");
foreach my $file (@input_files) {

    #Get the output file names;
    my $htab_file = $output_htab;
    if( -d $output_htab ) {
        my $base = $1 if($file =~ m|.*/([^/]+)\.[^/\.]+$| );
        $htab_file = $output_htab.$base.".htab";
    }

    #Generate the htab file.
    system( "rm -f $htab_file" ) if( -e $htab_file );
    my $hmm_data = &generate_htab( $file, $htab_file, \%hmm_info );
}


######################## SUB ROUTINES #######################################
sub generate_htab {
    my ( $file, $outfile, $hmm_db_info ) = @_;
    #If the hmmer output file is in multi sequence format, HmmTools.pm can't handle it.
    #So we will add the functionality here.
    my $tmp_dir = "/tmp/hmmpfam32htab/$options{prog_ver}/$$"; #append hmmer version and process id
    my @tmp_files = &write_tmp_files( $file, $tmp_dir );

    foreach my $tmp_file ( @tmp_files ) {
	my $data;
	if ($options{'prog_ver'} eq 'hmmpfam') {
		$data = &read_hmmer2_output( $tmp_file );
	} elsif ($options{'prog_ver'} eq 'hmmpfam3') {
        	$data = &read_hmmer3_output( $tmp_file );
        }
	my $htab_h;
        open( $htab_h, ">> $outfile") or die("Unable to open $outfile for writing ($!)");
	
	if ((keys %$data) <=0) {
        die("No results to be reported in .htab file. Empty $outfile file created");
		#$logger->logdie("No results to be reported in .htab file. Empty $outfile file created");
	} else {
        	&print_htab( $data, $hmm_db_info, $htab_h );
        }
	close( $htab_h );
    }

    #Remove the tmp directory
    system( "rm -rf $tmp_dir" );

}

sub write_tmp_files {
    my ($file, $outdir) = @_;
    my @files;
    my $header;

    open( RAW, "< $file" ) or die("Unable to open $file ($!)");
    system( "mkdir -p $outdir" );
    my ($oh, $flag);
    while( my $line = <RAW> ) {
        chomp($line);
        if((($line =~ /Query\:\s+(\S*)/) && ($options{'prog_ver'} eq 'hmmpfam3')) || (($line =~ /Query sequence\:\s+(.*)/) && ($options{'prog_ver'} eq 'hmmpfam'))) {
	    my $base = $1;
	    $base =~ s|/||g; # remove slashes from file name
            my $tmp_file = "$outdir/$base.tmp.raw";
            push( @files, $tmp_file );
            close($oh) if($oh);
            open( $oh, "> $tmp_file") or die("Can't open temp file for writing $tmp_file ($!)");
            print $oh "$header$line\n";
            $flag = 1;
        } elsif( $flag ) {
            print $oh $line."\n";
        } else {
            $header .= $line."\n";
        }
    }
    close($oh) if($oh);

    return @files;

}


sub check_parameters {
    my $options = shift;

    &_pod if($options{'help'});
    ## mldbm file must be passed
    if ( ! $options{mldbm_file} ) {
        die("Option mldbm_file is required\n");
        #$logger->logdie("Option mldbm_file is required\n");
    }
    
    if (! $options{'prog_ver'}) {
	$options{'prog_ver'} = 'hmmpfam3';
    }

    if($options{'input_file'}) {
        die("Option input_file ($options{'input_file'}) does not exist\n") unless(-e $options{'input_file'});
        #$logger->logdie("Option input_file ($options{'input_file'}) does not exist\n") unless(-e $options{'input_file'});
        my $infile = $options{'input_file'};
        open( IN, "< $infile") or die("Unable to open $infile ($!)");
        chomp( my $first_line = <IN> );
        if( -e $first_line ) {
            chomp( @input_files = <IN> );
            push( @input_files, $first_line );
        } else {
            push( @input_files, $infile );
        }
        close IN;

    } else {
        die("Option input_file is required\n");
        #$logger->logdie("Option input_file is required\n");
    }

    unless($options{'output_htab'}) {
        die("Option output_htab is required\n");
        #$logger->logdie("Option output_htab is required\n");
    } else {
        $output_htab = $options{'output_htab'};
        if( @input_files > 1 && -f $output_htab ) {
            $output_htab = $1 if($output_htab =~ m|^(.*)/[^/]+|);
            #$logger->warn("Using $output_htab as output directory because a list of files was passed in");
        }
        
    }

    if($options{'debug'}) {
        $debug = $options{'debug'};
    }
    
}

sub _pod {   
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
