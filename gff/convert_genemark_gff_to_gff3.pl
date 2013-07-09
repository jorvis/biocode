#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

    convert_genemark_gff_to_gff3.pl  - Format the annotation file according to gff3 format. 
                                       [Works for output from Genemark tools.]


=head1 SYNOPSIS

    convert_genemark_gff_to_gff3.pl   --a1 <annotation> [--o outdir]
                                      [--o outdir] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

    Assumptions : All features should be CDS. (Genemark)
                  Input file should be sorted by gene ids in column 9.


=head1 OPTIONS
    
    --a1 <annotation>      = Annotation file.

    --o <output dir>       = /path/to/output directory. Optional.[PWD]

    --v                    = generate runtime messages. Optional

=head1 DESCRIPTION

Format the annotation file from Genemark gene prediction tool according to gff3 format.

=head1 AUTHOR

 Priti Kumari
 Bioinformatics Software Engineer 
 Institute for Genome Sciences
 University of Maryland
 Baltimore, Maryland 21201

=cut

##############################################################################

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Spec;

##############################################################################
### Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

#use constant BIN_DIR => '/usr/local/bin';


use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM."\n";

GetOptions( \%hCmdLineOption,
            'outdir|o=s', 'annotation1|a1=s',
	    'verbose|v', 'debug',
            'help', 
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

## Define variables
my $sOutDir;
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;
my ($sRefID, $sSource, $sFeature, $sStart, $sEnd, $sScore, $sStrand, $sFrame, $sAttributes) ;
my ($count, $prefix, $i, $cmd, $t_count);
my ($exon_ID, $cds_ID, $sID, $stranscript);
my ($chr, $gene, $exon, ,$score, $frame, $source, $strand, $start, $end, $gene_start, $gene_stop, $trans);
my ($fh, $fout, $ftemp);
my (%gff, %ref, %att);
my (@exon , @cds, @atribs, @transcript, @temp);



################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'infile'} ...\n" : ();

$sOutDir = File::Spec->curdir();
if (defined $hCmdLineOption{'outdir'}) {
    $sOutDir = $hCmdLineOption{'outdir'};

    if (! -e $sOutDir) {
        mkdir($hCmdLineOption{'outdir'}) ||
            die "ERROR! Cannot create output directory\n";
    }
    elsif (! -d $hCmdLineOption{'outdir'}) {
            die "ERROR! $hCmdLineOption{'outdir'} is not a directory\n";
    }
}

$sOutDir = File::Spec->canonpath($sOutDir);


($_, $_, $prefix) = File::Spec->splitpath($hCmdLineOption{'annotation1'}) ;
$prefix =~ s/\.gtf//;
$prefix =~ s/\.gff//;

print "Reading annotation file ...\n\n" ;
open($fh, "<$hCmdLineOption{'annotation1'}") or die "Cannot open the annotation file.";
open( $ftemp , ">$sOutDir/temp.gff3") or die "Error.. Cannot open file for writing.";
open( $fout , ">$sOutDir/$prefix\.gff3") or die "Error.. Cannot open file for writing.";
print $fout "##gff-version 3\n";

$count = 1;
$t_count = 1;

while(<$fh>) {
    next if ($_ =~ /^\s+$/);
    next if ($_ =~ m/^\#/);
    chomp($_);
    ($sRefID, $sSource, $sFeature, $sStart, $sEnd, $sScore, $sStrand, $sFrame, $sAttributes) = split(/\t/, $_);
    next if ($sStrand !~ m/[-+]/);
 
    next if ($sFeature ne "CDS");

    ##Determining feature id.. 
    @atribs = split(/;/,$sAttributes);
    foreach (@atribs) {
        if ($_=~/\=/) {
            $att{(split (/\=/,$_)) [0]} = (split (/\=/,$_))[1];
        }
        elsif ($_=~/\s+/){
            $_=~s/\s+//g;
            @temp=split (/\"/,$_);
            $att{$temp[0]} = $temp[1];
        }
        else {die "Unrecognized annotation file";} 
    }

    if (exists $att{'gene_id'}) {
        $sID = $att{'gene_id'};
    }

    if (exists $att{'transcript_id'}) {
        $stranscript = $att{'transcript_id'};
    }

    if (exists $ref{$sRefID}) {$ref{$sRefID}[1] = $sEnd;}
    else {$ref{$sRefID} = [$sStart, $sEnd];}


    ##If gene exists add CDS to the hash
    if (exists $gff{$sRefID}{$sID}) {
	$gff{$sRefID}{$sID}{$stranscript}{$sStart}{$sEnd} = [$sSource, $sFeature, $sScore, $sStrand, $sFrame];
    }
    ##If new gene, print the last gene in the hash...
    else {    

	foreach $chr (keys %gff) {
	    foreach $gene (keys %{$gff{$chr}}) {
		foreach $trans (keys %{$gff{$chr}{$gene}}) {
		    foreach $start (sort {$a<=>$b} keys %{$gff{$chr}{$gene}{$trans}}) {
			foreach $end (keys %{$gff{$chr}{$gene}{$trans}{$start}}) {
			    $exon_ID = $gene.".exon.$count";
			    $source = $gff{$chr}{$gene}{$trans}{$start}{$end}[0] ;
			    if ($count == 1 ) {$gene_start = $start ;}
			    $gene_stop = $end ;
			    $score = $gff{$chr}{$gene}{$trans}{$start}{$end}[2] ;
			    $strand = $gff{$chr}{$gene}{$trans}{$start}{$end}[3];
			    $frame = $gff{$chr}{$gene}{$trans}{$start}{$end}[4];
			    push @exon, "$chr\t$source\texon\t$start\t$end\t$score\t$strand\t.\tID=$exon_ID;Parent=$gene.mRNA.$trans\n";
			    push @cds,  "$chr\t$source\tCDS\t$start\t$end\t$score\t$strand\t$frame\tID=$gene.cds.$t_count;Parent=$gene.mRNA.$trans\n";
			    $count ++;
			}
		    }
		    push (@transcript, "$chr\t$source\tmRNA\t$gene_start\t$gene_stop\t.\t$strand\t.\tID=$gene.mRNA.$trans;Parent=$gene\n");
		    $t_count ++;
		}

		print $ftemp  "$chr\t$source\tgene\t$gene_start\t$gene_stop\t.\t$strand\t.\tID=$gene\n";

		for ($i = 0; $i<scalar @transcript ;$i++) {
		    print $ftemp $transcript[$i] ;
		}

		for ($i = 0; $i<scalar @exon ;$i++) {
		   
		    print  $ftemp $exon[$i] ;
		    print  $ftemp $cds[$i] ;

		}	       
		    
	    }
	}

    	%gff = ();
	@transcript = ();
	@exon = ();
	@cds = ();
	$count = 1;
	$t_count = 1;}
	$gff{$sRefID}{$sID}{$stranscript}{$sStart}{$sEnd} = [$sSource, $sFeature, $sScore, $sStrand, $sFrame];
       
}	

##To print last gene in the hash...		
foreach $chr (keys %gff) {
    foreach $gene (keys %{$gff{$sRefID}}) {
	foreach $trans (keys %{$gff{$sRefID}{$gene}}) {
	    foreach $start (sort {$a<=>$b} keys %{$gff{$chr}{$gene}{$trans}}) {
		foreach $end (keys %{$gff{$chr}{$gene}{$trans}{$start}}) {
		    $exon_ID = $gene.".exon.$count";
		    $source = $gff{$chr}{$gene}{$trans}{$start}{$end}[0] ;
		    if ($count == 1 ) {$gene_start = $start ;}
		    $gene_stop = $end ;
		    $score = $gff{$chr}{$gene}{$trans}{$start}{$end}[2] ;
		    $strand = $gff{$chr}{$gene}{$trans}{$start}{$end}[3];
		    $frame = $gff{$chr}{$gene}{$trans}{$start}{$end}[4];
		    push @exon, "$chr\t$source\texon\t$start\t$end\t$score\t$strand\t.\tID=$exon_ID;Parent=$gene.mRNA.$trans\n";
		    push @cds,  "$chr\t$source\tCDS\t$start\t$end\t$score\t$strand\t$frame\tID=$gene.cds.$t_count;Parent=$gene.mRNA.$trans\n";
		    $count++;
		}
	    }

	    push (@transcript, "$chr\t$source\tmRNA\t$gene_start\t$gene_stop\t.\t$strand\t.\tID=$gene.mRNA.$trans;Parent=$gene\n");
	    $t_count++;
	}

	print $ftemp "$chr\t$source\tgene\t$gene_start\t$gene_stop\t.\t$strand\t.\tID=$gene\n";

	for ($i = 0; $i<scalar @transcript ;$i++) {
	    print $ftemp $transcript[$i] ;
	}

	for ($i = 0; $i<scalar @exon ;$i++) {
	    
	    print  $ftemp $exon[$i] ;
	    print  $ftemp $cds[$i] ;
	    
	}
    }
}


close $fh;
close $ftemp;

foreach (keys %ref) {
    print $fout "##sequence region\t$_\t$ref{$_}[0]\t$ref{$_}[1]\n";
}
close $fout;

$cmd = "cat $sOutDir/temp.gff3 >>$sOutDir/$prefix.gff3";
exec_command($cmd);
$cmd = "rm $sOutDir/temp.gff3";
exec_command($cmd);

	    



################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input files provided
    if (! (defined $phOptions->{'annotation1'}) ) {
	pod2usage( -msg => $sHelpHeader, -exitval => 1);
    }

}
    
sub exec_command {
        my $sCmd = shift;

        if ((!(defined $sCmd)) || ($sCmd eq "")) {
                die "\nSubroutine::exec_command : ERROR! Incorrect command!\n";
        }

        my $nExitCode;

        #print STDERR "\n$sCmd\n";
        $nExitCode = system("$sCmd");
        if ($nExitCode != 0) {
                die "\tERROR! Command Failed!\n\t$!\n";
        }
        print STDERR "\n";

        return;
}  

################################################################################
