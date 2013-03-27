#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

    convert_snap_Fgenesh_gff_to_gff3.pl  - Format the annotation file according to gff3 format. 
                                           [Works for output from SNAP and FGENESH tools.]


=head1 SYNOPSIS

    convert_snap_Fgenesh_gff_to_gff3.pl  --a1 <annotation> [--o outdir]
                                         [--o outdir] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

    Assumptions : All features should be exons. (SNAP, FGENESH)
                  Input file should be sorted gene ids in column 9.


=head1 OPTIONS
    
    --a1 <annotation>      = Annotation file.

    --o <output dir>       = /path/to/output directory. Optional.[PWD]

    --v                    = generate runtime messages. Optional

=head1 DESCRIPTION

Format the annotation file from SNAP and FGENESH gene prediction tool according to gff3 format.

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
my ($count, $prefix, $i, $cmd);
my ($exon_ID, $cds_ID, $mrna_ID);
my ($chr, $gene, $exon, ,$score, $frame, $source, $strand, $start, $end, $gene_start, $gene_stop);
my ($fh, $fout, $ftemp);
my (%gff, %ref);
my (@exon , @cds);



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
open($fh, "<$hCmdLineOption{'annotation1'}") or die "Cannot open the list file.";
open( $ftemp , ">$sOutDir/temp.gff3") or die "Error.. Cannot open file for writing.";
open( $fout , ">$sOutDir/$prefix\.gff3") or die "Error.. Cannot open file for writing.";
print $fout "##gff-version 3\n";

$count = 1;
while(<$fh>) {
    
    next if ($_ =~ /^\s+$/);
    next if ($_ =~ m/^\#/);
    chomp($_);
    ($sRefID, $sSource, $sFeature, $sStart, $sEnd, $sScore, $sStrand, $sFrame, $sAttributes) = split(/\t/, $_);
    next if ($sStrand !~ m/[-+]/);
    if ( ! ($sFeature =~ /exon/i) ) {die "Feature in column 3 is not exon.";}

    if (exists $ref{$sRefID}) {$ref{$sRefID}[1] = $sEnd;}
    else {$ref{$sRefID} = [$sStart, $sEnd];}

    if (exists $gff{$sRefID}{$sAttributes}) {
	$gff{$sRefID}{$sAttributes}{$sStart}{$sEnd} = [$sSource, $sFeature, $sScore, $sStrand, $sFrame];
    }
    else {

	foreach $chr (keys %gff) {
	    foreach $gene (keys %{$gff{$chr}}) {
		foreach $start (sort {$a<=>$b} keys %{$gff{$chr}{$gene}}) {
		    foreach $end (keys %{$gff{$chr}{$gene}{$start}}) {
			$exon_ID = $gene.".exon.$count";
			
			$source = $gff{$chr}{$gene}{$start}{$end}[0] ;
			if ($count == 1 ) {$gene_start = $start ;}
			$gene_stop = $end ;
			$score = $gff{$chr}{$gene}{$start}{$end}[2] ;
			$strand = $gff{$chr}{$gene}{$start}{$end}[3];
			$frame = $gff{$chr}{$gene}{$start}{$end}[4];
			
			push @exon, "$chr\t$source\texon\t$start\t$end\t$score\t$strand\t$frame\tID=$exon_ID;Parent=$gene.mRNA.1\n";
			push @cds,  "$chr\t$source\tCDS\t$start\t$end\t$score\t$strand\t0\tID=$gene.cds.1;Parent=$gene.mRNA.1\n";
			$count ++;
			
		    }
		}
		    
		$mrna_ID = $gene.".mRNA.1";
		print $ftemp "$chr\t$source\tgene\t$gene_start\t$gene_stop\t.\t$strand\t$frame\tID=$gene\n";
		print $ftemp "$chr\t$source\tmRNA\t$gene_start\t$gene_stop\t.\t$strand\t$frame\tID=$mrna_ID;Parent=$gene\n";
		
		for ($i = 0; $i<scalar @exon ;$i++) {
		    
		    print $ftemp $exon[$i] ;
		    print $ftemp $cds[$i] ;

		}
	    }
	}
	%gff = ();
	@exon =();
	@cds = ();
	$count = 1;
	    
    }
    $gff{$sRefID}{$sAttributes}{$sStart}{$sEnd} = [$sSource, $sFeature, $sScore, $sStrand, $sFrame];
    
}

foreach $chr (keys %gff) {
    foreach $gene (keys %{$gff{$chr}}) {
	foreach $start (sort {$a<=>$b} keys %{$gff{$chr}{$gene}}) {
	    foreach $end (keys %{$gff{$chr}{$gene}{$start}}) {
		$exon_ID = $gene.".exon.$count";
		
		$sSource = $gff{$chr}{$gene}{$start}{$end}[0] ;
		if ($count == 1 ) {$gene_start = $start ;}
		$gene_stop = $end ;
		$score = $gff{$chr}{$gene}{$start}{$end}[2] ;
		$strand = $gff{$chr}{$gene}{$start}{$end}[3];
		$frame = $gff{$chr}{$gene}{$start}{$end}[4];
		    
		push @exon, "$chr\t$source\texon\t$start\t$end\t$score\t$strand\t$frame\tID=$exon_ID;Parent=$gene.mRNA.1\n";
		push @cds,  "$chr\t$source\tCDS\t$start\t$end\t$score\t$strand\t0\tID=$gene.cds.1;Parent=$gene.mRNA.1\n";
		$count ++;
		
	    }
	}
	
	$mrna_ID = $gene.".mRNA.1";
	print $ftemp "$chr\t$source\tgene\t$gene_start\t$gene_stop\t.\t$strand\t$frame\tID=$gene\n";
	print $ftemp "$chr\t$source\tmRNA\t$gene_start\t$gene_stop\t.\t$strand\t$frame\tID=$mrna_ID;Parent=$gene\n";
	    
	for ($i = 0; $i<scalar @exon ;$i++) {
	     
	    print $ftemp $exon[$i] ;
	    print $ftemp $cds[$i] ;

	}
    }
}   

close $fh;
close $ftemp;

foreach (keys %ref) {
    print $fout "##sequence region\t$_\t$ref{$_}[0]\t$ref{$_}[1]\n";
}
close $ftemp;

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
