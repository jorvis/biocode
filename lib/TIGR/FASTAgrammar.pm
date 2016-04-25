# $Id: FASTAgrammar.pm,v 1.1 2004/04/28 15:03:43 aphillip Exp $

# Copyright @ 2002 - 2010 The Institute for Genomic Research (TIGR).
# All rights reserved.
# 
# This software is provided "AS IS".  TIGR makes no warranties, express or
# implied, including no representation or warranty with respect to the
# performance of the software and derivatives or their safety,
# effectiveness, or commercial viability.  TIGR does not warrant the
# merchantability or fitness of the software and derivatives for any
# particular purpose, or that they may be exploited without infringing the
# copyrights, patent rights or property rights of others.
# 
# This software program may not be sold, leased, transferred, exported or
# otherwise disclaimed to anyone, in whole or in part, without the prior
# written consent of TIGR.

package TIGR::FASTAgrammar;
{

=head1 NAME

FASTAgrammar - module for validating FASTA format records

=head1 SYNOPSIS

  use TIGR::FASTAgrammar ':public';

  $header = FASTA header here...
  $data = FASTA data here...
  $return_value = isValidFASTARecord($header, $data);
  ...

=head1 DESCRIPTION

This module provides functions for verifying compliance with TIGR's FASTA
file and record definitions.

=cut

   BEGIN {
      require 5.006_00;
   }

   use strict;
   require Exporter;

   ## internal variables and identifiers

   our @ISA = qw(Exporter);
   our $REVISION = (qw$Revision: 1.1 $)[-1];
   our $VERSION = '1.2';
   our $VERSION_STRING = "$VERSION (Build $REVISION)";
   our @DEPEND = ();


   ## Export methods

   our %EXPORT_TAGS = ( 'public'  => [ qw( isValidFASTArecord
                                           isValidFASTAheader
                                           isValidFASTAdata
                                           isValidFASTAlineLength
                                           setValidFASTAlineLength ) ] ,
                        'private' => [ qw( _headerToIdentifier
                                           _isNucleotideData
                                           _isPeptideData ) ] );

   our @EXPORT_OK = ( @{ $EXPORT_TAGS{'public'} },
                      @{ $EXPORT_TAGS{'private'} } );


   ## data structures
   
   # IUPAC extended codes acceptable for sequence data
   # see WU-BLAST format on http://tigrblast.tigr.org/html/fasta.html
   our $NA_IUPAC_CODES = 'ATUGCMRWSYKVHDBN\.\-';
   our $AA_IUPAC_CODES = 'ARNDBCQEZGHILKMFPSTUWXYV\*\-';

   # FASTA parameters
   our $FASTA_SEPARATOR = '^>';
   our $UNBOUND_FASTA_SEPARATOR = '>';

   # note: BLAST can accept 80 bases per line; this just emits a warning
   our $OUTPUT_LINE_LENGTH = 60;
   my $RECORD_LINE_LENGTH = 0;
 
   # the TIGR FASTA header parse
   our $FASTA_HEADER_PARSE = 
            '^>([[:graph:]]+)(?: +[ \ca[:graph:]]*$| *$)';
   ## prototypes
 
   sub isValidFASTArecord(@);
   sub isValidFASTAheader($);
   sub isValidFASTAdata($);
   sub isValidFASTAlineLength($);
   sub setValidFASTAlineLength($);
   sub _isNucleotideData($);
   sub _isPeptideData($);
   sub _headerToIdentifier($);

   ## implementation

=over

=item $result = isValidFASTArecord(@record_defn);

This method determines if a FASTA record, C<@record_defn>, fits the TIGR
definition for a FASTA record.  C<@record_defn> is an array of lines over
which the record is defined.  The first line should be the FASTA header, and
subsequent lines the data definition.  This method checks line width, 
character composition, and header format.  If the record parses correctly,
this method returns 1.  Otherwise, this method returns 0.

=cut


   sub isValidFASTArecord(@) {
     
      my $header = shift;
      my @data_lines = @_;
      my $valid_flag = 0;
      my $first_line_flag = 0;
      my $first_len = 0;
     
      # check conformance of header
      $valid_flag = isValidFASTAheader($header);

      # check conformance of data
      if ( ( $valid_flag != 0 ) &&
           ( scalar(@data_lines) > 0 ) ) {
         my $data_scalar = join "", @data_lines;
         $data_scalar =~ s/\n//g; # extract new lines from scalar data
         $valid_flag = isValidFASTAdata($data_scalar);
      }
      
      # check conformance of line length
      while ( ( $valid_flag != 0 ) &&
              ( defined ( my $line = shift @data_lines ) ) ) {
         chomp $line;
         
	 if($first_line_flag == 0) {
	    $first_len = setValidFASTAlineLength($line);
            $first_line_flag = 1;
	 }
        
         if(defined($first_len)) {
            my $line_len_flag = isValidFASTAlineLength($line);
            if ( ( $line_len_flag > 0 ) ||
               ( ( $line_len_flag < 0 ) &&
               ( $#data_lines != -1 ) ) ) {
               $valid_flag = 0;
            }
         }
      }

      return $valid_flag;
   }

=item $result = isValidFASTAheader($header);

This method determines if the FASTA header description, C<$header>, is
a valid description.  It checks for a leading carot symbol and trailing non
white space characters.  Any number of white space characters
may be interleaved throughout the text portion of the header, with the
exception that there may be no space between the carot and the first word.
If the header is valid, this method returns 1.  Otherwise, this method
returns 0.

=cut


   sub isValidFASTAheader($) {
      my $header = shift;
      my $return_val = 0;
      my $identifier = undef;
      
      if ( ( defined ($header) ) &&
           ( ($identifier) = $header =~ /$FASTA_HEADER_PARSE/ ) ) {
         
         if((defined $identifier) && ($identifier !~ /\//)) {
            $return_val = 1;
         }
         else {
            $return_val = 0;
	 }
      }
      else {
        $return_val = 0;
      }
      return $return_val;
   }

=item $result = isValidFASTAdata($data_def);

This method takes the scalar data definition of a FASTA record, C<$data_def>.
It tests the data and returns 1 if the data conforms to nucleotide data or if 
it conforms to peptide data or both. If the data is not recognizable or is 
undefined, it returns 0.

=cut

   
   sub isValidFASTAdata($) {
      my $data_definition = shift;
      my $return_val = undef;
      if(($data_definition =~ /^[$NA_IUPAC_CODES]+$/i) || 
            ($data_definition =~ /^[$AA_IUPAC_CODES]+$/i)) {
	  $return_val = 1;
      }
      else {
         $return_val = 0;
      }
   }

=item $result = isValidFASTAlineLength($line);

This method returns -1 if the data line, C<$line> is less than
the TIGR definition requirement for line length, 0 if the data
line meets the TIGR definition requirement for line length, and
1 if the data line is greater than the TIGR definition requirement
for line length.

=cut


   sub isValidFASTAlineLength($) {
      my $line = shift;
      my $line_len = undef;
      my $return_val = undef;

      if ( defined ($line) ) {
         chomp $line;
         $line_len = length($line);
         if ( $line_len > $RECORD_LINE_LENGTH ) {
            $return_val = 1;
         }
         elsif ( $line_len < $RECORD_LINE_LENGTH ) {
            $return_val = -1;
         }
         else {
            $return_val = 0;
         }
      }
   }

=item $result = setValidFASTAlineLength($);

This method takes in the first data line in the data portion of a FASTA record.
The function returns the length of this line if it is positive. This length 
determines the line length for all the data lines following this first line.  
The function returns undefined if unsuccessful.

=cut


   sub setValidFASTAlineLength($) {
      my $line = shift;
      my $line_len = undef;
      my $ret_len = undef;

      if(defined ($line)) {
         chomp $line;
         $line_len = length($line);

	 if($line_len > 0) {
	    $ret_len = $line_len ;
         }
      }
      $RECORD_LINE_LENGTH = $ret_len;
      return $ret_len;
   }
	      
# $result = _isNucleotideData($data_def);

#This method takes the scalar data definition of a FASTA record, C<$data_def>.
#It tests it for conformance to a  nucleotide data type.  If the data are 
#nucleotide IUPAC characters, this method returns 1.  If not, this method 
#returns 0.  This method returns 0 if C<$data_def> is undefined.


   sub _isNucleotideData($) {
      my $data_def = shift;
      my $return_val = undef;

      if ( ( defined ( $data_def ) ) &&
           ( $data_def =~ /^[$NA_IUPAC_CODES]+$/i ) ) {
         $return_val = 1;
      }
      else {
         $return_val = 0;
      }

      return $return_val;
   }


# $result = _isPeptideData($data_def);

#This method takes the scalar data definition of a FASTA record, C<$data_def>.
#It tests it for conformance to a peptide data type.  If the data are
#peptide IUPAC characters, this method returns 1.  If not, this method returns
#zero.  This method returns undefined if C<$data_def> is undefined.


   sub _isPeptideData($) {
      my $data_def = shift;
      my $return_val = undef;

      if ( ( defined ( $data_def ) ) &&
           ( $data_def =~ /^[$AA_IUPAC_CODES]+$/i ) ) {
         $return_val = 1;
      }
      else {
         $return_val = 0;
      }

      return $return_val;
   }

# $identifier = _headerToIdentifier($header);

#This function takes a FASTA header as a parameter, returning a parsed
#identifier for the record.  If the supplied header is invalid or undefined,
#this method returns undefined.

   
   sub _headerToIdentifier($) {
      my $header = shift;
      my $identifier = undef;
      my $return_val = undef;
    
      if ( ( defined ($header) ) &&
           ( ($identifier) = $header =~ /$FASTA_HEADER_PARSE/ ) &&
           ( defined ($identifier) ) ) {
         if($identifier =~ /\//) {
            $return_val = undef;
         }
         else {
            $return_val = $identifier;
         }
      }
      else {
         $return_val = undef;
      }
      
      return $return_val;
   }


=back

=head1 USAGE

This module is not intended for developer use.  Instead, use the front
end modules C<TIGR::FASTAreader> and C<TIGR::FASTArecord>.

=cut

}

1;
