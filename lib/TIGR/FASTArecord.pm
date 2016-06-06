# $Id: FASTArecord.pm,v 1.1 2004/04/28 15:03:43 aphillip Exp $

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

package TIGR::FASTArecord;
{

=head1 NAME

TIGR::FASTArecord - TIGR::FASTArecord class describing FASTA records

=head1 SYNOPSIS

  use TIGR::FASTArecord;
  my $obj_instance = new TIGR::FASTArecord ($record_header,
                                            $record_data);

=head1 DESCRIPTION

This module provides an object definition for a FASTA record.  It verifies
data entry on creation, and returns information queried on the record.

=cut

   BEGIN {
      require 5.006_00;
   }

   use strict;
   use TIGR::FASTAgrammar ':public';
   use TIGR::FASTAgrammar ':private';

   ## external variables

   my $OUTPUT_LINE_LENGTH_REF = \$TIGR::FASTAgrammar::OUTPUT_LINE_LENGTH;
   my $UNBOUND_FASTA_SEPARATOR = $TIGR::FASTAgrammar::UNBOUND_FASTA_SEPARATOR;
  
   ## internal variables and identifiers

   our $REVISION = (qw$Revision: 1.1 $)[-1];
   our $VERSION = '1.2';
   our $VERSION_STRING = "$VERSION (Build $REVISION)";
   our @DEPEND = ();

   ## prototypes

   sub new($$);
   sub equals($);
   sub getHeader();
   sub getIdentifier();
   sub getData();
   sub size();
   sub toString();
   sub reverseComplement($);
   sub reverseComplementData();
   sub subSequence($;$);
   sub unGap();

   ## implementation

=over

=item $obj_instance = new TIGR::FASTArecord ($header, $data_rec);

This method returns a new instance of a TIGR::FASTArecord object.  It takes
a record header, C<$header>, and the record data, C<$data_rec> as parameters.
The record header may optionally contain the leading header character, a 
C<<gt>> symbol.  Both parameters are parsed.  A new object instance is
returned on success.  If parsing fails or a record cannot be created, this
method returns undefined.

=cut


   sub new($$) {
      my $pkg = shift;
      my $header = shift;
      my $data_rec = shift;
      my $self = undef;
      my $identifier = undef;

      if ( ( defined ($header) )  &&
           ( $header =~ s[^($UNBOUND_FASTA_SEPARATOR){0,1}(.*)] # accept no
                         [$UNBOUND_FASTA_SEPARATOR$2] ) &&      # separator
           ( defined ($data_rec) )  &&
           ( defined($identifier = _headerToIdentifier($header)) ) && 
                                                             # 
                      # parsable id.
           ( isValidFASTAdata($data_rec) != 0 ) ) {          # valid data?
         $self = {};
         bless $self, $pkg;

         chomp($header);
         chomp($data_rec);

         $self->{header} = $header;
         $self->{identifier} = $identifier;
         $self->{data_rec} = $data_rec;
         $self->{rec_size} = length($data_rec);              # compute size
      }
      return $self;
   }


=item $result = $obj_instance->equals($fasta_obj);

This method takes a I<TIGR::FASTArecord> object as a parameter.  It compares
the passed object with the calling object's internal structures to test for
equality.  If the two objects are equal, this method returns true (1). 
Otherwise, this method returns false (undefined).

=cut


   sub equals($) {
      my $self = shift;
      my $other_obj = shift;
      my $return_val = undef;

      if ( ( defined ($other_obj) ) &&
           ( $self->toString() eq $other_obj->toString() )
         ) {
         $return_val = 1;
      }
      else {
         $return_val = undef;
      }

      return $return_val;
   }


=item $identifier = $obj_instance->getHeader();

This method returns the header string for the record.

=cut


   sub getHeader() {
      my $self = shift;
      my $header = $self->{header}; # should always be defined
      if (! isValidFASTAheader($header) ) {
         # need to prepend an $UNBOUND_FASTA_SEPARATOR
         $header = $UNBOUND_FASTA_SEPARATOR . $header;
      }
      return $header; 
   }


=item $identifier = $obj_instance->getIdentifier();

This method returns the identifier string for the record.

=cut


   sub getIdentifier() {
      my $self = shift;
      return $self->{identifier};  # should always be defined
   }


=item $data_contents = $obj_instance->getData();

This method returns the data contents of the record as an 
uninterrupted string.

=cut


   sub getData() {
      my $self = shift;
      return $self->{data_rec};  # should always be defined
   }


=item $size_of_rec = $obj_instance->size();

This method returns the size of the data string contained by the
record.

=cut


   sub size() {
      my $self = shift;
      return $self->{rec_size};  # should always be defined
   }


=item $str_rep = $obj_instance->toString();

This method returns a string representation of the FASTA record.
This string representation conforms to the TIGR definition of a
FASTA record.

=cut


   sub toString() {
      my $self = shift;
      my $data_copy = $self->{data_rec};
      my $str_rep = $self->getHeader() . "\n";
      my $seg = undef;
      my @segments = ();
      my $last_seg = undef;
      my $skip = undef;

      while((defined($seg = substr $data_copy, 0, 
                $$OUTPUT_LINE_LENGTH_REF,'')) && ($seg ne '')) {
         $str_rep .= $seg . "\n";
      }
      return $str_rep;
   }

=item $rev_compl_str = $obj_instance->reverseComplement($NA_strand);

This method takes in a string which represents a Nucleotide strand and
returns the reverse complement of the strand. If the string does not 
represent a Nucleotide strand or is undefined, the method returns undefined. 
When an empty string is passed into the method, we get an empty string on 
return.

=cut

   
   sub reverseComplement($) {
      my $self = shift;
      my $NA_strand = shift;
      my $error_condition = 1;
      my $rev_compl_str = undef;
      my $validData = 0;
      
      if((defined ($NA_strand)) && 
         ($validData = _isNucleotideData($NA_strand)) && ($validData == 1)) {
          
         $rev_compl_str = '';
         # taking the complement of the sequence
         if(!defined( $NA_strand =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBbXxNn.-/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVvXxNn.-/)) {
            $error_condition = undef;  
         }
         # taking the reverse of the complemented sequence
         elsif(!defined($rev_compl_str = reverse($NA_strand))) {
            $error_condition = undef;
	 }
      }
      else {
         $error_condition = undef;
      }
      return ( defined ( $error_condition ) ) ? $rev_compl_str : undef;
   }


=item $rev_compl_str = $obj_instance->reverseComplementData();

This method returns the reverse complement of the FASTA record data. If the 
FASTA record data does not represent a Nucleotide strand or is undefined, the 
method returns undefined.

=cut


   sub reverseComplementData() {
      my $self = shift;
      my $NA_strand = $self->{data_rec};
      my $error_condition = 1;
      my $rev_compl_str = undef;
      my $validData = 0;
      
      if((defined ($NA_strand)) && 
         ($validData = _isNucleotideData($NA_strand)) && ($validData == 1)) {
          
         $rev_compl_str = '';
         # taking the complement of the sequence
         if(!defined( $NA_strand =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBbXxNn.-/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVvXxNn.-/)) {
            $error_condition = undef;  
         }
         # taking the reverse of the complemented sequence
         elsif(!defined($rev_compl_str = reverse($NA_strand))) {
            $error_condition = undef;
	 }
      }
      else {
         $error_condition = undef;
      }
      return ( defined ( $error_condition ) ) ? $rev_compl_str : undef;
   }


=item $data_substr = $obj_instance->subSequence($startpos, $length);

This method behaves like substr(). This method takes in two numbers $startpos 
and $length and returns a substring of the record data. The $startpos for the 
first base in the data is 0. The $length is optional. The substring is 
extracted starting at $startpos characters from the front of the data string. 
If $startpos is negative, the substring starts that far from the end of the 
string instead. If $length is omitted, everything to the end of the string is 
returned. If $length is negative, the length is calculated to leave that many 
characters off the end of the string. Otherwise, $length indicates the length 
of the substring to extract. If $startpos is undefined the function returns 
undefined. If either $startpos or $length are greater than the data length, the
method returns undefined

=cut

   
   sub subSequence($;$) {
       my $self = shift;
       my ($startpos, $length) = @_;
       my $data_str = $self->getData();
       my $data_substr = undef;
       my $strlen = length($data_str);
       
       if((defined ($data_str)) && (defined ($startpos)) && 
         (defined ($length)) && (abs($startpos) < $strlen) && 
	 (abs($length) <= $strlen)) { 
          #both parameters are specified
          $data_substr = substr($data_str, $startpos, $length);
       }
       elsif((defined ($data_str)) && (defined ($startpos)) && 
             (!defined ($length)) && (abs($startpos) < $strlen) ) {
          #only one parameter is specified.
          $data_substr = substr($data_str, $startpos);  
       }
       else {
          #either $data_str is undefined or $startpos is undefined so do 
          #nothing.
       }
       return $data_substr;
   }

=item $result = $obj_instance->unGap();

This method removes all the gaps('-' characters) from the record data. The 
record size is changed after the gaps are removed. The method returns 1 on 
success and undefined otherwise.

=cut

   
   sub unGap() {
      my $self = shift;
      my $data_str = $self->getData();
      my $result = 1;
      if((defined ($data_str)) && ($data_str =~ /[-]+/)) {
         #removing the gaps in the data.
         $data_str =~ s/[-]+//g;
         $self->{data_rec} = $data_str;
         $self->{rec_size} = length($data_str);
         
      }
      elsif(!defined ($data_str)) {
	  $result = undef;
      }
      return $result;
   }

=back



=head1 USAGE

To use this module, load it using the C<use> function.  The object must
be initialized with a header and a data string.  An example follows.
Please refer to the C<TIGR::FASTAreader> package usage for more examples.

   #!/usr/local/bin/perl -w

   use strict;
   use TIGR::FASTArecord;

   MAIN:
   {

      # set up a simple example without using the leading carot.
      my $header1 = "ORF00001 The first ORF in the series";
      my $data1 = 
         "MEEISTPEGGVLVPISIETEVKRAYIDYSMSVIVSRALPDVRDGLKPVHRRILYAMEEKG" .
         "LRFSGPTRKCAKIVGDVLGSFHPHGDASVYDALVRLGQDFSLRYPVIHPQGNFGTIGGDP" .
         "PAAYRYTEAKMARIAESMVEDIKKETVSFVPNFDDSDVEPTVLPGRFPFLLANGSSGIAV" .
         "GMTTNMPPHNLREIAAAISAYIENPNLSIQELCDCINGPDFPTGGIIFGKNGIRQSYETG" .
         "RGKIVVRARFTIETDSKGRDTIIFTEVPYQVNTTMLVMRIGELARAKVIEGIANVNDETS" .
         "DRTGLRIVVELKKGTPAQVVLNHLFAKTPLQSSFNVINLALVEGRPRMLTLKDLVRYFVE" .
         "HRVDVVTRRAHFELRKAQERIHLVRALIRALDAIDKIITLIRHSQNTELAKQRLREQFDF" .
         "DNVQAQAIVDMQMKRLTGLEVESLRTELKDLTELISSLEELLTSPQKVLGVVKKETRDIA" .
         "DMFGDDRRTDIVSNEIEYLDVEDFIQKEEMVILISHLGYIKRVPVSAYRNQNRGGKGSSS" .
         "ANLAAHDFISQIFTASTHDYVMFVTSRGRAYWLKVYGIPESGRANRGSHIKSLLMVATDE" .
         "EITAIVSLREFSNKSYVFMATARGVVKKVTTDNFVNAKTRGIIALKLSGGDTLVSAVLVQ" .
         "DEDEVMLITRQGKALRMSGREVREMGRNSSGVIGIKLTSEDLVAGVLRVSEQRKVLIMTE" .
         "NGYGKRVSFSEFSVHGRGTAGQKIYTQTDRKGAIIGALAVLDTDECMCITGQGKTIRVDV" .
         "CAISVLGRGAQGVRVLDIEPSDLVVGLSCVMQG";

      my $fasta_record = new TIGR::FASTArecord $header1, $data1;
      if ( defined ( $fasta_record ) ) {
         print STDOUT "Sequence " . $fasta_record->getIdentifier() . " is " .
                      $fasta_record->size() . " residues.\n";
         print STDOUT $fasta_record->toString();
      }
      else {
         die "Invalid FASTA record 1";
      }
         
      # but this header is also valid.
      my $header2 = ">ORF00001 The first ORF in the series";
      my $data2 = 
         "MEEISTPEGGVLVPISIETEVKRAYIDYSMSVIVSRALPDVRDGLKPVHRRILYAMEEKG" .
         "LRFSGPTRKCAKIVGDVLGSFHPHGDASVYDALVRLGQDFSLRYPVIHPQGNFGTIGGDP" .
         "PAAYRYTEAKMARIAESMVEDIKKETVSFVPNFDDSDVEPTVLPGRFPFLLANGSSGIAV" .
         "GMTTNMPPHNLREIAAAISAYIENPNLSIQELCDCINGPDFPTGGIIFGKNGIRQSYETG" .
         "RGKIVVRARFTIETDSKGRDTIIFTEVPYQVNTTMLVMRIGELARAKVIEGIANVNDETS" .
         "DRTGLRIVVELKKGTPAQVVLNHLFAKTPLQSSFNVINLALVEGRPRMLTLKDLVRYFVE" .
         "HRVDVVTRRAHFELRKAQERIHLVRALIRALDAIDKIITLIRHSQNTELAKQRLREQFDF" .
         "DNVQAQAIVDMQMKRLTGLEVESLRTELKDLTELISSLEELLTSPQKVLGVVKKETRDIA" .
         "DMFGDDRRTDIVSNEIEYLDVEDFIQKEEMVILISHLGYIKRVPVSAYRNQNRGGKGSSS" .
         "ANLAAHDFISQIFTASTHDYVMFVTSRGRAYWLKVYGIPESGRANRGSHIKSLLMVATDE" .
         "EITAIVSLREFSNKSYVFMATARGVVKKVTTDNFVNAKTRGIIALKLSGGDTLVSAVLVQ" .
         "DEDEVMLITRQGKALRMSGREVREMGRNSSGVIGIKLTSEDLVAGVLRVSEQRKVLIMTE" .
         "NGYGKRVSFSEFSVHGRGTAGQKIYTQTDRKGAIIGALAVLDTDECMCITGQGKTIRVDV" .
         "CAISVLGRGAQGVRVLDIEPSDLVVGLSCVMQG";

      my $fasta_record2 = new TIGR::FASTArecord $header2, $data2;
      if ( defined ( $fasta_record2 ) ) {
         print STDOUT "Sequence " . $fasta_record2->getIdentifier() . " is " .
                      $fasta_record2->size() . " residues.\n";
         print STDOUT $fasta_record2->toString();
      }
      else {
         die "Invalid FASTA record 2";
      }

      # this entry fails; note the 'J' in the second line of data (8th char.)
      my $header3 = "ORF00001 The first ORF in the series";
      my $data3 = 
         "MEEISTPEGGVLVPISIETEVKRAYIDYSMSVIVSRALPDVRDGLKPVHRRILYAMEEKG" .
         "LRFSGPTJKCAKIVGDVLGSFHPHGDASVYDALVRLGQDFSLRYPVIHPQGNFGTIGGDP" .
         "PAAYRYTEAKMARIAESMVEDIKKETVSFVPNFDDSDVEPTVLPGRFPFLLANGSSGIAV" .
         "GMTTNMPPHNLREIAAAISAYIENPNLSIQELCDCINGPDFPTGGIIFGKNGIRQSYETG" .
         "RGKIVVRARFTIETDSKGRDTIIFTEVPYQVNTTMLVMRIGELARAKVIEGIANVNDETS" .
         "DRTGLRIVVELKKGTPAQVVLNHLFAKTPLQSSFNVINLALVEGRPRMLTLKDLVRYFVE" .
         "HRVDVVTRRAHFELRKAQERIHLVRALIRALDAIDKIITLIRHSQNTELAKQRLREQFDF" .
         "DNVQAQAIVDMQMKRLTGLEVESLRTELKDLTELISSLEELLTSPQKVLGVVKKETRDIA" .
         "DMFGDDRRTDIVSNEIEYLDVEDFIQKEEMVILISHLGYIKRVPVSAYRNQNRGGKGSSS" .
         "ANLAAHDFISQIFTASTHDYVMFVTSRGRAYWLKVYGIPESGRANRGSHIKSLLMVATDE" .
         "EITAIVSLREFSNKSYVFMATARGVVKKVTTDNFVNAKTRGIIALKLSGGDTLVSAVLVQ" .
         "DEDEVMLITRQGKALRMSGREVREMGRNSSGVIGIKLTSEDLVAGVLRVSEQRKVLIMTE" .
         "NGYGKRVSFSEFSVHGRGTAGQKIYTQTDRKGAIIGALAVLDTDECMCITGQGKTIRVDV" .
         "CAISVLGRGAQGVRVLDIEPSDLVVGLSCVMQG";

      my $fasta_record3 = new TIGR::FASTArecord $header3, $data3;
      if ( defined ( $fasta_record3 ) ) {
         print STDOUT "Sequence " . $fasta_record3->getIdentifier() . " is " .
                      $fasta_record3->size() . " residues.\n";
         print STDOUT $fasta_record3->toString();
      }
      else {
         die "Invalid FASTA record 3";
      }
   }

=cut

}

1;
