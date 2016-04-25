# $Id: FASTAreader.pm,v 1.1 2004/04/28 15:03:43 aphillip Exp $

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


package TIGR::FASTAreader;
{
   
=head1 NAME

TIGR::FASTAreader - TIGR::FASTAreader class for parsing and navigating
FASTA format files.

=head1 SYNOPSIS

   use TIGR::FASTAreader;  
   my $obj_instance = new TIGR::FASTAreader ($foundation_obj_ref, 
                      $error_array_ref, $fasta_file_name);

=head1 DESCRIPTION

This module iterates over a FASTA formatted database file.  It provides
data extraction and simple analysis routines.  This module utilizes 
acceptance validation of FASTA formatted files via the TIGR::FASTAgrammar
module.

=cut

   BEGIN {
      require 5.006_00;
   }

   use strict;
   use IO::File;
   use TIGR::Foundation;
   use TIGR::FASTAgrammar ':public';
   use TIGR::FASTAgrammar ':private';
   use TIGR::FASTArecord;


   ## internal variables and identifiers

   our $REVISION = (qw$Revision: 1.1 $)[-1];
   our $VERSION = '1.21';
   our $VERSION_STRING = "$VERSION (Build $REVISION)";
   our @DEPEND = ();

   my $SYS_ERR = 0;           # this flag specifies non-user related error
   my $USR_ERR = 1;           # this flag specifies user related error

   ## external variables

   my $UNBOUND_FASTA_SEPARATOR = $TIGR::FASTAgrammar::UNBOUND_FASTA_SEPARATOR;
   
   # debugging scheme
   #
   #   Debugging via the TIGR Foundation uses increasing log levels based on
   #   nesting.  'MAIN' starts at level 1.  Every nest increments the level by
   #   1.  
   #   Subroutines always start nesting at level 2.  As debugging levels
   #   increase, logging is more verbose.  This makes sense as you log at
   #   greater depth (ie. deeper branching).
   #
   #   The following definitions help emphasize the debugging in the program.
   #
   my $DEBUG_LEVEL_1 = 1;
   my $DEBUG_LEVEL_2 = 2;
   my $DEBUG_LEVEL_3 = 3;
   my $DEBUG_LEVEL_4 = 4;
   my $DEBUG_LEVEL_5 = 5;
   my $DEBUG_LEVEL_6 = 6;
   my $DEBUG_LEVEL_7 = 7;
   my $DEBUG_LEVEL_8 = 8;
   my $DEBUG_LEVEL_9 = 9;

   ## prototypes

   sub new(;$$$);
   sub open($;$);
   sub close();
   sub index();
   sub seekIndex($);
   sub getRecordByIdentifier($);
   sub seekIdentifer($);
   sub get();
   sub next();
   sub hasNext();
   sub count();
   sub path();
   sub _initialize();
   sub _parseDBfile();
   sub _nullRecordHandler($$);
   sub _errorHandler($$$);


   ## implementation

=over

=item $obj_instance = new TIGR::FASTAreader ($foundation_object,
                             $error_array_ref, $db_file);

This method returns a new instance of a TIGR::FASTAreader object.  It takes
three optional parameters: a TIGR::Foundation object (C<$foundation_object>),
a reference to an array for logging user error messages (C<$error_array_ref>),
and FASTA file (C<$db_file>).  The new instance is returned on success.  If
the file supplied cannot be opened or is invalid, this method returns 
undefined. This method also returns undefined if the parameters supplied are 
invalid. Errors in parsing are written to the array at C<$error_array_ref> 
and the log file.

=cut

   sub new(;$$$) {
      my $pkg = shift;
      my @method_args = @_;

      my $error_condition = 0;
      my $self = {};
      bless $self, $pkg;
      $self->_initialize();  # set up internal variables;

      if ( ( scalar (@method_args) > 0 ) &&
           ( ( ref ($method_args[0]) ) =~ /foundation/i ) ) {
         $self->{foundation} = shift @method_args;
         $self->_errorHandler("Got TIGR::Foundation in new()", $DEBUG_LEVEL_3, 
                               $SYS_ERR);
      }
      else {
         $self->{foundation} = undef;
         $self->_errorHandler("No TIGR::Foundation in new()", $DEBUG_LEVEL_3, 
                               $SYS_ERR);
      }

      if ( ( scalar (@method_args) > 0 ) &&
           ( ( ref ($method_args[0]) ) =~ /array/i ) ) {
         $self->{error_ref} = shift @method_args;
         $self->_errorHandler("Got Error ARRAY in new()", $DEBUG_LEVEL_3, 
                               $SYS_ERR);
      }
      else {
         $self->{error_ref} = undef;
         $self->_errorHandler("No Error ARRAY in new()", $DEBUG_LEVEL_3, 
                               $SYS_ERR);
      }

      if ( ( scalar (@method_args) > 0 ) &&
           ( ! ref ($method_args[0]) ) ) {
         my $filename = shift @method_args;
         if(defined($filename)) {
            $self->{db_file_name} = $filename ;
            $self->_errorHandler("Got file name in new()", $DEBUG_LEVEL_4, 
                                  $SYS_ERR);
         }
         else {
            $self->_errorHandler("undef passed as filename", $DEBUG_LEVEL_4, 
				  $USR_ERR);
	 }
      }
      else {
         $self->{db_file_name} = undef;
         $self->_errorHandler("No file name in new()", $DEBUG_LEVEL_3, 
                               $SYS_ERR);
      }

      # check for invocation errors
      if ( ( scalar (@method_args) > 0 ) ) {
         $error_condition = 1;
         $self->_errorHandler("Too many parameters passed to new() method",
                               $DEBUG_LEVEL_3, $SYS_ERR);
      }
      elsif (   defined ( $self->{db_file_name} ) &&
           ! defined ( $self->open($self->{db_file_name}, "r") ) ) {
         # the error message is logged via the open() routine
         $self = undef;
      }
      return ( $error_condition == 0 ) ? $self : undef;
   }


=item $result = $obj_instance->open($file_name, $flag);

This method opens a FASTA file for reading.  This method also parses the file
for correctness.  The file, C<$file_name>, is opened using the C<open()> flags
specified by C<$flag>.  On success, this method returns 1.  If the file cannot
be opened or parsing fails, this method returns undefined.

=cut

   sub open($;$) {
      my $self = shift;
      my $db_file_name = shift;
      my $open_flags = shift;

      my $error_condition = 0;

      if ( ( ! defined ($open_flags) ) ||
           ( $open_flags !~ /^r$/i ) ) {
         $open_flags = "r";
      }
      $self->_errorHandler("Open flags = \'$open_flags\' in open()", 
                            $DEBUG_LEVEL_3, $SYS_ERR);

      # close a previously open file
      if ( defined ($self->{db_handle}) ) {
         $self->_errorHandler("Closing old handle in open()", $DEBUG_LEVEL_3, 
                               $SYS_ERR);
         $self->close();
      }

      if (!(
            ( defined ( $db_file_name ) ) &&
            ( $self->{db_file_name} = $db_file_name ) &&
            ( defined ( $self->{db_file_name} )) &&
            ( defined ( $self->{db_handle} = 
                     new IO::File $self->{db_file_name}, $open_flags ))
         ) ) {
         $error_condition = 1;
         $self->_errorHandler(
            "Cannot open file \'$self->{db_file_name}\'", $DEBUG_LEVEL_3, 
             $USR_ERR);
      }
      elsif ( ( defined ( $self->{db_file_name} ) ) &&
              ( defined ( $self->{db_handle} ) ) &&
              ( $self->_parseDBfile() == 0 ) ) {
         $error_condition = 1;
         $self->_errorHandler("Encountered errors in file " .
            "\'$self->{db_file_name}\'.", $DEBUG_LEVEL_3, $USR_ERR);
      }

      if ( $error_condition == 1 ) {
         $self->_initialize(); # reset object state
      }

      return ($error_condition == 1) ? undef : 1;
   }


=item $result = $obj_instance->close();

This method closes the object file stream and resets all internal data
structures.  The result of the operation is returned.  If the file stream
is closed successfully, this object returns true (1), otherwise false
(undefined).

=cut

   sub close() {
      my $self = shift;
      my $return_val = undef;

      if ( defined ( $self->{db_handle} ) ) {
         $return_val = $self->{db_handle}->close();
         if (!$return_val) {
            $return_val = undef;
            $self->_errorHandler(
               "Error closing FASTA file: $self->{db_file_name}", 
                $DEBUG_LEVEL_4, $USR_ERR);
         }
      }
      $self->_initialize();
      return $return_val;
   }


=item $record_num = $obj_instance->index();

This method returns the record number of the active record.  If no record has 
been selected (ie. made active), then this method returns undefined.  If
the active record pointer is before the first record, this method returns
'-1'.

=cut

   sub index() {
      my $self = shift;
      my $return_val = undef;

      if ( defined ($self->{active_record}) ) {
         $return_val = $self->{active_record};
      }
      else {
         $return_val = undef;
      }
      return $return_val;
   }


=item $result = $obj_instance->seekIndex($num);

This method selects a record by record order index.  The C<$num> ordered
record is selected.  If C<$num> is out of range for the database or not -1
(indicating to seek one record before the first record), this
function returns undefined and the active record pointer is not changed.
Otherwise, the requested record is made active and the method returns 1.

=cut

   sub seekIndex($) {
      my $self = shift;
      my $active_num = shift;
      my $return_val;

      if ( (defined ($active_num) ) &&
           ( ($active_num =~ /^\d+$/) ||
             ($active_num == -1 ) ) &&
           ($active_num < $self->count()) ) {
         $self->_errorHandler(
            "Setting active record num to $active_num.", 
             $DEBUG_LEVEL_3, $SYS_ERR);
         $self->{active_record} = $active_num;
         $return_val = 1;
      }
      else {
         if ( ! defined ($active_num) ) {
            $active_num = "<undefined>";
         }
         $self->_errorHandler(
            "Cannot set active record num to $active_num, out of " .
            $self->count(), $DEBUG_LEVEL_3, $SYS_ERR);
         $return_val = undef;
      }
      return $return_val;
   }
 

=item $result = $obj_instance->next();

This method selects the next record in numerical order to be the active
record.  It returns the record on success, undefined on failure.  If the active
record is equal to -1, the first record is selected.

=cut

   sub next() {
      my $self = shift;
      my $return_val = undef;

      if ( defined ( $self->hasNext() ) ) {
         $self->{active_record}++;
         $return_val = $self->get();
      }
      # set undefined if no more records left or no records at all
      else {
         $return_val = undef;
      }

      return $return_val;
   }


=item $result = $obj_instance->hasNext();

This method returns true (1) if there are more elements beyond the
current element.  If not, this method returns false (undefined).

=cut

   sub hasNext() {
      my $self = shift;
      my $return_val = undef;

      if ( ( defined ($self->{active_record}) ) &&
           ( $self->{active_record} >= -1 ) &&
           ( $self->{active_record} < ( $self->count() - 1 ) )
         ) {
         $return_val = 1;
      }
      else {
         $return_val = undef;
      }
      return $return_val;
   }


=item $result = 
         $obj_instance->getRecordByIdentifier($identifier);

This method selects a record by record minimal identifier.
If C<$identifier> does not exist in the set of records, this function
returns undefined and the previously active record remains active.  Otherwise,
the requested record is made active and the method returns a 
C<TIGR::FASTArecord> object representation of the current(active) record.

=cut

   sub getRecordByIdentifier($) {
    
      my $self = shift;
      my $identifier = shift;
      my $fasta_record = undef;
      my $seek_result = undef;
      if((defined $identifier)) {
         
         my $seek_result = $self->seekIdentifier($identifier);
         if( (defined ($seek_result)) && ($seek_result == 1)) {
	    $fasta_record = $self->get();
         }   
      }
      else {
         $self->_errorHandler("undefined identifier passed", $DEBUG_LEVEL_3, 
			       $USR_ERR);
      }
      return $fasta_record;
   }
 

=item $result = $obj_instance->seekIdentifier($identifier);

This method selects a record by record minimal identifier.
If C<$identifier> does not exist in the set of records, this function
returns undefined and the previously active record remains active.  Otherwise,
the requested record is made active and the method returns 1.

=cut

   sub seekIdentifier($) {
      my $self = shift;
      my $identifier = shift;
      my $number = undef;

      if ( (defined $identifier) &&
           (exists $self->{identifier_to_number_hash}->{$identifier}) ) {
         $number = $self->{identifier_to_number_hash}->{$identifier};
         $self->_errorHandler(
            "Got record number $number for $identifier.", $DEBUG_LEVEL_3, 
             $SYS_ERR);
      }
      return $self->seekIndex($number);
   }
      

=item $record_contents = $obj_instance->get();

This method returns a C<TIGR::FASTArecord> object representation of the current
(active) record.  If no defined record is active, this method returns 
undefined.

=cut

   sub get() {
      my $self = shift;
      my $db_handle = defined ($self->{db_handle}) ?
         $self->{db_handle} : undef;
      my $header = "";
      my $data = "";
      my $f_obj = undef;
      my $pos_str = undef;
      my @pos_arr = ();
      my $record = undef;

      # open FASTA file for reading
      if (! defined ($db_handle) ) {
         $self->_errorHandler("db_handle not defined: " . 
            "cannot access FASTA file \'$self->{db_file_name}\'.", 
             $DEBUG_LEVEL_3);
      }

      # search and extract the FASTA record
      if (( defined ( $self->{active_record} )) &&
          ( $self->{active_record} > -1 ) &&
          ( $self->{active_record} < $self->count() ) &&
          ( defined ( $db_handle )) &&
          ( defined ( $pos_str = 
                      $self->{number_to_fp_array}->{$self->{active_record}} ))
         ) {
	 @pos_arr = split " ", $pos_str;

         # seek to the start position of the record in the file and read the 
         # record information
         if((defined $pos_arr[0]) && 
            (defined($db_handle->seek($pos_arr[0],SEEK_SET))) &&
            (defined $pos_arr[1]) && 
            (defined (read $db_handle, $record,
                                          (($pos_arr[1]-$pos_arr[0])+1)))) {
	    
	    ($header) = $record =~ /(>.*)\n/;
            $record =~ s/(>.*)\n//;
            $record =~ s/[\s\n]+//g;
            $f_obj = new TIGR::FASTArecord $header, $record;
         }
      }
      return $f_obj;
   }
      

=item $db_name = $obj_instance->path();

This method returns the path to the file used for processing.

=cut

   sub path() {
      my $self = shift;
      # the existence of db_file_name is checked in new()
      return $self->{db_file_name};
   }


=item $cnt = $obj_instance->count();

This method returns the number of records in the database file.

=cut

   sub count() {
      my $self = shift;
      return $self->{num_records};
   }


# $obj_instance->_initialize();

#This method resets the object to its initial state.  Internal data structures
#are reset.  This method does not return.

   sub _initialize() {
      my $self = shift;

      $self->{num_records} = 0;                  # number of recs

      # look up methods for records here
      # the active record is stored as a sequence number
      $self->{active_record} = undef;            # current working record
      $self->{number_to_fp_array} = ();          # map seq# to file loc
      $self->{number_to_identifier_array} = ();  # map seq# to identifier
      $self->{identifier_to_number_hash} = ();   # map seq identifier to seq#
      $self->{error_cnt} = 0;                    # parse error tabulator
      $self->{db_file_name} = "";
      $self->{db_handle} = undef;
   }


# $obj_instance->_parseDBfile();

#This method parses the FASTA database file passed via the C<open()> method.
#It defines all of the sequence look-ups and validates every record.  This
#method finds the number of sequences and maximum sequence length.  This
#method is called from the C<open()> method.  The active record is un-selected
#by this method.

   sub _parseDBfile() {
      my $self = shift;
      my $last_line_length_lt_std_flag = 0;
      my $line_number = 0;
      my $record_identifier = "<undefined>";
      my $preceding_header_flag = 0;
      my $first_data_line_length = undef;
      my $empty_line_found = 0;
      my $db_handle = defined ( $self->{db_handle} ) ? 
         $self->{db_handle} : undef;
      #the file position where the record starts
      my $pos1 = 0;
      #the file position where the record ends
      my $pos2 = undef;
      #the start and end positions of a record separated with a space
      my $string = "";
      # variable to give the length of each line
      my $line_len = undef;
      # the sum of the length of all the lines in a file  
      my $sum_len = 0;
      # loop through FASTA file
      while ( ( defined ( $db_handle ) ) &&
              ( defined ( my $line = <$db_handle> ) ) &&
              ( ++$line_number ) ) {
         chomp $line;
         $line_len = length($line);
         $sum_len += $line_len;
         $sum_len++;
         
         # check FASTA data
         if ( ( defined ( $record_identifier ) ) &&
              ( $record_identifier !~ /<undefined>/ ) &&
              ( ( isValidFASTAdata($line) ) != 0 ) ) {
	    
            # check if previous line was empty
            if ( $empty_line_found == 1 ) {
               $self->{error_cnt}++;
               $self->_errorHandler("ERROR: Empty line found at line ".
                                    ($line_number - 1). " - empty lines are ".
                                    "allowed only at the end of a file", 
                                    $DEBUG_LEVEL_5, $USR_ERR);
               $empty_line_found = 0; 
            }
            
	    if($preceding_header_flag == 1) {
	       $first_data_line_length = setValidFASTAlineLength($line);
            }
	   
            # check $last_line_length_lt_std_flag for an error on previous line
	    if(defined ($first_data_line_length)) {
               if ( $last_line_length_lt_std_flag == 1 ) {
                  $self->{error_cnt}++;
                  $self->_errorHandler("Expected: FASTA data definition " .
                     "lines should be $first_data_line_length bases " .
                     "(characters) across. Only the last line of a sequence ".
                     "data definition may be less than " .
                     "$first_data_line_length bases (characters) " .
                     "across, if applicable.  See line " . 
                     ($line_number - 1) . '.', $DEBUG_LEVEL_5, $USR_ERR);
               }
               $last_line_length_lt_std_flag = 0;
	    
               # check current line for over-length problem
               if ( $line_len > $first_data_line_length ) {
                  $self->{error_cnt}++;
                  $self->_errorHandler("Expected: FASTA data definition " .
                     "lines should be $first_data_line_length bases " .
                     "(characters) across. Only the last line of a sequence ".
                     "data definition may be less than " .
                     "$first_data_line_length bases (characters) ".
                     "across, if applicable.  See line " . $line_number . '.',
                      $DEBUG_LEVEL_5,$USR_ERR);
               }
            
               #check current line for under-length problem; report only if not
               #the last line in the data definition
               elsif ( $line_len < $first_data_line_length ) {
                  $last_line_length_lt_std_flag = 1;
               }
	    }
            $preceding_header_flag = 0;
         }
         # check for FASTA header
         elsif ( ( isValidFASTAheader($line) ) != 0 ) {
            if ( ! defined ( $self->{active_record} ) ) {
               $self->{active_record} = -1;
            }
            
            $self->{active_record}++;
            
            if( (defined $line_number) &&
                ($line_number > 1) ) {
               $pos2 = (($sum_len - $line_len)-2);  
	       $pos1 = $pos2+1;

               if((defined $pos1) && ( defined $string)) { 
                  $string .= "$pos2";
                  # store the start and end of a fasta record in a hash
                  $self->{number_to_fp_array}->{($self->{active_record})-1} = 
                     $string;
	       }
            }
     
            if(defined $pos1) {
               $string = "$pos1 ";
            }
 
            # check if previous line was a FASTA header
            if ( $preceding_header_flag == 1 ) {
               $self->_nullRecordHandler($self->{active_record} - 1,
                  $line_number);
            }
            
            # check if previous line was empty
            if ( $empty_line_found == 1 ) {
               $self->{error_cnt}++;
               $self->_errorHandler("ERROR: Empty line found at line ".
                                    ($line_number - 1). " - empty lines are ".
                                    "allowed only at the end of a file", 
                                    $DEBUG_LEVEL_5, $USR_ERR);
               $empty_line_found = 0;
            }
            # if it's a valid FASTA header, then don't need to check again
            # extract the record IDENTIFIER
            $record_identifier = _headerToIdentifier($line); 
            if ( defined (
                 $self->{identifier_to_number_hash}->{$record_identifier} ) ) {
               $self->{error_cnt}++;
               $self->_errorHandler("Expected: unique FASTA " .
                  "identifier.  \'$record_identifier\' is a duplicate at " .
                  "line $line_number.", $DEBUG_LEVEL_5, $USR_ERR);
            }
            else {
               $self->{identifier_to_number_hash}->{$record_identifier} =
                  $self->{active_record};
               $self->{number_to_identifier_array}->{$self->{active_record}} =
                  $record_identifier;
            }

            # set up the variables for parsing a new record
            $last_line_length_lt_std_flag = 0;
            $preceding_header_flag = 1;
            $self->{num_records}++;
         }
         # handle empty space
         # empty space after the last record is allowed
         elsif($line eq "") {
            $empty_line_found = 1;
	    next;
	 }
         # handle error data types
         else {
	    $self->{error_cnt}++;
            # check if previous line was empty
            if ( $empty_line_found == 1 ) {
               $self->{error_cnt}++;
               $self->_errorHandler("ERROR: Empty line found at line ".
                                    ($line_number - 1). " - empty lines are ".
                                    "allowed only at the end of a file", 
                                    $DEBUG_LEVEL_5, $USR_ERR);
               $empty_line_found = 0;
            }
            
            # line has a separator token in it, so it may be header
            if ( $line =~ /$UNBOUND_FASTA_SEPARATOR/ ) {
               $self->_errorHandler("Expected: record header " .
                  "information in FASTA record header.  Got: \'$line\' at " .
                  "line $line_number.", $DEBUG_LEVEL_6, $USR_ERR);
               $last_line_length_lt_std_flag = 0;
            }
            # if last data line was small, expect this to be a header too
            elsif ( $last_line_length_lt_std_flag == 1 ) {
               $self->_errorHandler("Expected: FASTA record header " .
                  "beginning with \'>\'.  Got: \'$line\' at line ".
                  "$line_number.",$DEBUG_LEVEL_6, $USR_ERR);
               $last_line_length_lt_std_flag = 0;
            }
            elsif ( ( defined ( $record_identifier ) ) && 
                    ( $record_identifier !~ /<undefined>/ ) ) {
               $self->_errorHandler("Expected: valid FASTA data " .
                  "definition for record identifier \'$record_identifier\'. " .
                  "Check sequence content at line $line_number for invalid " .
                  "bases (data type: invalid data).", $DEBUG_LEVEL_6, 
                   $USR_ERR);
            }
            else {
               $self->_errorHandler("Expected: FASTA record header " .
                  "followed by definition of sequence.  Invalid input at " .
                  "line $line_number.", $DEBUG_LEVEL_6, $USR_ERR);
            }
         }
      }  # end while
      
      $pos2 = $sum_len;
      if((defined $pos1) && ( defined $string)) {  
         $string .= "$pos2";
         # store the start and end of a fasta record in a hash
         $self->{number_to_fp_array}->{$self->{active_record}} = $string;
      }
      
      # check terminal case data definition
      if ( $preceding_header_flag == 1 ) {
         $self->_nullRecordHandler($self->{active_record}, $line_number);
      }
      
      $self->{active_record} = -1; # set counter to the beginning
      
      return ( $self->{error_cnt} == 0 ) ? 1 : 0;
   }


# $obj_instance->_nullRecordHandler($$);

#This method handles the case of a null or equivalently empty record
#encountered during parsing.  It logs the appropriate message to the 
#TIGR Foundation object.  The only arguments are the record number
#and the line number.

   sub _nullRecordHandler($$) {
      my $self = shift;
      my $active_num = shift;
      my $line_number = shift;
      my $preceding_rec_line_number = undef;
      my $record_identifier = undef;

      if ( defined ( $self->{number_to_identifier_array}->
              {$active_num} ) ) {
         $record_identifier = $self->{number_to_identifier_array}->
              {$active_num};
      }
      else {
         $record_identifier = "<unknown>";
      }
      
      if(defined $line_number) {
         if((defined $active_num) && 
           ($active_num < ($self->{active_record}))) {
            $preceding_rec_line_number = $line_number-1;
	 }
         elsif((defined $active_num) && 
               ($active_num == ($self->{active_record}))) {
            $preceding_rec_line_number = $line_number;
	 }
      }
      else {
         $preceding_rec_line_number = "<unknown>";
      }

      $self->{error_cnt}++;
      if ( $self->{db_handle}->eof() == 1 ) {
         $self->_errorHandler("Expected: FASTA record header " .
            "followed by definition of sequence.  Record identifier " .
            "\'" . $record_identifier . "\' is undefined from line " .
            $preceding_rec_line_number . ".  Got end of file after line " . 
            $line_number . ".", $DEBUG_LEVEL_5, $USR_ERR);
      }
      else {
         $self->_errorHandler("Expected: FASTA record header " .
            "followed by definition of sequence.  Record identifier " .
            "\'" . $record_identifier . "\' is undefined from line " .
            $preceding_rec_line_number . ".  Got FASTA header at line " . 
            $line_number . ".", $DEBUG_LEVEL_5, $USR_ERR);
      }
   }


# $message = $obj_instance->_errorHandler($message, $tf_level,
# $internal_log_flag);

#This method handles logging to the TIGR::Foundation module and
#internal error record reference array.  The C<$message> argument is logged
#to the appropriate service.  The C<$tf_level> parameter specifies the
#logging level for TIGR::Foundation, while the C<$internal_log_flag> parameter
#specifies if C<$message> should be written to the internal array reference
#specified in C<new()>.  If a TIGR::Foundation instance does not exist,
#no logging to that facility occurs.  This method returns C<$message>.

   sub _errorHandler($$$) {
      
      my $self = shift;

      my ( $message, $tf_level, $log_facility ) = @_;

      if ( defined ($message) &&
           defined ($tf_level) &&
           defined ($log_facility) ) {

         if ( defined ($self->{foundation}) ) {
            if ( $log_facility != $USR_ERR ) { # all user errors go to .error
               $self->{foundation}->logLocal($message, $tf_level);
            }
            else {
               $self->{foundation}->logError($message);
            }
         }

         if ( ( defined ($self->{error_ref}) ) &&
              ( $log_facility == $USR_ERR ) ) {
            push @{$self->{error_ref}}, $message;
         }
      }
      return $message;
   }

=head1 USAGE

To use this module, load the C<TIGR::FASTAreader> package via the
C<use> function.  Then, create a new instance of the object via the
C<new()> method, as shown below.  There are several invocations possible
for this method since all parameters to C<new()> are optional.
To access records from the C<TIGR::FASTAreader> instance, the 
C<TIGR::FASTArecord> package must be loaded via the C<use> function.
An example script using this module follows.  The C<TIGR::Foundation>
module is included for completeness but does not have to be used.

   #!/usr/local/bin/perl -w

   # This script accepts FASTA files with the '-i' option
   # on the command line and validates every one in turn.
   # Parse errors are collected to the '@errors_list' array.
   # This program concatenates all of the records together to 
   # one output file specified with the '-o' option.
   # NOTE: The '-i' option must be specified before every input file.
   # NOTE: The 'TIGR::FASTAwriter' module is intended for writing 
   #       FASTA records.

   use strict;
   use TIGR::FASTAreader;
   use TIGR::FASTArecord;

   MAIN:
   {
      my $tf_object = new TIGR::Foundation;
      my @errors_list = ();
      my @input_files = ();
      my $output_file = undef;

      # Capture the return code from the TIGR::Foundation method
      my $result = $tf_object->TIGR_GetOptions('i=s' => \@input_files,
                                               'o=s' => \$output_file);
      if ( $result != 1 ) {
         $tf_object->bail("Invalid command line options.");
      }

      # Create a TIGR::FASTAreader instance using TIGR::Foundation and
      # an error message list.
      my $fasta_reader = new TIGR::FASTAreader $tf_object, \@errors_list;

      if ( !(  defined ( $output_file ) &&
               open OUTFILE, ">$output_file" ) ) {
         $tf_object->bail("Cannot open output file for writing.");
      }

      foreach my $in_file ( @input_files ) {
         $fasta_reader->open($in_file) or
            $tf_object->logLocal("Cannot open or read file $in_file", 2);

         if ( scalar(@errors_list) > 0 ) { # are there parse errors?
            while ( @errors_list ) { # get the messages from the list
               my $message = shift @errors_list; 
               print STDERR $message, "\n";
            }
         }

         while ( $fasta_reader->hasNext() ) {
            # print each record to OUTFILE
            print OUTFILE $fasta_reader->next()->toString();
         }
      }
   }

=cut

}

1;
