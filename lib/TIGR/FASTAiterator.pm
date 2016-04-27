# $Id: FASTAiterator.pm,v 1.1 2004/04/28 15:03:43 aphillip Exp $

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

package TIGR::FASTAiterator;
{

=head1 NAME

TIGR::FASTAiterator - TIGR::FASTAiterator class for parsing and navigating
FASTA format files and streams. An object of this class can parse FASTA
records from STDIN and from a pipe.

=head1 SYNOPSIS

  use TIGR::FASTAiterator;
  my $obj_instance = new TIGR::FASTAiterator ($foundation_obj_ref,
                                              $error_array_ref,
                                              $fasta_file_name);

=head1 DESCRIPTION

This module iterates over a FASTA formatted file stream.  It provides
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
   our $VERSION = '1.11';
   our $VERSION_STRING = "$VERSION (Build $REVISION)";
   our @DEPEND = 
   (
    "TIGR::Foundation",
    "TIGR::FASTAgrammar",
    "TIGR::FASTArecord",
   );

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
   sub open($);
   sub close();
   sub hasNext();
   sub next();
   sub get();
   sub _initialize();
   sub _parse();
   sub _nullRecordHandler($);
   sub _errorHandler($$$);


   ## implementation

=over

=item $obj_instance = new TIGR::FASTAiterator ($foundation_object,
   $error_array_ref, $db_file);

This method returns a new instance of a TIGR::FASTAiterator object.  It takes
three optional parameters: a TIGR::Foundation object (C<$foundation_object>),
a reference to an array for logging user error messages (C<$error_array_ref>),
and a FASTA file (C<$db_file>) or stream. The filename "-" describes stdin. 
The new instance is returned on success.  If the file supplied cannot be 
opened or is invalid, this method returns undefined. This method also returns 
undefined if the parameters supplied are invalid. Errors in parsing are written
to the array C<$error_array_ref>, the error file and the log file.

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
           ! defined ( $self->open($self->{db_file_name}) ) ) {
         # the error message is logged via the open() routine
         $self = undef;
      }
      return ( $error_condition == 0 ) ? $self : undef;
   }


=item $result = $obj_instance->open($file_name);

This method opens a FASTA file or pipe for reading. It takes in the filename to
be opened. If the file name is "-" the input is taken from stdin. On success, 
this method returns 1.  If the file cannot be opened or parsing fails, this 
method returns undefined.

=cut


   sub open($) {
      my $self = shift;
      my $db_file_name = shift;
     
      my $error_condition = 0;

      # close a previously open file
      if ( defined ($self->{db_handle}) ) {
         $self->_errorHandler("Closing old handle in open()", $DEBUG_LEVEL_3, 
                               $SYS_ERR);
         $self->close();
      }
      my $name = $self->{db_file_name};
     
      if (!(
            ( defined ( $db_file_name ) ) &&
            ( $self->{db_file_name} = $db_file_name ) &&
            ( defined ( $self->{db_file_name} ))
         ) ) {
	  
         $error_condition = 1;
         $self->_errorHandler(
            "File name does not exist", $DEBUG_LEVEL_3, $USR_ERR);
      }
      elsif(!defined ( $self->{db_handle} = 
                     new IO::File $self->{db_file_name})) {
         $error_condition = 1;
         $self->_errorHandler(
            "Cannot open file \'$self->{db_file_name}\'", $DEBUG_LEVEL_3, 
             $USR_ERR); 
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
   

=item $result = $obj_instance->hasNext();

This method returns true (1) if there are more elements beyond the current 
element in the filestream. If not, this method returns false (undefined).

=cut

   sub hasNext() {
      my $self = shift;
      my $next_header = $self->{rec_header};
      $self->_errorHandler(
                "Checking to see if the header of the next record is set", 
                $DEBUG_LEVEL_2, $SYS_ERR);
      my $result = undef;
      my $newline = undef;
      my $line_number = $self->{line_number};
      if ($line_number == 0) {
         $self->_errorHandler(
                "No record has been parsed", 
                $DEBUG_LEVEL_3, $SYS_ERR);
         my $db_handle = defined ( $self->{db_handle} ) ? 
                          $self->{db_handle} : undef;
        
         if(defined $db_handle) {
            $newline = <$db_handle>;
         }
         #reading the first line from a file.
         if(defined($newline)) {
	    $line_number++;
            $self->{line_number} = $line_number;
            $next_header = $newline;
            $self->_errorHandler(
                "Assigned the header of the next record", 
                $DEBUG_LEVEL_4, $SYS_ERR);
            $self->{rec_header} = $next_header;
            $result = 1;
	 }
      }
          
      if((defined ($next_header)) && (($next_header) ne "")) {
         $result = 1;
      }
      return $result;
   }


=item $result = $obj_instance->next();

This method selects the next record in the file stream for parsing. If the 
record parses, it is returned, else the method returns undefined. If there is 
no record in the file stream, the method returns undefined.  

=cut  

   sub next() {
      my $self = shift;
      my $record = undef;
      if(($self->_parse()) == 1) {
          $self->_errorHandler(
                "The record parsed", 
                 $DEBUG_LEVEL_3, $SYS_ERR);
         #obtaining the stored record.
         my $recordarray_ref = $self->{recordinfo};
          
         if( (defined $recordarray_ref) && 
             (( ref ($recordarray_ref) ) =~ /array/i) ) {
	     
            my @recordarray = @$recordarray_ref;
            my $array_length = 0;

            if((defined $recordarray_ref) && 
               ($array_length = @recordarray) && 
               ($array_length > 0)) {
         
               my $header = shift @recordarray;
               $self->_errorHandler(
                  "Got the record header", 
                  $DEBUG_LEVEL_5, $SYS_ERR);
               my $data = undef;
               if ( scalar(@recordarray) > 0 )  {
                  $data = join "", @recordarray;
                  $data =~ s/\n//g; # extract new lines from scalar data
               }
               $self->_errorHandler(
                  "Got the record data", 
                  $DEBUG_LEVEL_5, $SYS_ERR);
               $record = new TIGR::FASTArecord ($header, $data);
               
               if(defined($record)) {
                  $self->_errorHandler(
                     "Created new record", 
                     $DEBUG_LEVEL_6, $SYS_ERR);
	       }
	    }
         }
      }
      return $record;
   }


=item $record_contents = $obj_instance->get();

This method returns the current TIGR::FASTArecord object (active record). If 
the current object (active record) is undefined, this method returns undefined.

=cut
   
   sub get() {
      my $self = shift;
      my $record = undef;
      #obtaining the stored record information.
      my $recordarray_ref = $self->{recordinfo};
      if( (defined $recordarray_ref) && 
        (( ref ($recordarray_ref) ) =~ /array/i) ) {
         my @recordarray = @$recordarray_ref;
         my $array_length = 0;

         if((defined $recordarray_ref) && 
           ($array_length = @recordarray) && 
           ($array_length > 0)) {
         
            my $header = shift @recordarray;
            $self->_errorHandler(
                    "Got the record header", 
                     $DEBUG_LEVEL_4, $SYS_ERR);
            my $data = undef;
            if ( scalar(@recordarray) > 0 )  {
                 $data = join "", @recordarray;
                 $data =~ s/\n//g; # extract new lines from scalar data
            }
            $self->_errorHandler(
                   "Got the record data", 
                    $DEBUG_LEVEL_4, $SYS_ERR);
            $record = new TIGR::FASTArecord ($header, $data);
            if(defined($record)) {
               $self->_errorHandler(
                   "Created new record", 
                    $DEBUG_LEVEL_5, $SYS_ERR);
	    }
	 
         }
      }
      return $record;
   }


# $obj_instance->_initialize();

#This method resets the object to its initial state.  Internal data structures
#are reset.  This method does not return.


   sub _initialize() {
      my $self = shift;
      # look up methods for records here
      $self->{error_cnt} = 0;      # parse error tabulator
      $self->{db_file_name} = "";
      $self->{db_handle} = undef;
      $self->{rec_header} = undef; # the next record header
      $self->{recordinfo} = undef; # reference to the record contents
      $self->{line_number} = 0;   #  the line number in the FASTA file.
   }

# $obj_instance->_parse();

#This method parses a FASTA record from the file stream. 
#All the parsing errors for this record are recorded in 
#the logfile. If the record parses correctly, the method returns 1, else it 
#returns 0.
      
 
   sub _parse() {
      
      my $self = shift;
      my $last_line_length_lt_std_flag = 0;
      my $record_identifier = "<undefined>";
      my $preceding_header_flag = 0;
      my $first_data_line_length = 0;
      my @recarray;
      my $preceding_record_flag = 0;

      my $line_number = $self->{line_number};
      
      $self->{error_cnt} = 0; 
      my $db_handle = defined ( $self->{db_handle} ) ? 
         $self->{db_handle} : undef;
     
      # check for FASTA header of next record
      my $header = $self->{rec_header};
      
      #when in the beginning and end of the filestream.
      if(!defined ($header)) {
	  my $newline  = <$db_handle>;
          if((defined($newline)) && ($newline ne "")) {
             $line_number++;
	     $header = $newline;
	  }
      }
      #parsing the header
      if ( ( isValidFASTAheader($header) ) != 0 ) {
         # set up the variables for parsing a new record
         $last_line_length_lt_std_flag = 0;
         $preceding_header_flag = 1;
         $preceding_record_flag = 1;
         $self->{rec_header} = undef;

	 # if it's a valid FASTA header, then don't need to check again
         # extract the record IDENTIFIER
         $record_identifier = _headerToIdentifier($header); 
         push @recarray,$header;
      }
      else { #the header is not valid.
         if((defined $header) && (defined $line_number)) {
            $self->_errorHandler("Expected: record header " .
               "information in FASTA record header.  Got: \'$header\' ".
               "at line $line_number.", $DEBUG_LEVEL_3, $USR_ERR);
         }           
         $preceding_record_flag = 1;
         $self->{rec_header} = undef;
      }
      
      while ( ( defined ( $db_handle ) ) &&
              ( defined ( my $line = <$db_handle> )) &&
              (++$line_number)) {
         chomp $line;
         
         # check FASTA data
         if ( ( defined ( $record_identifier ) ) &&
              ( $record_identifier !~ /<undefined>/ ) &&
              ( ( isValidFASTAdata($line) ) != 0 ) ) {
	     
            push @recarray,$line;
	    if($preceding_header_flag == 1) {
	       $first_data_line_length = setValidFASTAlineLength($line);
            }
	    
            # check $last_line_length_lt_std_flag for an error on 
            # previous line
	    if(defined ($first_data_line_length)) {
               if ( $last_line_length_lt_std_flag == 1 ) {
                    $self->{error_cnt}++;
                    $self->_errorHandler("Expected: FASTA data ".
                       "definition lines should be ".
                       "$first_data_line_length bases (characters) ".
                       "across. Only the last line of a sequence ".
                       "data definition may be less than " .
                       "$first_data_line_length bases (characters) " .
                       "across, if applicable.  See line " . 
                       ($line_number - 1) . '.', $DEBUG_LEVEL_6, $USR_ERR);
               }
               $last_line_length_lt_std_flag = 0;
	    
               # check current line for over-length problem
               if ( length($line) > $first_data_line_length ) {
                  $self->{error_cnt}++;
                  $self->_errorHandler("Expected: FASTA data ".
                     "definition lines should be $first_data_line_length ".
                     "bases (characters) across. Only the last line of a ".
                     "sequence data definition may be less than " .
                     "$first_data_line_length bases (characters) ".
                     "across, if applicable.  See line " . $line_number . 
                     '.',
                     $DEBUG_LEVEL_6,$USR_ERR);
               }  
            
               #check current line for under-length problem; 
               #report only if not
               #the last line in the data definition
               elsif ( length($line) < $first_data_line_length ) {
                  $last_line_length_lt_std_flag = 1;
               }
	    } 
            $preceding_header_flag = 0;
	 }
         elsif($line =~ /^>/) { #the next header
            if ($preceding_record_flag == 1) { #its the next record.
	       $self->{rec_header} = $line;
               $self->_errorHandler(
                   "Assigned the header of the next record", 
                   $DEBUG_LEVEL_5, $SYS_ERR);
               last;
            }
         }
         # handle data error types
         else {
            $self->{error_cnt}++;

            # line has a separator token in it, so it may be header
            if ( $line =~ /$UNBOUND_FASTA_SEPARATOR/ ) {
               $self->_errorHandler("Expected: record header " .
                  "information in FASTA record header.  Got: \'$line\' ".
                  "at line $line_number.", $DEBUG_LEVEL_6, $USR_ERR);
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
                  "definition for record identifier ".
                  "\'$record_identifier\' Check sequence content at line ".
                  "$line_number for invalid bases ".
                  "(data type: invalid data).", $DEBUG_LEVEL_6, $USR_ERR);
            }
            else {
               $self->_errorHandler("Expected: FASTA record header " .
                  "followed by definition of sequence.  Invalid input at " .
                  "line $line_number.", $DEBUG_LEVEL_6, $USR_ERR);
            }
         }
      } # end while

      # check terminal case data definition
      if ( $preceding_header_flag == 1 ) {
         $self->_nullRecordHandler($line_number);
      }
      $self->{recordinfo} = \@recarray;
      $self->{line_number} = $line_number;
      return ( $self->{error_cnt} == 0 ) ? 1 : 0;
   }


# $obj_instance->_nullRecordHandler($);

#This method handles the case of a null or equivalently empty record
#encountered during parsing.  It logs the appropriate message to the 
#TIGR Foundation object.  The only argument is the line number.


   sub _nullRecordHandler($) {
      my $self = shift;
      my $line_number = shift;
 
      if ( ! defined ($line_number) ) {
         $line_number = "<unknown>";
      }

      $self->{error_cnt}++;
      if ( $self->{db_handle}->eof() == 1 ) {
           $self->_errorHandler("Expected: FASTA record header " .
              "followed by definition of sequence. " .
              "Got end of file after line " . 
              ($line_number) . ".", $DEBUG_LEVEL_5, $USR_ERR);
      }
      else {
         $self->_errorHandler("Expected: FASTA record header " .
            "followed by definition of sequence " .
            "Got FASTA header at line " . 
            ($line_number-1) . ".", $DEBUG_LEVEL_5, $USR_ERR);
      }
   }



# $message = $obj_instance->_errorHandler($message, $tf_level,
#   $internal_log_flag);

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

To use this module, load the C<TIGR::FASTAiterator> package via the
C<use> function.  Then, create a new instance of the object via the
C<new()> method, as shown below.  There are several invocations possible
for this method since all parameters to C<new()> are optional.

To access records from the C<TIGR::FASTAiterator> instance, the 
C<TIGR::FASTArecord> package must be loaded via the C<use> function.

An example script using this module follows.  The C<TIGR::Foundation>
module is included for completeness but does not have to be used.

   #!/usr/local/bin/perl -w

   # This script accepts FASTA files with the '-i' option
   # on the command line and validates every record in the file.
   # Parse errors for each record are collected to the 
   # '@errors_list' array and written to the .error file.
   # This program concatenates all of the correct records together to 
   # one output file specified with the '-o' option.
   # NOTE: The '-i' option must be specified before every input file.
   

   use strict;
   use TIGR::FASTAiterator;
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

      # Create a TIGR::FASTAiterator instance using TIGR::Foundation and
      # an error message list.
      my $fasta_iterator = new TIGR::FASTAiterator $tf_object, \@errors_list;

      if ( !(  defined ( $output_file ) &&
               open OUTFILE, ">$output_file" ) ) {
         $tf_object->bail("Cannot open output file for writing.");
      }

      foreach my $in_file ( @input_files ) {
         $fasta_iterator->open($in_file) or
         $tf_object->logLocal("Cannot open or read file $in_file", 2);

         if ( scalar(@errors_list) > 0 ) { # are there parse errors?
            while ( @errors_list ) { # get the messages from the list
               my $message = shift @errors_list; 
               print STDERR $message, "\n";
            }
         }
         #
         while ( $fasta_iterator->hasNext() ) {
            # print each record to OUTFILE
            my $record = $fasta_iterator->next();

            # print each record to OUTFILE
            if(defined $record) {
               print OUTFILE $record->toString();
	    }
         }
      }
   }

=cut

}

1;
