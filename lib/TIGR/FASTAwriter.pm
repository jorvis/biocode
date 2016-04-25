# $Id: FASTAwriter.pm,v 1.1 2004/04/28 15:03:43 aphillip Exp $

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

package TIGR::FASTAwriter;
{

=head1 NAME

TIGR::FASTAwriter - TIGR::FASTAwriter class for writing TIGR::FASTArecord
objects to a file

=head1 SYNOPSIS

  use TIGR::FASTAwriter;
  my $obj_instance = new TIGR::FASTAwriter ($tigr_foundation_obj,
                                            $output_file_name);

=head1 DESCRIPTION

This module provides an object definition for a TIGR::FASTAwriter. 
The TIGR::FASTAwriter object accepts TIGR::FASTArecord objects for
printing to an output file.

=cut

   BEGIN {
      require 5.006_00;
   }

   use strict;
   use IO::File;
   use TIGR::FASTArecord;



   ## internal variables and identifiers

   our $REVISION = (qw$Revision: 1.1 $)[-1];
   our $VERSION = '1.0';
   our $VERSION_STRING = "$VERSION (Build $REVISION)";
   our @DEPEND = ();
   
   my $SYS_ERR = 0;           # this flag specifies non-user related error
   my $USR_ERR = 1;           # this flag specifies user related error
   
   #debugging scheme
   #
   #   Debugging via the TIGR Foundation uses increasing log levels based on
   #   nesting. 'MAIN' starts at level 1. Every nest increments the level by 1.
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
   sub write($);
   sub _errorHandler($$$);

   ## implementation

=over

=item $obj_instance = new TIGR::FASTAwriter ($foundation_object,
   $error_array_ref, $output_file);

This method returns a new instance of a TIGR::FASTAwriter object. It takes
three optional parameters: a TIGR::Foundation object (C<$foundation_object>),
a reference to an array for logging user error messages (C<$error_array_ref>),
and an output file name, C<$output_file>, as parameters. A new object instance 
is returned on success and successful opening of a specified output 
file.
If the file supplied cannot be opened, this method returns undefined.
This method also returns undefined if the parameters supplied are invalid.
Writing errors are written to the array at C<$error_array_ref> and the 
log file.

=cut


   sub new(;$$$) {
      
      my $pkg = shift;
      my @method_args = @_;

      my $error_condition = 0;
      my $self = {};
      bless $self, $pkg;
      
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
            $self->{db_file_name} = "" ;
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
      elsif ((defined ( $self->{db_file_name} )) &&
           (! defined ( $self->open($self->{db_file_name}, "w") ) )) {
         # the error message is logged via the open() routine
         $self = undef;
      }
      
      return ( $error_condition == 0 ) ? $self : undef;
     
   }


=item $result = $obj_instance->open($file_name, $flag);

This method opens a FASTA file for writing or appending.  The file, 
F<$file_name>, is opened using the C<open()> flags specified by C<$flag>.
Supported flags include 'w' and 'a'.  On success, this method returns 1.
The default C<open()> method is 'w', or truncated open.  If the file cannot
be opened, this method returns undefined.

=cut


   sub open($;$) {
      my $self = shift;
      my $db_file_name = shift;
      my $open_flags = shift;

      my $error_condition = 0;

      if ( ( ! defined ($open_flags) ) ||
           ( ( $open_flags !~ /^w$/i ) &&
             ( $open_flags !~ /^a$/i ) ) ) {
         $open_flags = "w";
      }

      # close a previously open file
      if ( defined ($self->{db_handle}) ) {
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
         if ( ! defined ($return_val) ) {
            $return_val = undef;
            $self->_errorHandler(
               "Error closing FASTA file: $self->{db_file_name}", 
                $DEBUG_LEVEL_4, $USR_ERR);
         }
      }
      $self->{db_file_name} = undef;
      $self->{db_handle} = undef;
      return $return_val;
   }
 

=item $result = $obj_instance->write($fasta_obj);

This method takes a TIGR::FASTArecord object, C<$fasta_obj>, and writes it
to the file specified in C<new()> or C<open()>.  On success, this method 
returns true (1).  On error, this method returns false (undefined) and logs
an error message.

=cut


   sub write($) {
      my $self = shift;
      my $fasta_obj = shift;
      my $return_val = 1;
      
      if ( ( defined ($fasta_obj) ) &&
           ( ( ref($fasta_obj) ) =~ /fastarecord/i ) &&
           ( defined ($self->{db_handle}) ) ) {
	 
         if ( $self->{db_handle}->print($fasta_obj->toString()) ) {
	    $return_val = 1;
         }
         else {
	    
            $return_val = undef;
            $self->_errorHandler(
               "Error printing to FASTA file: $self->{db_file_name}", 
                $DEBUG_LEVEL_3, $USR_ERR);
         }
      }
      else {
         $return_val = undef;
         $self->_errorHandler(
               "Invalid method of initialization for " .
                     "TIGR::FASTAwriter", $DEBUG_LEVEL_3, $USR_ERR);
      }
      return $return_val;
   }

# $message = $obj_instance->_errorHandler($message, $tf_level,
#    $internal_log_flag);

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


=back

=head1 USAGE

To use this module, load the C<TIGR::FASTArecord> and C<TIGR::FASTAwriter>
modules with the C<use> function. Then, create a new instance of the object 
via the C<new()> method, as shown below. There are several invocations 
possible for this method since all parameters to C<new()> are optional.
An example script using this module follows. The C<TIGR::Foundation>
module is included for completeness but does not have to be used.

   #!/usr/local/bin/perl -w

   # This example uses the TIGR::FASTAwriter object to write 
   # a simple TIGR::FASTArecord object to a file specified with
   # the '-o' option to this script.
   # Writing errors are collected to the '@errors_list' array.
 
   use strict;
   use TIGR::Foundation;
   use TIGR::FASTArecord;
   use TIGR::FASTAwriter;

   MAIN:
   {
      my $tf_object = new TIGR::Foundation;
      my @errors_list = ();
      my $output_file = undef;

      my $getopts_result = undef;

      $getopts_result = $tf_object->TIGR_GetOptions( "o=s" => \$output_file );
      
      if ( $getopts_result != 1 ) {
         $tf_object->bail("Invalid command line option.");
      }

      if ( ! defined ( $output_file ) ) {
         $tf_object->bail("Must specify an output file with the '-o' option");
      }

      my $header = "ORF00001";
      my $data = "ATGC";

      my $fasta_record = new TIGR::FASTArecord $header, $data;
      if ( ! defined ( $fasta_record ) ) {
         $tf_object->bail("Cannot create TIGR::FASTArecord object");
      }
      
      # Create a TIGR::FASTAwriter instance using TIGR::Foundation and
      # an error message list.

      my $fasta_writer = new TIGR::FASTAwriter $tf_object, \@errors_list;

      $fasta_writer->open($output_file) or 
         $tf_object->logLocal("Cannot open output file $output_file", 
                               $DEBUG_LEVEL_1);

      if ( scalar(@errors_list) > 0 ) { # are there parse errors?
         while ( @errors_list ) { # get the messages from the list
            my $message = shift @errors_list; 
            print STDERR $message, "\n";
         }
      }

      $fasta_writer->write($fasta_record ) or 
         $tf_object->logLocal("Cannot write FASTA record to $output_file", 
                               $DEBUG_LEVEL_1);
   }

=cut

}

1;
