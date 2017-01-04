#!/usr/bin/env python3

'''
Description:

FASTA is a relatively easy format to get right, but sometimes accidents happen.
This script checks the more common expectations of a FASTA file but, given that
there's no real 'official' specification that I'm aware of, isn't perfect.

Here's what it does check for:

- ERROR: Entry with 0-length sequence
- ERROR: No > symbols embedded in sequence residues (suggests merged records)
- ERROR: Nonzero-length identifier after > symbol until first whitespace
- ERROR: Multiple records within an individual file with same identifier
- ERROR: Incorrect record count (if --expected_record_count is passed)
- ERROR: Make sure the file ends with a blank line (else commonly used tools like cat will cause errors)
- WARNING: Residue lines > 60bp
- WARNING: Comment lines with # (most parsers will fail to handle this)
- WARNING: Homopolymers of length > --homopolymer_limit which are not Ns
- WARNING: Internal stops (*) detected in protein FASTA (if --check_internal_stops is passed)

Input: 

- One or more FASTA files passed positionally
- A list of FASTA files passed via --input_list
- Any combination of these


Output:


This will all be printed to STDOUT or a file if you pass the -o option.

'''

import argparse
import re
import sys

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Performs selected validation of a FASTA file')

    RESIDUE_LINE_LENGTH_LIMIT = 60

    ## output file to be written
    parser.add_argument('fasta_files', metavar='N', type=str, nargs='*', help='Pass one or more FASTA files')
    parser.add_argument('-o', '--output_file', type=str, required=False, default=None, help='Optional path to an output file to be created' )
    parser.add_argument('-l', '--input_list', type=str, required=False, default=None, help='Optional path to a list of list of input files' )
    parser.add_argument('-erc', '--expected_record_count', type=int, required=False, default=None, help='Optional count of records expected in the input.  An exception is raised if this is not matched.' )
    parser.add_argument('-hl', '--homopolymer_limit', type=int, required=False, default=None, help='Issues a warning for any sequences with a homopolymer of length > this' )
    parser.add_argument('-cis', '--check_internal_stops', dest='check_internal_stops', action='store_true', help='Meant when processing protein files, checks for internal * characters')
    parser.set_defaults(check_internal_stops=False)
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output_file is not None:
        fout = open(args.output_file, 'wt')

    ## gather input
    input_files = args.fasta_files

    if args.input_list is not None:
        input_files.extend(utils.read_list_file(args.input_list))

    if len(input_files) == 0:
        raise Exception("ERROR: No input files defined")
        
    total_records = 0
    total_residues = 0
    error_count = 0
    warning_count = 0

    for ifile in input_files:
        fout.write("INFO: validating file {0}\n".format(ifile))
        line_number = 0
        ids_found = list()
        long_line_count = 0
        current_seq_length = None
        current_seq_id = None
        current_seq = ''
        current_homopolymer_base = None
        current_homopolymer_length = 0
        
        for line in open(ifile):
            line = line.rstrip()
            line_number += 1
            
            if line.startswith('>'):
                total_records += 1

                # found a new sequence entry - make sure the last one (if there was one) had residues
                if current_seq_length == 0:
                    fout.write("ERROR: Entry without residues found just before line {0} of {1}\n".format(line_number, ifile))
                    error_count += 1

                if args.check_internal_stops == True:
                    if '*' in current_seq[0:-1]:
                        fout.write("WARNING: Internal stops found within {0} of {1}\n".format(current_seq_id, ifile))
                        warning_count += 1
                
                current_seq = ''
                current_seq_length = 0
                m = re.match(">(\S+)", line)
                if m:
                    current_seq_id = m.group(1)
                    if current_seq_id in ids_found:
                        fout.write("ERROR: Duplicate ID ({2}) found on line {0} of {1}\n".format(line_number, ifile, current_seq_id))
                        error_count += 1
                    else:
                        ids_found.append(current_seq_id)
                else:
                    fout.write("ERROR: Record without ID on line {0} of {1}\n".format(line_number, ifile))
                    error_count += 1
                
            elif line.startswith('#'):
                # warn about a comment line
                fout.write("WARNING: Comment detected on line {0} of {1}\n".format(line_number, ifile))
                warning_count += 1

            else:
                # residue line
                total_residues += len(line)
                current_seq_length += len(line)
                current_seq += line

                if args.homopolymer_limit is not None:
                    for base in line:
                        if base == current_homopolymer_base:
                            current_homopolymer_length += 1
                        else:
                            if current_homopolymer_length > args.homopolymer_limit and current_homopolymer_base != 'N':
                                fout.write("WARNING: Sequence ID {0} in file {1} contains a homopolymer run ({2}) of length {3}\n".format(current_seq_id, ifile, current_homopolymer_base, current_homopolymer_length))
                                warning_count += 1

                            current_homopolymer_base = base
                            current_homopolymer_length = 1

                if '>' in line:
                    fout.write("ERROR: > character embedded in sequence residues on line {0} of {1}\n".format(line_number, ifile))
                    error_count += 1

                # not practical to print warnings for each line here, so we'll just do it once per file
                if len(line) > RESIDUE_LINE_LENGTH_LIMIT:
                    long_line_count += 1

        if current_seq_length == 0:
            fout.write("ERROR: Entry without residues found on line {0} of {1}\n".format(line_number, ifile))
            error_count += 1
        
        if long_line_count > 0:
            fout.write("WARNING: {0} residue line(s) were detected longer than {1} in file {2}\n".format(long_line_count, RESIDUE_LINE_LENGTH_LIMIT, ifile))
            warning_count += long_line_count

    if args.expected_record_count is not None and args.expected_record_count != total_records:
        fout.write("ERROR: Expected record count:{0} does not match what was found:{1}\n".format(args.expected_record_count, total_records))
        error_count += 1

    if error_count == 0:
        if warning_count == 0:
            fout.write("INFO: All files in input set appear to be valid\n")
        else:
            fout.write("INFO: All files in input set appear to be valid, but with {0} warning(s)\n".format(warning_count))
    else:
        fout.write("ERROR: total errors found in all files from input set: {0}, warning(s): {1}\n".format(error_count, warning_count))
                  



if __name__ == '__main__':
    main()







