#!/usr/bin/env python3

"""
A memory-efficient script for subsampling single or paired-end FASTQ files
randomly (rather than taking the first N, or some other non-random approach)
while keeping mate pairs together in order (if running with paired-end options)

In terms of compute complexity this is at best O(n) and at worst O(3n).  
The input file (or left file in the case of paired-end data) must be read once to
get a read count (unless you pass the --input_read_count option), then a
second time to do the subsampling, but no more than one record (4 lines) is ever
held in memory at a time.

The best other implementation of this I could find was in C, but required the
storage of an array of index positions.  I wrote this script to subsample 200
million reads out of a fileset of 2.5 billion, and this took too much memory.

Referenced C version:
  http://homes.cs.washington.edu/~dcjones/fastq-tools/fastq-sample.html

Python implementation which has the same issue:
  http://pythonforbiologists.com/index.php/randomly-sampling-reads-from-a-fastq-file/

Overview of steps:

1.  Parse --single or --left files to get a read count (skipped if --input_read_count provided)
2.  Pass through files once, randomly selecting reads based on the probability
    needed for each to get the requested --output_read_count
3.  A partial third pass through the files *might* be needed if not quite enough
    were pulled randomly in the first selection pass.  This will usually be < 5%,
    if needed at all.

So, in practice, if you know the read count beforehand, the complexity will most often
be around O(1.05N)

Output:

The files created are defined using the --output_base option.  If processing a single
file, this will just append .fastq to the value passed.  If processing paired-end
files, the following will be created:

   {output_base}.R1.fastq
   {output_base}.R2.fastq

This is true regardless of how many input (comma-delimited) files are passed, only
one output file will be written per direction.
  
Assumptions:

Input reads can be in raw FASTQ format or compressed with a .gz extension.  If
passing paired-end data, all pairwise entries must be in the same order within
their files (this is very common.)

Author:
Joshua Orvis
"""

import argparse
import gzip
import os
import random


def main():
    parser = argparse.ArgumentParser( description='Randomized subsampling of single or paired-end FASTQ files')

    ## output file to be written
    parser.add_argument('-l', '--left', type=str, required=False, help='Comma-separated list of left reads in an paired set' )
    parser.add_argument('-r', '--right', type=str, required=False, help='Comma-separated list of right reads in an paired set' )
    parser.add_argument('-s', '--single', type=str, required=False, help='Use this if you have only one file of unpaired reads' )
    parser.add_argument('-irc', '--input_read_count', type=int, required=False, help='If passed, skips a step of reading over the file(s) once just to get the read count' )
    parser.add_argument('-orc', '--output_read_count', type=int, required=True, help='Subsample size - number of reads to export')
    parser.add_argument('-ob', '--output_base', type=str, required=True, help='Base name of output files to be created')
    args = parser.parse_args()

    ######################
    ## Validate, then gather files and get read counts
    ######################
    
    check_input_options(args)
    input_read_count = args.input_read_count
    left_files = list()
    right_files = list()

    if args.single is not None:
        left_files = args.single.split(',')
    elif args.left is not None:
        left_files = args.left.split(',')
        right_files = args.right.split(',')

    if input_read_count is None:
        input_read_count = get_read_count(left_files)
        print("INFO: parsed input read count is: {0}".format(input_read_count))
    else:
        print("INFO: user-passed input read count is: {0}".format(input_read_count))

    #######################
    ## Do a first random selection pass
    #######################

    # Calculate the probability of selecting any individual read, given
    #  the proportion we want to select.
    individual_read_prop = args.output_read_count / input_read_count
    print("INFO: The probability of any individual read being selected is: {0}".format(individual_read_prop))

    output_reads_written = 0
    
    if args.single is not None:
        ofh = open("{0}.fastq".format(args.output_base), 'wt')
        ofh_right = None
    else:
        ofh = open("{0}.R1.fastq".format(args.output_base), 'wt')
        ofh_right = open("{0}.R2.fastq".format(args.output_base), 'wt')

    lindex = 0
    for lpath in left_files:
        if lpath.endswith('.gz'):
            ifh = gzip.open(lpath, 'rb')
            left_is_compressed = True

            if args.right is not None:
                ifh_right = gzip.open(right_files[lindex], 'rb')
        else:
            ifh = open(lpath, 'rU')
            left_is_compressed = False
            if args.right is not None:
                ifh_right = open(right_files[lindex], 'rU')

        for line1 in ifh:
            line2 = ifh.readline()
            line3 = ifh.readline()
            line4 = ifh.readline()

            if args.right is not None:
                rline1 = ifh_right.readline()
                rline2 = ifh_right.readline()
                rline3 = ifh_right.readline()
                rline4 = ifh_right.readline()

            if random.random() < individual_read_prop:
                if left_is_compressed:
                    ofh.write(line1.decode())
                    ofh.write(line2.decode())
                    ofh.write(line3.decode())
                    ofh.write(line4.decode())

                    if args.right is not None:
                        ofh_right.write(rline1.decode())
                        ofh_right.write(rline2.decode())
                        ofh_right.write(rline3.decode())
                        ofh_right.write(rline4.decode())
                else:
                    ofh.write(line1)
                    ofh.write(line2)
                    ofh.write(line3)
                    ofh.write(line4)

                    if args.right is not None:
                        ofh_right.write(rline1)
                        ofh_right.write(rline2)
                        ofh_right.write(rline3)
                        ofh_right.write(rline4)

                output_reads_written += 1

                if output_reads_written == args.output_read_count:
                    break
        lindex += 1

        if output_reads_written == args.output_read_count:
            # This one shouldn't really happen, since it means the record count was reached
            #  before all files were opened, but we still need to add it.
            break

    print("INFO: {0} reads written to output file after first pass".format(output_reads_written))
            

def check_input_options(args):
    # user must pass either singleton option, right AND left options, but not both
    if args.single is not None:
        if args.left is not None or args.right is not None:
            raise Exception("ERROR: Cannot use --left or --right options if --single is passed")
    elif args.left is not None:
        if args.right is None:
            raise Exception("ERROR: Must pass --right option if --left option is used")

        left_files = args.left.split(',')
        right_files = args.right.split(',')
        
        if len(left_files) != len(right_files):
            raise Exception("ERROR: Count of left files ({0}) must be the same as right files ({1})".format(
                len(left_files), len(right_files)))
    else:
        raise Exception("ERROR: You must pass either --single or --left AND --right to define input")
        
def get_read_count(files):
    line_count = 0

    for path in files:
        if path.endswith('.gz'):
            fh = gzip.open(path, 'rb')
        else:
            fh = open(path, 'rU')

        for line in fh:
            line_count += 1

    return int(line_count / 4)
    

if __name__ == '__main__':
    main()







