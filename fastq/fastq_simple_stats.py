#!/usr/bin/env python3

"""
Reads one or more FASTQ files and reports some basic statistics.  Example:

$ ./fastq_simple_stats.py my.fastq

File: my.fastq
Total sequence entries: 277971801
Total bases: 28353123702
Avg sequence length: 102.0
GC percentage: 44.5

Contact: jorvis@gmail.com
"""


import argparse
import codecs
import gzip
import os
import sys


def main():
    parser = argparse.ArgumentParser( description='Provides simple quantitative statistics for a given FASTQ file')

    ## output file to be written
    parser.add_argument('input_files', metavar='N', type=str, nargs='+', help='Path to one or more input files, separated by spaces' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Optional path to an output file to be created, else prints on STDOUT' )
    args = parser.parse_args()

    ## open the output file
    fout = None
    if args.output_file is None:
        fout = codecs.getwriter('utf8')(sys.stdout.buffer)
    else:
        fout = open(args.output_file, "w")

    ## values that will be reported
    entry_count = 0
    total_bases = 0
    gc_count    = 0
    
    line_number = 0

    for input_file in args.input_files:
        if input_file.endswith('.gz'):
            fh = gzip.open( input_file, 'rb')
            is_compressed = True
        else:
            fh = open( input_file, 'rU' )
            is_compressed = False
        
        for line in fh:
            if is_compressed:
                line = line.decode()
            
            line = line.rstrip()
            line_number += 1

            ## every 4th line is the start of a sequence entry (check still)
            if line_number % 4 == 1:
                if line.startswith('@'):
                    entry_count += 1
                else:
                    raise Exception("Error, expected every 4th line to start with @ and this one didn't: {0}".format(line) )
            elif line_number % 4 == 2:
                total_bases += len(line)
                gc_count += line.count('G')
                gc_count += line.count('C')


    avg_entry_length = total_bases / entry_count

    if len(args.input_files) > 1:
        fout.write("File: multiple\n")
    else:
        fout.write("File: {0}\n".format(args.input_files[0]))

    fout.write("Total sequence entries: {0}\n".format(entry_count))
    fout.write("Total bases: {0}\n".format(total_bases))
    fout.write("Avg sequence length: {0:.1f}\n".format(avg_entry_length))
    fout.write("GC percentage: {0:.1f}\n".format( (gc_count / total_bases) * 100))
    

if __name__ == '__main__':
    main()







