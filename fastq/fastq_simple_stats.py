#!/usr/bin/env python3.2

import argparse
import os

## CONTACT: jorvis@gmail.com

def main():
    parser = argparse.ArgumentParser( description='Provides simple quantitative statistics for a given FASTQ file')

    ## output file to be written
    parser.add_argument('input_files', metavar='N', type=str, nargs='+', help='Path to one or more input files, separated by spaces' )
    args = parser.parse_args()

    ## values that will be reported
    entry_count = 0
    total_bases = 0
    
    line_number = 0

    for input_file in args.input_files:
        for line in open(input_file, 'r'):
            line.rstrip()
            line_number += 1

            ## every 4th line is the start of a sequence entry (check still)
            if line_number % 4 == 1:
                if line.startswith('@'):
                    entry_count += 1
                else:
                    raise Exception("Error, expected every 4th line to start with @ and this one didn't: {0}".format(line) )
            elif line_number % 4 == 2:
                total_bases += len(line)


    avg_entry_length = total_bases / entry_count
    
    print("Total sequence entries: {0}".format(entry_count))
    print("Total bases: {0}".format(total_bases))
    print("Avg sequence length: {0:.1f}".format(avg_entry_length))

if __name__ == '__main__':
    main()







