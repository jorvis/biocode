#!/usr/bin/env python3

"""

The use case here was that I generated sets of large SAM/BAM files where the contig
IDs had / characters within them, which wasn't legal for downstream loading into
JBrowse.  So I had to change the source FASTA files, but also needed to change the
downstream IDs as well.

This script allows for renaming these using a regex while keeping the other columns
of the SAM file intact.

"""

import argparse
import os


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    ofh = open(args.output_file, 'wt')
    
    # Your code goes here
    for line in open(args.input_file):
        if line.startswith('@SQ'):
            ofh.write(line.replace('/', '__'))
        else:
            cols = line.split("\t")
            if len(cols) >= 11:
                cols[2] = cols[2].replace('/', '__')
                cols[6] = cols[6].replace('/', '__')
                ofh.write("\t".join(cols))
            else:
                ofh.write(line)

if __name__ == '__main__':
    main()







