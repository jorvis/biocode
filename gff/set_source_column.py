#!/usr/bin/env python3

"""
Sometimes you just need to overwrite the source (2nd) column in a GFF3 file.  This script does
that and leaves all other comment or FASTA lines in place.
"""

import argparse
import os
import sys


def main():
    parser = argparse.ArgumentParser( description='Overwrite source (2nd) column in GFF3 file')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    parser.add_argument('-s', '--source', type=str, required=True, help='The new source value to be set' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    ofh = sys.stdout
    if args.output_file is not None:
        ofh = open(args.output_file, 'wt')

    for line in open(args.input_file):
        cols = line.split("\t")

        if len(cols) == 9:
            cols[1] = args.source
            ofh.write("\t".join(cols))
        else:
            ofh.write(line)


if __name__ == '__main__':
    main()







