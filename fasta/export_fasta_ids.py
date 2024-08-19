#!/usr/bin/env python3

"""
Takes as input a FASTA file and outputs a list of the IDs in that file.

Yes, this could be done with a shell script or command-line Fu, but this is a Python script.
"""

import argparse
import re
import sys

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    parser.add_argument('infile', nargs='?', help='Path to an input file to be read' )
    parser.add_argument('outfile', nargs='?', help='Path to an output file to be created' )
    args = parser.parse_args()

    if not args.infile:
        print("No input file specified.")
        sys.exit(1)

    if args.outfile:
        ofh = open(args.outfile, 'w')
    else:
        ofh = sys.stdout

    with open(args.infile, 'r') as fh:
        for line in fh:
            if line[0] == ">":
                m = re.search(r'^>(\S+)', line)
                if m:
                    ofh.write(m.group(1) + '\n')

if __name__ == '__main__':
    main()







