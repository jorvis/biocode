#!/usr/bin/env python3

"""
The GFF3 format supports embedding of the FASTA sequence after the feature table,
but this is often undesired when using that GFF3 file for things like track data
within IGV.js. 

This script removes the FASTA section of a GFF3 file and exports only the comments
and feature table.

"""

import argparse
import os

def main():
    parser = argparse.ArgumentParser( description='Remove FASTA from GFF3')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    in_fasta_section = False

    ofh = open(args.output_file, 'wt')

    for line in open(args.input_file):
        if line.startswith('##FASTA'):
            in_fasta_section = True

        if not in_fasta_section:
            ofh.write(line)

if __name__ == '__main__':
    main()







