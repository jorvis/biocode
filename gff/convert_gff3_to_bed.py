#!/usr/bin/env python3

"""
This script can be used to transform a GFF3 file into a BED file.  

LIMITATIONS:

- The only currently-tested input GFF3 consists of coding gene models.
- Its encoding is meant to include prokaryotes and eukaryotes, though it was
  written/tested first for eukaryotes.

The options:

--input_file:  Must follow the GFF3 specification
--output_file: Output BED file to be created

"""

import argparse

from biocode import utils, gff, things, bed


def main():
    parser = argparse.ArgumentParser( description='Create a BED file from GFF3')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to output file to be created' )
    
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_file)
    
    ofh = open(args.output_file, 'wt')
    bed.print_bed_from_assemblies(assemblies=assemblies, ofh=ofh)

if __name__ == '__main__':
    main()







