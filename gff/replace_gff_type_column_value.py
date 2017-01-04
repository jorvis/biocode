#!/usr/bin/env python3

"""
This simple script allows you to replace the type (3rd) column in a GFF with another value.  For
example, WebApollo writes 'transcript' features rows, and many of our scripts assume these will
be 'mRNA' instead:

./replace_gff_type_column_value.py -i in.gff -it transcript -o out.gff -ot mRNA

Author: Joshua Orvis
"""

import argparse


def main():
    parser = argparse.ArgumentParser( description='Adds gene features for RNAs which lack them')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input GFF3 file' )
    parser.add_argument('-it', '--input_type', type=str, required=True, help='Type of feature you want to replace' )
    parser.add_argument('-o', '--output', type=str, required=True, help='Output GFF3 file to write' )
    parser.add_argument('-ot', '--output_type', type=str, required=True, help='New feature type value' )
    args = parser.parse_args()

    infile = open(args.input)
    ofh = open(args.output, 'wt')

    for line in infile:
        
        if line.startswith('#'):
            ofh.write(line)
            continue
        
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) != 9:
            ofh.write("{0}\n".format(line) )
            continue

        if cols[2] == args.input_type:
            cols[2] = args.output_type
            ofh.write("{0}\n".format("\t".join(cols)) )
        else:
            ofh.write("{0}\n".format(line) )


if __name__ == '__main__':
    main()







