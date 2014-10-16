#!/usr/bin/env python3

"""
This script came about because I had a collection of individual Trinity RNA-Seq assemblies
and each had the same series of identifiers within them.  I wanted to create a merged version
of these files, but without the duplicate headers.

This script solves this by allowing you to prepend a string into the existing identifiers.

For example, if you start with a header lines like this:

    >transcript1_comp4_path3
    >transcript1_comp4_path5

And you run this script like this:

    ./prepend_to_fasta_header.py -i foo.fna -o bar.fna -p 'bmicroti_'

Then the new file will have the following headers:

    >bmicroti_transcript1_comp4_path3
    >bmicroti_transcript1_comp4_path5

Author: Joshua Orvis (jorvis AT gmail)

"""

import argparse
import gzip
import re
import sys

def main():
    parser = argparse.ArgumentParser( description='Prepend strings to the beginning of FASTA read headers')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read (gzip supported)' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    parser.add_argument('-p', '--prepend', type=str, required=False, help='String to be prepended to each sequence header' )
    args = parser.parse_args()

    linenum = 0

    if args.output_file is None:
        ofh = sys.stdout
    else:
        ofh = open( args.output_file, 'wt' )

    if args.input_file.endswith('.gz'):
        fh = gzip.open( args.input_file, 'rb')
        is_compressed = True
    else:
        fh = open(args.input_file, 'rU')
        is_compressed = False

    for line in fh:
        if is_compressed:
            line = line.decode()

        linenum += 1
            
        if line.startswith('>'):
            m = re.match('>(\S+)(.*)', line)
            if m:
                ofh.write(">{0}{1}{2}\n".format(args.prepend, m.group(1), m.group(2)))
            else:
                raise Exception("ERROR: Found a possible > entry with no identifier on line {0}".format(linenum))

        else:
            ofh.write(line)

    ofh.close()


if __name__ == '__main__':
    main()
