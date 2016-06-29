#!/usr/bin/env python3

"""
The initial use case for this script was needing a tool to append version numbers to all
entries within the FASTA file.  

For example, if you start with a header lines like this:

    >transcript1_comp4_path3
    >transcript1_comp4_path5

And you run this script like this:

    ./append_to_fasta_header.py -i foo.fna -o bar.fna -s '.1'

Then the new file will have the following headers:

    >transcript1_comp4_path3.1
    >transcript1_comp4_path5.1

Author: Joshua Orvis (jorvis AT gmail)

"""

import argparse
import gzip
import re
import sys

def main():
    parser = argparse.ArgumentParser( description='Append strings to the end of FASTA read headers')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read (gzip supported)' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    parser.add_argument('-s', '--suffix', type=str, required=True, help='String to be appended to each sequence header' )
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
                ofh.write(">{0}{1}{2}\n".format(m.group(1), args.suffix, m.group(2)))
            else:
                raise Exception("ERROR: Found a possible > entry with no identifier on line {0}".format(linenum))

        else:
            ofh.write(line)

    ofh.close()


if __name__ == '__main__':
    main()
