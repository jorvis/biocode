#!/usr/bin/env python3

'''
Description:

This script reads a FASTA file and matches each header against a regular
expression passed by the user, stripping the rest of the header sequence 
after (or including) the regex portion, depending on passed options.

The use case for creating this was that I had a FASTA file with headers like this:

    >HiC_scaffold_361:1-14571
    >HiC_scaffold_352:1-14780
    >HiC_scaffold_350:1-14861
    >HiC_scaffold_1977:1-1518

But I wanted only the identifier portions, like this:

    >HiC_scaffold_361
    >HiC_scaffold_352
    >HiC_scaffold_350
    >HiC_scaffold_1977

So I passed --regex '\:' and --include=no

Input: 

- A (multi-)FASTA file.
- The > symbol is NOT included in what is matched by the regex, so don't include
  it in your pattern.
- If the pattern isn't found at all, the entire header will still be exported.

Output:

A FASTA file with headers appropriately modified.

'''

import argparse
import os
import re
import sys

def main():
    parser = argparse.ArgumentParser( description='Modified the headers of a FASTA file based on user options')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    parser.add_argument('-r', '--regex', type=str, required=True, help='Regular expression to match against each header' )
    parser.add_argument('-n', '--include', type=str, required=True, help='Include the regex portion in the output? Values are yes or no' )
    args = parser.parse_args()

    if args.include not in ['yes', 'no']:
        raise Exception("ERROR: The only valid values for --include are 'yes' or 'no'")

    if args.output_file is None:
        ofh = sys.stdout
    else:
        ofh = open( args.output_file, 'wt' )

    for line in open(args.input_file):
        if line.startswith('>'):
            line = line.rstrip()
            
            m = re.match("(.*)({0})(.*)".format(args.regex), line[1:])
            if m:
                if args.include == 'yes':
                    line = ">{0}{1}\n".format(m.group(1), m.group(2))
                elif args.include == 'no':
                    line = ">{0}\n".format(m.group(1))

        ofh.write(line)

if __name__ == '__main__':
    main()







