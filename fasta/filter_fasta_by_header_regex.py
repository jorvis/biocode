#!/usr/bin/env python3

'''
Description:

This script reads a FASTA file and matches each header against a regular
expression passed by the user, only exporting those records which matched
it.

The use case for creating this was we had a mixed FASTA file of transcripts
from several different RNA-Seq assemblers.  They each had different naming
conventions, like:

   >Locus_1_Transcript_1/803056_Confidence_0.000_Length_397
   >TR2|c0_g1_i1 len=241 path=[462:0-30 463:31-54 464:55-240] [-1, 462, 463, 464, -2]

And we wanted to separate these out.  It's also conceivable you might want to use
this to filter FASTA files by annotated properties, for example.

Input: 

- A (multi-)FASTA file.
- The > symbol is NOT included in what is matched by the regex, so don't include
  it in your pattern.

Output:

A FASTA file containing the matching records.

'''

import argparse
import os
import re
import sys

def main():
    parser = argparse.ArgumentParser( description='Split multi-FASTA file into separate protein and nucleotide files')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    parser.add_argument('-r', '--regex', type=str, required=True, help='Regular expression to match against each header' )
    args = parser.parse_args()

    if args.output_file is None:
        ofh = sys.stdout
    else:
        ofh = open( args.output_file, 'wt' )

    keep = False
    
    for line in open(args.input_file):
        if line.startswith('>'):
            m = re.search(args.regex, line[1:])
            if m:
                keep = True
            else:
                keep = False

        if keep:
            ofh.write(line)

if __name__ == '__main__':
    main()







