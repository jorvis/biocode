#!/usr/bin/env python3

"""
This housekeeping script reads a GFF3 file and writes a new one, correcting the Parent
attribute of any exon features to ensure that they point to their mRNA/tRNA/rRNA.

The initial use-case here was a GFF file dumped from WebApollo which had this issue.

WARNING: this script is extremely dependent on feature order, such that the RNA feature
appears before the child exons.  This shouldn't be an issue, but it could attach features
to the wrong parent if they were out of order.

Author: Joshua Orvis
"""

import argparse

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Updates exon Parent attributes to point at the correct RNA feature')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input GFF3 file' )
    parser.add_argument('-o', '--output', type=str, required=True, help='Output GFF3 file to write' )
    args = parser.parse_args()

    infile = open(args.input)
    ofh = open(args.output, 'wt')

    last_rna_id = None
    
    for line in infile:
        
        if line.startswith('#'):
            ofh.write(line)
            continue
        
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) != 9:
            ofh.write("{0}\n".format(line) )
            continue

        id     = gff.column_9_value(cols[8], 'ID')
        parent = gff.column_9_value(cols[8], 'Parent')

        if cols[2].endswith('RNA'):
            last_rna_id = id
            ofh.write("{0}\n".format(line) )

        elif cols[2] == 'exon':
            if parent != last_rna_id:
                print("INFO: correcting unexpected parentage for feature ({0}) type {2}.  Expected ({1})".format(id, last_rna_id, cols[2]) )
                cols[8] = gff.set_column_9_value(cols[8], 'Parent', last_rna_id)
                ofh.write("{0}\n".format("\t".join(cols)) )
            else:
                ofh.write("{0}\n".format(line) )
        else:
            ofh.write("{0}\n".format(line) )
            

        

if __name__ == '__main__':
    main()







