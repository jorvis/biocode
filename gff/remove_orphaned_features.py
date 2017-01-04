#!/usr/bin/env python3

"""
This housekeeping script reads a GFF3 file and writes a new one, excluding any features
which aren't related to another feature.

The initial use-case here was a GFF file with gene features that had no mRNA children.

Warning: Because feature parentage can be anywhere in the file, it can be difficult to
do some things like this without holding the whole file into memory.  Instead, this makes
two passes over the file and just stores line indexes to keep, along with a few IDs to
track parentage.  A bit more processing time, much less memory needed.

Author: Joshua Orvis
"""

import argparse

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Removes orphaned features in a GFF3 file')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input GFF3 file' )
    parser.add_argument('-o', '--output', type=str, required=True, help='Output GFF3 file to write' )
    #parser.add_argument('-t', '--type', type=str, required=False, help='Type of features to remove' )
    args = parser.parse_args()

    # going to try saving memory by tracking line numbers instead of storing all of it
    #  true means keep the line, false means to omit it
    # doing tracking this way since it's technically legal for a feature to have no identifier at all.
    lines = list()
    parents = dict()
    current_line_num = -1

    infile = open(args.input)
    
    for line in infile:
        current_line_num += 1
        
        if line.startswith('#'):
            lines.append(True)
            continue
        
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) != 9:
            lines.append(True)
            continue

        id     = gff.column_9_value(cols[8], 'ID')
        parent = gff.column_9_value(cols[8], 'Parent')

        if parent is None:
            # this might be overwritten later
            lines.append(False)
            
            if id is not None:
                if parent not in parents:
                    parents[parent] = False
        else:
            lines.append(True)
            parents[parent] = True
                
    infile.seek(0)
    current_line_num = -1

    outfh = open(args.output, 'wt')
    
    for line in infile:
        current_line_num += 1

        if lines[current_line_num] == True:
            outfh.write(line)
        else:
            line = line.rstrip()
            cols = line.split("\t")

            if len(cols) == 9:
                id     = gff.column_9_value(cols[8], 'ID')

                if id is not None and id in parents and parents[id] == True:
                    outfh.write("{0}\n".format(line))
                else:
                    print("WARN: removing this line: {0}".format(line))
        
            

        

if __name__ == '__main__':
    main()







