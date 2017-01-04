#!/usr/bin/env python3

"""
This housekeeping script reads a GFF3 file and reverses all CDS coordinates where stop < start.

The initial use-case here was a GFF file dumped from WebApollo which had this issue.

INPUT EXAMPLE:
unitig_0|quiver IGS gene    2622416 2625146 .   +   .   Name=598BBB62C7AAA2833D58F914778A9B48;ID=598BBB62C7AAA2833D58F914778A9B48
unitig_0|quiver IGS mRNA    2622416 2625146 .   +   .   Name=unitig_0|quiver:2623455-2623886;Parent=598BBB62C7AAA2833D58F914778A9B48;ID=CC226B225CA7A622FE45C9F664899D9E
unitig_0|quiver IGS exon    2624443 2625146 .   +   .   Name=C2E998B335DEAC334DE78741E6E8D536;Parent=CC226B225CA7A622FE45C9F664899D9E;ID=C2E998B335DEAC334DE78741E6E8D536
unitig_0|quiver IGS exon    2622416 2623082 .   +   .   Name=D45C00C52F9C5556594F632A93270FF6;Parent=CC226B225CA7A622FE45C9F664899D9E;ID=D45C00C52F9C5556594F632A93270FF6
unitig_0|quiver IGS exon    2623100 2624381 .   +   .   Name=5E2F4885E4EF1C9D4B2CA97206797637;Parent=CC226B225CA7A622FE45C9F664899D9E;ID=5E2F4885E4EF1C9D4B2CA97206797637
unitig_0|quiver IGS CDS 2624620 2622890 .   +   0   Name=CC226B225CA7A622FE45C9F664899D9E-CDS;Parent=CC226B225CA7A622FE45C9F664899D9E;ID=CC226B225CA7A622FE45C9F664899D9E-CDS;Note=Manually set translation start

Author: Kyle Tretina
"""

import argparse

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Reverses CDS coodinates where stop < start')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input GFF3 file' )
    parser.add_argument('-o', '--output', type=str, required=True, help='Output GFF3 file to write' )
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

        if cols[2] == 'CDS' and int(cols[4]) < int(cols[3]):
            temp = cols[3]
            cols[3] = cols[4]
            cols[4] = temp
            id = gff.column_9_value(cols[8], 'ID')
            print("CDS reversed: {0}".format(id))
            ofh.write("{0}\n".format("\t".join(cols)) )
        else:
            ofh.write("{0}\n".format(line) )


if __name__ == '__main__':
    main()







