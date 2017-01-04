#!/usr/bin/env python3

"""
This housekeeping script reads a GFF3 file and writes a new one, adding a 'gene'
row for any RNA feature which doesn't have one.  The coordinates of the RNA will
be copied.

The initial use-case here was a GFF file dumped from WebApollo which had this issue.

In this particular use case, the orphan mRNAs have ID attributes but no Parent
though this is corrected.

INPUT EXAMPLE:
###
ChromosomeII_BmicrotiR1	IGS	mRNA	1467897	1468187	.	+	.	Name=ChromosomeII_BmicrotiR1:1467871-1468187;ID=101D714C468A44840D49A6FAAD27AFE5
ChromosomeII_BmicrotiR1	IGS	exon	1467897	1468187	.	+	.	Name=DE1443B2DABA5DEDBDEBE79EB433EEB8;Parent=101D714C468A44840D49A6FAAD27AFE5;ID=DE1443B2DABA5DEDBDEBE79EB433EEB8
ChromosomeII_BmicrotiR1	IGS	CDS	1467897	1468187	.	+	0	Name=101D714C468A44840D49A6FAAD27AFE5-CDS;Parent=101D714C468A44840D49A6FAAD27AFE5;ID=101D714C468A44840D49A6FAAD27AFE5-CDS

Author: Joshua Orvis
"""

import argparse

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Adds gene features for RNAs which lack them')

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

        id     = gff.column_9_value(cols[8], 'ID')
        parent = gff.column_9_value(cols[8], 'Parent')

        if cols[2].endswith('RNA') and parent is None:
            gene_cols = list(cols)
            gene_cols[2] = 'gene'
            gene_cols[8] = gff.set_column_9_value(gene_cols[8], 'ID', "{0}.gene".format(id))
            ofh.write("{0}\n".format("\t".join(gene_cols)) )

            cols[8] = gff.set_column_9_value(cols[8], 'Parent', "{0}.gene".format(id))
            ofh.write("{0}\n".format("\t".join(cols)) )
        else:
            ofh.write("{0}\n".format(line) )


if __name__ == '__main__':
    main()







