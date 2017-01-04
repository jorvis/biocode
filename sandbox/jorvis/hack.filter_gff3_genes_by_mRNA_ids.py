#!/usr/bin/env python3

"""
There is a filter_gff3_by_id_list.py script, but it doesn't currently have support for
features which share IDs (like segmented CDS) and I need something quick.  The use
case here is that you have a list of mRNA.ids and this well export them along with their
child and parent features.

WARNING: This script very much assumes a feature-type-sorted order where the gene appears
first, then mRNA, then any child features of that.

"""

import argparse
import sys

from biocode import gff


def main():
    parser = argparse.ArgumentParser('Filter the genes of a GFF3 file by mRNA child IDs')

    ## output file to be written
    parser.add_argument('-i', '--input_gff3', type=str, required=True, help='GFF3 file of source molecules' )
    parser.add_argument('-l', '--id_list', type=str, required=True, help='List file of mRNA IDs to keep' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Optional output file path (else STDOUT)' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output_file is not None:
        fout = open(args.output_file, 'wt')

    ids_to_keep = list()

    for line in open(args.id_list):
        line = line.rstrip()
        if len(line) > 2:
            ids_to_keep.append(line)
        
    fout.write("##gff-version 3\n")

    current_gene_lines = list()
    current_gene_id = None
    keep = False
        
    for line in open(args.input_gff3):
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) != 9:
            continue

        # grab the ID and Parent columns if any
        id = gff.column_9_value(cols[8], 'ID')
        parent = gff.column_9_value(cols[8], 'Parent')

        type = cols[2]

        if type == 'gene':
            # purge the current gene, if any
            if len(current_gene_lines) > 1:
                for li in current_gene_lines:
                    fout.write("{0}\n".format(li) )

            # reset
            current_gene_lines = list()
            current_gene_lines.append( line )
            current_gene_id = id

        else:
            if type == 'mRNA':
                if id in ids_to_keep:
                    keep = True
                else:
                    keep = False

            if keep == True:
                current_gene_lines.append(line)


        





if __name__ == '__main__':
    main()







