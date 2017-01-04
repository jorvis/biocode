#!/usr/bin/env python3

"""
The default output from glimmerHMM (when using the -g option) produces something LIKE
GFF but has only mRNA and CDS features.  This script transforms this into the canonical
genes described in the specification.


Example input:
  Supercontig_3.1 GlimmerHMM      mRNA    15493   15926   .       +       .       ID=Supercontig_3.1.path1.gene7;Name=Supercontig_3.1.path1.gene7
  Supercontig_3.1 GlimmerHMM      CDS     15493   15562   .       +       0       ID=Supercontig_3.1.cds7.1;Parent=Supercontig_3.1.path1.gene7;Name=Supercontig_3.1.path1.gene7;Note=initial-exon
  Supercontig_3.1 GlimmerHMM      CDS     15853   15926   .       +       2       ID=Supercontig_3.1.cds7.2;Parent=Supercontig_3.1.path1.gene7;Name=Supercontig_3.1.path1.gene7;Note=final-exon

Example output:
  Supercontig_3.1 GlimmerHMM      gene    15493   15926   .       +       .       ID=Supercontig_3.1.path1.gene7;Name=Supercontig_3.1.path1.gene7
  Supercontig_3.1 GlimmerHMM      mRNA    15493   15926   .       +       .       ID=Supercontig_3.1.path1.gene7.mRNA;Name=Supercontig_3.1.path1.gene7.mRNA
  Supercontig_3.1 GlimmerHMM      CDS     15493   15562   .       +       0       ID=Supercontig_3.1.path1.gene7.cds;Parent=Supercontig_3.1.path1.gene7.mRNA;Name=Supercontig_3.1.path1.gene7;Note=initial-exon
  Supercontig_3.1 GlimmerHMM      exon     15493   15562   .       +       0       ID=Supercontig_3.1.path1.gene7.exon.1;Parent=Supercontig_3.1.path1.gene7.mRNA;Name=Supercontig_3.1.path1.gene7;Note=initial-exon
  Supercontig_3.1 GlimmerHMM      CDS     15853   15926   .       +       2       ID=Supercontig_3.1.path1.gene7.cds;Parent=Supercontig_3.1.path1.gene7.mRNA;Name=Supercontig_3.1.path1.gene7;Note=final-exon
  Supercontig_3.1 GlimmerHMM      exon     15853   15926   .       +       2       ID=Supercontig_3.1.path1.gene7.exon.2;Parent=Supercontig_3.1.path1.gene7.mRNA;Name=Supercontig_3.1.path1.gene7;Note=final-exon



Author: Joshua Orvis
Contact: jorvis AT gmail
"""

import argparse
from collections import defaultdict

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Converts glimmerHMM GFF output to GFF3')

    # output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to parse' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    
    args = parser.parse_args()

    fout = open(args.output_file, 'w')

    current_gene = None
    current_mRNA = None

    next_exon_num = defaultdict(int)

    for line in open(args.input_file, 'r'):
        if line.startswith('#'):
            fout.write(line)
            continue

        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) != 9:
            continue

        mol_id = cols[0]
        feat_type = cols[2]
        feat_fmin = int(cols[3]) - 1
        feat_fmax = int(cols[4])

        id = gff.column_9_value(cols[8], 'ID')
        parent = gff.column_9_value(cols[8], 'Parent')

        if feat_type == 'mRNA':
            gene_cols = list(cols)
            gene_cols[2] = 'gene'

            cols[8] = gff.set_column_9_value(cols[8], 'ID', "{0}.mRNA".format(id))
            cols[8] = gff.set_column_9_value(cols[8], 'Name', "{0}.mRNA".format(id))
            cols[8] = gff.order_column_9(cols[8])
            
            # print the gene and mRNA
            fout.write( "{0}\n".format("\t".join(gene_cols)) )
            fout.write( "{0}\n".format("\t".join(cols)) )
            
        elif feat_type == 'CDS':
            exon_cols = list(cols)

            cols[8] = gff.set_column_9_value(cols[8], 'ID', "{0}.cds".format(parent))
            cols[8] = gff.set_column_9_value(cols[8], 'Name', "{0}.cds".format(parent))
            cols[8] = gff.set_column_9_value(cols[8], 'Parent', "{0}.mRNA".format(parent))
            cols[8] = gff.order_column_9(cols[8])

            exon_id = "{0}.exon.{1}".format(parent, next_exon_num[parent] )
            next_exon_num[parent] += 1
            
            exon_cols[2] = 'exon'
            exon_cols[7] = '.'
            exon_cols[8] = gff.set_column_9_value(exon_cols[8], 'ID', exon_id)
            exon_cols[8] = gff.set_column_9_value(exon_cols[8], 'Name', exon_id)
            exon_cols[8] = gff.set_column_9_value(exon_cols[8], 'Parent', "{0}.mRNA".format(parent))
            exon_cols[8] = gff.order_column_9(exon_cols[8])

            fout.write( "{0}\n".format("\t".join(exon_cols)) )
            fout.write( "{0}\n".format("\t".join(cols)) )


if __name__ == '__main__':
    main()





















