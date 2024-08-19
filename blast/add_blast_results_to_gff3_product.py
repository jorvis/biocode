#!/usr/bin/env python3

"""

Imagine you're doing a quick and dirty annotation of a new GFF3 file, such as the
output of a gene predictor, and you want to use the top BLAST match to assign
gene products.

This script does that, assuming tabular BLAST output.

Also assumes that the IDs in the BLAST output correspond to the mRNA IDs in the GFF

By default, output is assumed to have been created with this option:

    -outfmt '6 qseqid sseqid evalue bitscore pident stitle' 

Yes, this is messy, but it can be an educational exercise for students.
"""

import argparse
import os
import re

from biocode import utils, gff

def main():
    parser = argparse.ArgumentParser( description='Updates GFF3 files with gene products from BLAST')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input GFF3 file to be read' )
    parser.add_argument('-b', '--blast', type=str, required=True, help='Path to an input BLAST tab file to be read' )
    parser.add_argument('-s', '--source', type=str, required=True, help='Becomes the source (2nd) column of the GFF output' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output GFF3 file to be created' )
    
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_file)

    for line in open(args.blast):
        if line.startswith(">"):
            continue
        
        line = line.rstrip()
        cols = line.split("\t")

        ## TODO: Make this configurable
        qid = cols[0]
        product = cols[5]

        ## SWISS-PROT matches look like these.  Let's process it:
        #   RecName: Full=Threonine synthase; Short=TS [Escherichia coli K-12]
        #   RecName: Full=Trp operon repressor [Escherichia coli O139:H28 str. E24377A]
        products = product.split(';')
        m = re.match("RecName: Full=(.+)", products[0])
        if m:
            product = m.group(1)

            # if it has a bracketed genus/species, remove it.
            m = re.match("(.+)\[.*", product)
            if m:
                product = m.group(1)

        mrna = features[qid]
        polypeptide = mrna.children['polypeptide'][0]
        polypeptide.annotation.product_name = product

    ofh = open(args.output_file, 'wt')
    ofh.write("##gff-version 3\n\n")

    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            gene.print_as(fh=ofh, source=args.source, format='gff3')
    
    ofh.close()
        
            

if __name__ == '__main__':
    main()







