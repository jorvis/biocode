#!/usr/bin/env python3

'''
This script reports some basic statistics about a GFF3 file.  It was
written with an initial focus on gene-structure containing content,
though can be expanded as desired.

The output is a tab-delimited file where the first column is the

Follow the GFF3 specification!

Author:  Joshua Orvis
'''

import argparse
import os
import biocodegff
import sys
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser( description='Checks the CDS features against a genome sequence to report/correct phase columns.')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    args = parser.parse_args()

    (assemblies, features) = biocodegff.get_gff3_features( args.input_file )

    ## output will either be a file or STDOUT
    ofh = sys.stdout
    if args.output_file is not None:
        ofh = open(args.output_file, 'wt')

    type_counts = defaultdict(int)
    type_lengths = defaultdict(int)
        
    for assembly_id in assemblies:
        type_counts['assembly'] += 1
        type_lengths['assembly'] += assemblies[assembly_id].length
        
        for gene in assemblies[assembly_id].genes():
            type_counts['gene'] += 1
            type_lengths['gene'] += gene.length
            
            for mRNA in gene.mRNAs():
                type_counts['mRNA'] += 1
                type_lengths['mRNA'] += mRNA.length

                for exon in mRNA.exons():
                    type_counts['exon'] += 1
                    type_lengths['exon'] += exon.length
                
                for CDS in mRNA.CDSs():
                    type_counts['CDS fragments'] += 1
                    type_lengths['CDS fragments'] += CDS.length
                    

    ofh.write("Assembly count\t{0}\n".format(type_counts['assembly']))
    ofh.write("Gene count\t{0}\n".format(type_counts['gene']))
    ofh.write("mRNA count\t{0}\n".format(type_counts['mRNA']))
    ofh.write("exon count\t{0}\n".format(type_counts['exon']))
    ofh.write("CDS count\t{0}\n".format(type_counts['CDS fragments']))

    ofh.write("Assembly length\t{0}\n".format(type_lengths['assembly']))
    ofh.write("Gene length\t{0}\n".format(type_lengths['gene']))
    ofh.write("mRNA length\t{0}\n".format(type_lengths['mRNA']))
    ofh.write("exon length\t{0}\n".format(type_lengths['exon']))
    ofh.write("CDS fragment length\t{0}\n".format(type_lengths['CDS fragments']))
        

if __name__ == '__main__':
    main()






