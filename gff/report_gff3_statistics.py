#!/usr/bin/env python3

'''
This script reports some basic statistics about a GFF3 file.  It was
written with an initial focus on gene-structure containing content,
though can be expanded as desired.

The output is a tab-delimited file where the first column is the description
of a statistic and the second is the value.

Warning:  If you parse this output you'll need to skip blank lines and any
which begin with the # symbol.

Follow the GFF3 specification!

Author:  Joshua Orvis
'''

import argparse
import sys
from collections import defaultdict

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Checks the CDS features against a genome sequence to report/correct phase columns.')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_file)

    ## output will either be a file or STDOUT
    ofh = sys.stdout
    if args.output_file is not None:
        ofh = open(args.output_file, 'wt')

    type_counts = defaultdict(int)
    type_lengths = defaultdict(int)
    assembly_lengths_found = False

    # key is number of exons, value is number of mRNAs with that many
    CDS_profile = defaultdict(int)
        
    for assembly_id in assemblies:
        type_counts['assembly'] += 1

        if assemblies[assembly_id].length is not None:
            type_lengths['assembly'] += assemblies[assembly_id].length
            assembly_lengths_found = True
        
        for gene in assemblies[assembly_id].genes():
            type_counts['gene'] += 1
            type_lengths['gene'] += gene.length
            
            for mRNA in gene.mRNAs():
                type_counts['mRNA'] += 1
                type_lengths['mRNA'] += mRNA.length
                CDS_profile[mRNA.CDS_count()] += 1

                for exon in mRNA.exons():
                    type_counts['exon'] += 1
                    type_lengths['exon'] += exon.length
                
                for CDS in mRNA.CDSs():
                    type_counts['CDS fragments'] += 1
                    type_lengths['CDS fragments'] += CDS.length
                    

    ofh.write("Assembly count\t{0}\n".format(type_counts['assembly']))
    if assembly_lengths_found:
        ofh.write("Assembly length\t{0}\n".format(type_lengths['assembly']))
    else:
        ofh.write("Assembly length\tN/A (no FASTA data in GFF?)\n")

    gene_length_mean = type_lengths['gene'] / type_counts['gene']
    mRNA_length_mean = type_lengths['mRNA'] / type_counts['mRNA']
    exon_length_mean = type_lengths['exon'] / type_counts['exon']
    CDS_length_mean = type_lengths['CDS fragments'] / type_counts['CDS fragments']

    mRNAs_per_gene_mean = type_counts['mRNA'] / type_counts['gene']
    exons_per_mRNA_mean = type_counts['exon'] / type_counts['mRNA']
    CDS_per_mRNA_mean = type_counts['CDS fragments'] / type_counts['mRNA']
    
    ofh.write("\nGene count\t{0}\n".format(type_counts['gene']))
    ofh.write("Gene length (mean)\t{0:.1f}\n".format(gene_length_mean))
    ofh.write("Gene length (sum)\t{0}\n".format(type_lengths['gene']))
    
    
    ofh.write("\nmRNA count\t{0}\n".format(type_counts['mRNA']))
    ofh.write("mRNA length (mean)\t{0:.1f}\n".format(mRNA_length_mean))
    ofh.write("mRNA length (sum)\t{0}\n".format(type_lengths['mRNA']))
    ofh.write("mRNAs per gene (mean)\t{:.1f}\n".format(mRNAs_per_gene_mean) )
    
    ofh.write("\nexon count\t{0}\n".format(type_counts['exon']))
    ofh.write("exon length (mean)\t{0:.1f}\n".format(exon_length_mean))
    ofh.write("exon length (sum)\t{0}\n".format(type_lengths['exon']))
    ofh.write("exons per mRNA (mean)\t{:.1f}\n".format(exons_per_mRNA_mean) )

    ofh.write("\nCDS count\t{0}\n".format(type_counts['CDS fragments']))
    ofh.write("CDS length (mean)\t{0:.1f}\n".format(CDS_length_mean))
    ofh.write("CDS fragment length (sum)\t{0}\n".format(type_lengths['CDS fragments']))
    ofh.write("CDS per mRNA (mean)\t{:.1f}\n".format(CDS_per_mRNA_mean) )
    
    ofh.write("\n# CDS fragment composition profile: count<tab>percentage\n")
    for cds_count in sorted(CDS_profile):
        perc = (CDS_profile[cds_count] / type_counts['mRNA']) * 100
        ofh.write("mRNAs with {0} CDS\t{1}\t{2:.3}\n".format(cds_count, CDS_profile[cds_count], perc) )

if __name__ == '__main__':
    main()






