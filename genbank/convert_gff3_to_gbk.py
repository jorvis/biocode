#!/usr/bin/env python3

"""
Converts GFF3 representing gene models to Genbank flat-file format.

GFF3 specification:
http://www.sequenceontology.org/gff3.shtml

Genbank flat file specification:
https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html

"""

import argparse
import biocodegenbank
import biocodegff
import os
import sys

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input GFF3 file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to a Genbank flat file to be created' )
    args = parser.parse_args()

    (assemblies, features) = biocodegff.get_gff3_features( args.input_file )
    ofh = sys.stdout
    if args.output_file is not None:
        ofh = open(args.output_file, 'wt')

    for assembly_id in assemblies:
        assembly = assemblies[assembly_id]
        
        for gene in assemblies[assembly_id].genes():
            biocodegenbank.print_biogene( gene=gene, fh=ofh, on=assembly )


if __name__ == '__main__':
    main()







