#!/usr/bin/env python3

"""
The official GFF3 specification is great about describing the document structurall and provides
some examples of different gene graphs.  What it does NOT do almost at all is provide conventions
for storing functional annotation on genes.  In the absence of this, I decided what I thought made
the most sense for our annotations at TIGR, JCVI and later at IGS.  A main component of this was
the decision to put annotation of protein-coding genes on the polypeptide features.

Recently, NCBI published specifications for doing genome submission more directly using GFF3 as
input, including conventions for storing functional annotation.  NCBI being a much more influential
source, this will supercede what I decided to do within biocode.  For now, this script converts
biocode-conventional GFF3 to NCBI's expectations for prokaryotes and eukaryotes.

GFF3 specification:
http://www.sequenceontology.org/gff3.shtml

NCBI GFF3 usage conventions:
https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/

MAJOR assumptions:
- Functional annotation is attached to a polypeptide feature in the source GFF3

"""

import argparse
import os
import sys

from biocode import utils, gff, ncbigff

def main():
    parser = argparse.ArgumentParser( description='Converts biocode GFF3 into NCBI-spec GFF3')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input GFF3 file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to a Genbank GFF3 file to be created.' )
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_file)
    ofh = sys.stdout
    if args.output_file is not None:
        ofh = open(args.output_file, 'wt')

    for assembly_id in assemblies:
        assembly = assemblies[assembly_id]

        for gene in assemblies[assembly_id].genes():
            ncbigff.print_biogene(gene=gene, fh=ofh, on=assembly, source='INCOMPLETE')

        

if __name__ == '__main__':
    main()







