#!/usr/bin/env python3.2

'''
This script is used if you want to filter a FASTA file by IDs but also maintain any relationships 
defined within the file.  That is, if you have a list of mRNA identifiers and want to keep them, 
but also wanted to keep any associated exons, CDS, etc.

Because common GFF3 usage doesn't guarantee that features are defined before they are referenced
as a parent, this creates a graph of objects for all features on a molecule first, then iterates
through to filter them.  It doesn't matter on what feature level you pass IDs, the entire graph
of vertical relatives for that feature will be retained.

What happens to everything else?  Obviously all the other feature graphs are ommitted, but currently
so are any associated comments, sequence data, etc.

Follow the GFF3 specification!

Author:  Joshua Orvis
'''

import argparse
import os
import biocodegff
from operator import itemgetter



def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-l', '--list_file', type=str, required=True, help='A file with one identifier per line')
    args = parser.parse_args()

    molgraph = biocodegff.parse_gff3_by_relationship( args.input_file )

    fout = open(args.output_file, mode='wt', encoding='utf-8')

    for mol_id in molgraph:
        mol = molgraph[mol_id]
        for feat in sorted(mol.items(), key=lambda mol: mol['fmin']):
            row = '\t'.join(feat['cols'])
            fout.write( row )
            fout.write( '\n' )

if __name__ == '__main__':
    main()







