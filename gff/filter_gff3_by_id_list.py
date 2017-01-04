#!/usr/bin/env python3

'''
This script is used if you want to filter a GFF3 file by IDs but also maintain any relationships 
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

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Filters the features of a GFF3 file by IDs while keeping related features.')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-l', '--list_file', type=str, required=True, help='A file with one identifier per line')
    args = parser.parse_args()

    molgraph = gff.parse_gff3_by_relationship(args.input_file)
    id_list  = parse_id_list( args.list_file )

    fout = open(args.output_file, mode='wt', encoding='utf-8')

    for mol_id in molgraph:
        mol = molgraph[mol_id]

        for (feat_id, feat) in sorted(mol.items(), key=lambda mol: int(mol[1]['fmin'])):
            ## used to filter whether this one is exported
            buffer = list()
            keep = False
            
            buffer.append( feat['cols'] )
            if feat_id in id_list:
                keep = True
                id_list[feat_id] = True

            for child in feat['children']:
                buffer.append( child['cols'] )
                if child['id'] is not None and child['id'] in id_list:
                    keep = True
                    id_list[ child['id'] ] = True

            if keep:
                for cols in buffer:
                    row = '\t'.join(cols)
                    fout.write( row )
                    fout.write( '\n' )



def parse_id_list(file):
    ids = dict()

    for line in open(file):
        ## value initialized to False, later set True if the script actually finds it
        ids[line.rstrip()] = False

    return ids

        

if __name__ == '__main__':
    main()






