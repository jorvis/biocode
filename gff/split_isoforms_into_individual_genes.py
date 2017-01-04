#!/usr/bin/env python3

'''
There are some scripts and utilities out there which can't handle genes with
multiple mRNA children.  This script splits each of these (and children) and
instead creates new gene features.

The ID of the newly-generated gene feature duplicates the source one but with
the suffix "_N" added, where N increases from 2..X for each isoform present for
a given gene.

Follow the GFF3 specification!

Author:  Joshua Orvis
'''

import argparse
import sys

from biocode import gff, things


def main():
    parser = argparse.ArgumentParser( description='Checks for genes with multiple mRNA children and creates new genes for each.')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_file)

    ## output will either be a file or STDOUT
    ofh = sys.stdout
    if args.output_file is not None:
        ofh = open(args.output_file, 'wt')

    ofh.write("##gff-version 3\n")

    for assembly_id in assemblies:
        current_assembly = assemblies[assembly_id]
        
        for gene in assemblies[assembly_id].genes():
            rnas_found = 0
            mRNAs = gene.mRNAs()
            
            for mRNA in mRNAs:
                mRNA_loc = mRNA.location_on(current_assembly)
                rnas_found += 1

                if rnas_found > 1:
                    gene.remove_mRNA(mRNA)
                    
                    print("INFO: splitting mRNA off gene {0}".format(gene.id))
                    new_gene = things.Gene(id="{0}_{1}".format(gene.id, rnas_found))
                    new_gene.locate_on(target=current_assembly, fmin=mRNA_loc.fmin, fmax=mRNA_loc.fmax, strand=mRNA_loc.strand)
                    new_gene.add_RNA(mRNA)
                    new_gene.print_as(fh=ofh, format='gff3')

            if len(mRNAs) > 1:
                gene_loc = gene.location_on(current_assembly)
                mRNA_loc = mRNAs[0].location_on(current_assembly)
                gene_loc.fmin = mRNA_loc.fmin
                gene_loc.fmax = mRNA_loc.fmax
                gene_loc.strand = mRNA_loc.strand

            gene.print_as(fh=ofh, format='gff3')


if __name__ == '__main__':
    main()






