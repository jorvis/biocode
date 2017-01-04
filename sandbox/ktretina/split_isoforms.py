#!/usr/bin/env python3

"""

split_isoforms.py - splits all GFF3 mRNA isoforms into their own gene models (adds _# to old gene ID to get new gene ID)

"""

import argparse

from biocode import gff, things


def main():
    parser = argparse.ArgumentParser( description='Splits all GFF3 mRNA isoforms into their own gene models')

    ## Get the variables
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Input GFF3 file' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output GFF3 file' )
    args = parser.parse_args()
    ofh = open(args.output_file, 'wt')

    print("INFO: Parsing GFF3 features\n")
    (assemblies, ref_features) = gff.get_gff3_features(args.input_file)

    print("INFO: Finding genes with isoforms and splitting them\n")
    ofh.write("##gff-version 3\n")
    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            # only changing the gene features with isoforms
            if len(gene.mRNAs()) > 1:
                counter = 1
                for mRNA in gene.mRNAs():
                    new_gene_id = str(gene.id) + "_" + str(counter)
                    counter += 1
                    mRNA_loc = mRNA.location() 
                    print("Splitting " + gene.id)
                    # create a new gene model, correcting the gene coords to the mRNA coords
                    new_gene = things.Gene(id = new_gene_id)
                    new_gene.locate_on( target=assemblies[assembly_id], fmin=mRNA_loc.fmin, fmax=mRNA_loc.fmax, strand=mRNA_loc.strand )
                    mRNA.parent.id = new_gene_id
                    #Now add the mRNA to the gene model
                    new_gene.add_mRNA(mRNA)
                    # print out the new gene model
                    new_gene.print_as(fh=ofh, source='IGS', format='gff3')
            else:
                gene.print_as(fh=ofh, source='IGS', format='gff3')


if __name__ == '__main__':
    main()


