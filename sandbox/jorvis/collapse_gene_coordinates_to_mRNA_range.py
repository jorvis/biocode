#!/usr/bin/env python3

'''
This script was needed when some tools (such as WebApollo) were generating gene models where the
gene coordinates extended past the boundaries of the mRNA child feature, and it was required
to trim them to a matching range.

INPUT GFF3

The encoding convention for gene models is:

Supercontig_3.10        EVM     gene    23723   28269   .       +       .       ID=evm.TU.Supercontig_3.10.9;Name=EVM%20prediction%20Supercontig_3.10.9
Supercontig_3.10        EVM     mRNA    23723   28269   .       +       .       ID=evm.model.Supercontig_3.10.9;Parent=evm.TU.Supercontig_3.10.9
Supercontig_3.10        EVM     exon    23723   24696   .       +       .       ID=evm.model.Supercontig_3.10.9.exon1;Parent=evm.model.Supercontig_3.10.9
Supercontig_3.10        EVM     CDS     23723   24696   .       +       0       ID=cds.evm.model.Supercontig_3.10.9;Parent=evm.model.Supercontig_3.10.9
Supercontig_3.10        EVM     exon    24756   28269   .       +       .       ID=evm.model.Supercontig_3.10.9.exon2;Parent=evm.model.Supercontig_3.10.9
Supercontig_3.10        EVM     CDS     24756   28269   .       +       1       ID=cds.evm.model.Supercontig_3.10.9;Parent=evm.model.Supercontig_3.10.9


Follow the GFF3 specification!

Author:  Joshua Orvis
'''

import argparse

import gff


def main():
    parser = argparse.ArgumentParser( description='Shortens gene feature coordinates to their longest child mRNA')

    ## output file to be written
    parser.add_argument('-i', '--input_gff3', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-o', '--output_gff3', type=str, required=True, help='Path to GFF3 output file to be created')
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_gff3)
    gff_out = open(args.output_gff3, 'wt')

    gff_out.write("##gff-version 3\n")

    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            gene_loc = gene.location()

            # loop through the mRNAs and store the outer boundaries of those found
            min_coord = None
            max_coord = None

            mRNAs = gene.mRNAs()

            if len(mRNAs) >= 1:
                for mRNA in mRNAs:
                    mRNA_loc = mRNA.location()

                    if min_coord is None or mRNA_loc.fmin < min_coord:
                        min_coord = mRNA_loc.fmin

                    if max_coord is None or mRNA_loc.fmax > max_coord:
                        max_coord = mRNA_loc.fmax

                if min_coord != gene_loc.fmin or max_coord != gene_loc.fmax:
                    print("DEBUG: Changed gene {0} from {1}-{2} to {3}-{4}".format(gene.id, gene_loc.fmin, gene_loc.fmax, min_coord, max_coord))
                    gene_loc.fmin = min_coord
                    gene_loc.fmax = max_coord
                
            gene.print_as(fh=gff_out, source='IGS', format='gff3')



if __name__ == '__main__':
    main()






