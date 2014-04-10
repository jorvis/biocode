#!/usr/bin/env python3

'''
If you've generated GFF3 which represents gene models and then later performed repeatmasking
this script allows you to filter your models based on overlap of those predicted repeats.

INPUT GFF3

The encoding convention for gene models is:

Supercontig_3.10        EVM     gene    23723   28269   .       +       .       ID=evm.TU.Supercontig_3.10.9;Name=EVM%20prediction%20Supercontig_3.10.9
Supercontig_3.10        EVM     mRNA    23723   28269   .       +       .       ID=evm.model.Supercontig_3.10.9;Parent=evm.TU.Supercontig_3.10.9
Supercontig_3.10        EVM     exon    23723   24696   .       +       .       ID=evm.model.Supercontig_3.10.9.exon1;Parent=evm.model.Supercontig_3.10.9
Supercontig_3.10        EVM     CDS     23723   24696   .       +       0       ID=cds.evm.model.Supercontig_3.10.9;Parent=evm.model.Supercontig_3.10.9
Supercontig_3.10        EVM     exon    24756   28269   .       +       .       ID=evm.model.Supercontig_3.10.9.exon2;Parent=evm.model.Supercontig_3.10.9
Supercontig_3.10        EVM     CDS     24756   28269   .       +       1       ID=cds.evm.model.Supercontig_3.10.9;Parent=evm.model.Supercontig_3.10.9

Any embedded FASTA sequence within the GFF3 is ignored in favor of the FASTA passed via --masked_fasta

INPUT FASTA

The FASTA needs to have molecule IDs corresponding to column 1 in the GFF3 file, and
contained sequence masked with N characters.  (soft-masking not currently supported.)


Follow the GFF3 specification!

Author:  Joshua Orvis
'''

import argparse
import os
import biocodegff
import biocodeutils


def main():
    parser = argparse.ArgumentParser( description='Removes gene models whose sequence has been masked.')

    ## output file to be written
    parser.add_argument('-i', '--input_gff3', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-m', '--masked_fasta', type=str, required=True, help='FASTA with sequence masked with N characters')
    parser.add_argument('-o', '--output_gff3', type=str, required=False, help='Path to GFF3 output file to be created')
    args = parser.parse_args()

    (assemblies, features) = biocodegff.get_gff3_features( args.input_gff3 )
    biocodeutils.add_assembly_fasta(assemblies, args.masked_fasta)

    gff_out = open(args.output_gff3, 'wt')
    gene_count = 0
    kept_count = 0
        
    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            keep = True
            gene_count += 1
            
            for mRNA in gene.mRNAs():
                
                coding_seq = mRNA.get_CDS_residues()
                n_count = coding_seq.count('N')
                perc_repeat = (n_count / len(coding_seq)) * 100

                if perc_repeat >= 10:
                    keep = False

            if keep == True:
                kept_count += 1
                gene.print_as(fh=gff_out, source='IGS', format='gff3')


    print("INFO: {0} genes kept out of {1} ({2:.1f}%)".format(kept_count, gene_count, ((kept_count/gene_count) * 100)))

if __name__ == '__main__':
    main()






