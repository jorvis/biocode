#!/usr/bin/env python3

"""

Put some general, high-level documentation here

"""

import argparse
import biocodegff
import biocodeutils
import os
import sys


def main():
    parser = argparse.ArgumentParser( description='Extracts the protein or CDS seqeunces from a GFF3 file')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input GFF3 file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output FASTA file to be created' )
    parser.add_argument('-t', '--type', type=str, required=False, default='protein', choices=['protein', 'cds'], help='Type of features to export')
    parser.add_argument('-f', '--fasta', type=str, required=False, help='If the FASTA entries for the underlying assemblies is absent from the GFF3 document passed, you will need to specify this option' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output_file is not None:
        fout = open(args.output_file, 'wt')

    (assemblies, features) = biocodegff.get_gff3_features(args.input_file)

    ## add sequence residues from external FASTA file if the user passed one
    if args.fasta is not None:
        biocodeutils.add_assembly_fasta(assemblies, args.fasta)
    
    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            for mRNA in gene.mRNAs():

                ## initial values of id and header to export (can be overridden by available annotation)
                export_id = mRNA.id
                export_header = None

                if mRNA.locus_tag is not None:
                    export_id = mRNA.locus_tag

                ## Add the gene product name if there is one
                for polypeptide in mRNA.polypeptides():
                    if polypeptide.annotation is not None:
                        if polypeptide.annotation.product_name is not None:
                            export_header = polypeptide.annotation.product_name
                            break
                
                fout.write(">{0}".format(export_id))
                if export_header is not None:
                    fout.write(" {0}\n".format(export_header))
                else:
                    fout.write("\n")
                
                coding_seq = mRNA.get_CDS_residues()

                if args.type == 'cds':
                    fout.write("{0}\n".format(biocodeutils.wrapped_fasta(coding_seq)))
                else:
                    translated_seq = biocodeutils.translate(coding_seq)
                    fout.write("{0}\n".format(biocodeutils.wrapped_fasta(translated_seq)))

if __name__ == '__main__':
    main()







