#!/usr/bin/env python3

"""

This utility reads an annotated GFF3 file and writes FASTA of different feature types,
either the CDS or protein sequences (protein seqs are translated from the CDS).  In
either case, the genomic FASTA must be available, either embedded at the end of the
GFF3 or as a separate file (passed via -f).

The exported IDs will be the ID attributes of the mRNA features unless these have
locus_tag values in the 9th column.  Also, if the polypeptides associated with each
mRNA have product_name values, that will also be written to the FASTA header.

If the --check_ends option is passed, each CDS extracted is checked for canonical start
and stop sites.  When this happens, a warning is printed on STDERR like:

   WARN: Non-canonical start codon (AAA) in mRNA G678.mRNA.9856.1

Warnings are also printed if translating to protein and undefined codons are found:

   WARN: Encountered unknown codon during translation: NTA
   WARN: Non-canonical stop codon (CCA) in mRNA G678.mRNA.5532.1


Author:  Joshua Orvis (jorvis AT gmail)
"""

import argparse
import sys

from biocode import utils, gff


def main():
    parser = argparse.ArgumentParser( description='Extracts the protein or CDS seqeunces from a GFF3 file')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input GFF3 file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output FASTA file to be created' )
    parser.add_argument('-t', '--type', type=str, required=False, default='protein', choices=['protein', 'cds'], help='Type of features to export')
    parser.add_argument('-f', '--fasta', type=str, required=False, help='If the FASTA entries for the underlying assemblies is absent from the GFF3 document passed, you will need to specify this option' )
    parser.add_argument('-ft', '--feature_type', type=str, required=False, default='mRNA', choices=['mRNA', 'polypeptide'], help='IDs and coordinates will come from this feature type' )
    parser.add_argument('--check_ends', dest='check_ends', action='store_true')
    parser.add_argument('--check_internal_stops', dest='check_internal_stops', action='store_true')
    parser.set_defaults(check_ends=False, check_internal_stops=False)
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output_file is not None:
        fout = open(args.output_file, 'wt')

    # sanity option check
    if args.check_internal_stops == True and args.type == 'cds':
        raise Exception("Error:  Checking internal stops for CDS features not currently supported.")

    (assemblies, features) = gff.get_gff3_features(args.input_file)

    # only doing the standard codon table for now
    start_codons = ['ATG', 'GTG', 'TTG']
    stop_codons  = ['TAG', 'TAA', 'TGA']

    ## add sequence residues from external FASTA file if the user passed one
    if args.fasta is not None:
        utils.add_assembly_fasta(assemblies, args.fasta)
    
    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            if args.feature_type == 'mRNA':
                feats = gene.mRNAs()
            elif args.feature_type == 'polypeptide':
                feats = gene.polypeptides()
            
            for feat in feats:
                ## initial values of id and header to export (can be overridden by available annotation)
                export_id = feat.id
                export_header = None

                ## Add the gene product name if there is one
                if args.feature_type == 'mRNA':
                    for polypeptide in feat.polypeptides():
                        if polypeptide.annotation is not None:
                            if polypeptide.annotation.product_name is not None:
                                export_header = polypeptide.annotation.product_name
                                break

                    coding_seq = feat.get_CDS_residues(for_translation=True)
                    if feat.locus_tag is not None:
                        export_id = feat.locus_tag
                        
                elif args.feature_type == 'polypeptide':
                    export_header = feat.annotation.product_name
                    coding_seq = feat.parent.get_CDS_residues(for_translation=True)
                    if feat.parent.locus_tag is not None:
                        export_id = feat.parent.locus_tag

                if len(coding_seq) > 0:
                    fout.write(">{0}".format(export_id))
                    if export_header is not None:
                        fout.write(" {0}\n".format(export_header))
                    else:
                        fout.write("\n")

                    if args.check_ends == True:
                        # check the starting codon
                        start_codon = coding_seq[0:3].upper()
                        if start_codon not in start_codons:
                            sys.stderr.write("WARN: Non-canonical start codon ({0}) in mRNA {1}\n".format(start_codon, feat.id))

                        stop_codon = coding_seq[-3:].upper()
                        if stop_codon not in stop_codons:
                            sys.stderr.write("WARN: Non-canonical stop codon ({0}) in mRNA {1}\n".format(stop_codon, feat.id))                        

                    if args.type == 'cds':
                        fout.write("{0}\n".format(utils.wrapped_fasta(coding_seq)))
                    else:
                        translated_seq = utils.translate(coding_seq)

                        if args.check_internal_stops == True:
                            internal_stop_count = translated_seq[:-1].count('*')
                            if internal_stop_count > 0:
                                sys.stderr.write("Found {0} internal stops in mRNA {1}\n".format(internal_stop_count, feat.id))

                        fout.write("{0}\n".format(utils.wrapped_fasta(translated_seq)))
                else:
                    print("WARNING: Skipped feature {0} because it had no associated CDS features".format(export_id), file=sys.stderr)

if __name__ == '__main__':
    main()







