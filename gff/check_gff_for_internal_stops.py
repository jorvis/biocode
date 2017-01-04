#!/usr/bin/env python3

'''
This script reads a GFF3 file and FASTA file (or FASTA embedded in the GFF) and
checks the translation of all CDS for internal stops.  The output is like:

    Total mRNAs found:4091
    mRNAs with embedded stops: 895

The script checks and removes terminal stop codons, so these will not be reported
if stops are found at the very end of the reading frame.

If you're trying to debug and want to see some of the CDS which are reported to contain
internal stops, then just use the --print_n_with_stops option and it will print that
many sequences which were found to contain stops.  

Follow the GFF3 specification!

Author:  Joshua Orvis
'''

import argparse

from biocode import utils, gff


def main():
    parser = argparse.ArgumentParser( description='Checks the CDS features against a genome sequence report non-terminal internal stops.')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-g', '--genome_fasta', type=str, required=False, help='Optional.  You must specify this unless the FASTA sequences for the molecules are embedded in the GFF')
    parser.add_argument('-p', '--print_n_with_stops', type=int, required=False, default=0, help='Optional.  Pass the number of sequences with internal stops you want printed (usually for debugging purposes)' )
    parser.add_argument('-o', '--output_fasta', type=str, required=False, help='Optional.  Writes an output (translated) FASTA file for all those features which had internal stops')
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_file)

    # deal with the FASTA file if the user passed one
    if args.genome_fasta is not None:
        utils.add_assembly_fasta(assemblies, args.genome_fasta)

    total_mRNAs = 0
    mRNAs_with_stops = 0

    # If this is set to the ID of any particular mRNA feature, the CDS and translation will be printed for it.
    debug_mRNA = None

    fasta_out_fh = None
    
    if args.output_fasta is not None:
        fasta_out_fh = open(args.output_fasta, 'wt')
        
    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            for mRNA in gene.mRNAs():
                coding_seq = mRNA.get_CDS_residues()
                total_mRNAs += 1

                if debug_mRNA is not None and mRNA.id == debug_mRNA:
                    print("CDS:{0}".format(coding_seq))

                if utils.translate(coding_seq).rstrip('*').count('*') > 0:
                    mRNAs_with_stops += 1
                    translated_seq = utils.translate(coding_seq)

                    if fasta_out_fh is not None:
                        loc = mRNA.location_on(assemblies[assembly_id])
                        fasta_out_fh.write(">{0} {1} {2}-{3} ({4})\n".format(mRNA.id, assembly_id, loc.fmin + 1, loc.fmax, loc.strand) )
                        fasta_out_fh.write("{0}\n".format(utils.wrapped_fasta(translated_seq)))
                    
                    if debug_mRNA is not None and mRNA.id == debug_mRNA:
                        print("TRANSLATION WITH STOP ({1}): {0}".format(translated_seq, mRNA.id) )

                    if mRNAs_with_stops <= args.print_n_with_stops:
                        print("\nmRNA id: {0}".format(mRNA.id) )
                        print("\tCDS:{0}".format(coding_seq))
                        print("\tTRANSLATION WITH STOP ({1}): {0}".format(translated_seq, mRNA.id) )


    print("\nTotal mRNAs found:{0}".format(total_mRNAs))
    print("mRNAs with embedded stops: {0}".format(mRNAs_with_stops))


if __name__ == '__main__':
    main()






