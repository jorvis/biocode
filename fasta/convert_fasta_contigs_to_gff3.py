#!/usr/bin/env python3

"""

This script assumes that the input FASTA represents genomic sequence and creates
a GFF3 file for these, devoid of anything other than placeholder rows for each
molecule.  It is assumed other processes will add to this.

Input example (fake):

>AL009126.2 Bacillus subtilis contig 2
ATCTTTTTCGGCTTTTTTTAGTATCCACAGAGGTTATCGACAACATTTTCACATTACCAACCCCTGTGGACAAGGTTTTT
TCAACAGGTTGTCCGCTTTGTGGATAAGATTGTGA
>AL009126.3 Bacillus subtilis contig 3
ATCTTTTTCGGCTTTTTTTAGTATCCACAGAGGTTATCGACAACATTTTCACATTACCAACCCCTGTGGACAAGGTTTTT
TCAACAGGTTGTCCGCTTTGTGGATAAGATTGTGACAACCATTGCAAGCTCTCGTTTATTTTGGTATTATATTTGTGTTT
TAACTCTTGATTACTAATCCTACCTTTCCTCTTTATCCACAAAGTGTGGATAAGTTGTGGATTGATTTCACACAGCTTGT

Output:
AL009126.2    IGS     contig      1       100613  .       .       .       ID=AL009126.2;Name=Bacillus subtilis contig 2
AL009126.3    IGS     contig      1       230218  .       .       .       ID=AL009126.3;Name=Bacillus subtilis contig 3

The value of the 2nd column is up to you and defined using the --source parameter.

"""

import argparse

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Creates a GFF3 file from a genomic FASTA')

    ## output file to be written
    parser.add_argument('-i', '--input_fasta', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_gff3', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-s', '--source', type=str, required=True, help='Source, fills column 2 of the output GFF3 file' )
    parser.add_argument('-t', '--molecule_term', type=str, required=False, default='contig', help='SO term to represent genomic sequence type')
    args = parser.parse_args()

    ofh = open(args.output_gff3, 'wt')
    seqs = utils.fasta_dict_from_file(args.input_fasta)

    # header
    ofh.write("##gff-version 3\n")

    for seq_id in seqs:
        ofh.write("{0}\t{1}\t{2}\t1\t{3}\t.\t.\t.\tID={0}".format(seq_id, args.source, args.molecule_term, len(seqs[seq_id]['s'])))
        
        if len(seqs[seq_id]['h']) > 0:
            ofh.write(";Name={0}".format(seqs[seq_id]['h']))

        ofh.write("\n")
            


if __name__ == '__main__':
    main()







