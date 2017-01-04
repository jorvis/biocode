#!/usr/bin/env python3

'''
Description:

This script reformats a FASTA file such that there are no more than N-characters
of sequence residues per line.  The default is 60.

Input: 

- A (multi-)FASTA file.  
- Optional width (-w) of the residue lines.  Default = 60

'''

import argparse
import sys

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Reformats a FASTA file such that there are no more than -w characters of sequence residues per line.')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-w', '--width', type=int, required=False, default=60, help='Width - number of residues per line' )
    parser.add_argument('-o', '--output', type=str, required=False, help='Output file to be created.  Default = STDOUT' )
    parser.add_argument('-uc', '--upper_case', action='store_true', required=False, help='Forces all bases to be upper-case' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output is not None:
        fout = open(args.output, 'wt')

    seqs = utils.fasta_dict_from_file(args.input)

    for seq_id in seqs:
        if args.upper_case == True:
            seqs[seq_id]['s'] = seqs[seq_id]['s'].upper()
            
        seq_wrapped = utils.wrapped_fasta(seqs[seq_id]['s'], every=args.width)
        fout.write(">{0} {1}\n{2}\n".format(seq_id, seqs[seq_id]['h'], seq_wrapped))


if __name__ == '__main__':
    main()







