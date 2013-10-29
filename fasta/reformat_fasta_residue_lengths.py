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
import os
import sys
import biocodeutils

def main():
    parser = argparse.ArgumentParser( description='Reformats a FASTA file such that there are no more than -w characters of sequence residues per line.')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-w', '--width', type=int, required=False, default=60, help='Width - number of residues per line' )
    parser.add_argument('-o', '--output', type=str, required=False, help='Output file to be created.  Default = STDOUT' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output is not None:
        fout = open(args.output, 'wt')

    seqs = biocodeutils.fasta_dict_from_file( args.input )

    for seq_id in seqs:
        seq_wrapped = wrapped(seqs[seq_id]['s'], every=args.width)
        fout.write(">{0} {1}\n{2}\n".format(seq_id, seqs[seq_id]['h'], seq_wrapped))


def wrapped(string, every=60):
    ''' this runs orders of magnitude faster than using the textwrap module '''
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))


if __name__ == '__main__':
    main()







