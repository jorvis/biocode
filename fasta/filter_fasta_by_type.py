#!/usr/bin/env python3

'''
Description:

This script is for those (rare?) occasions where you have a FASTA
file with mixed protein and nucleotide sequences.  It inspects
each entry in a FASTA file and uses the composition to predict the
type of sequence.

How?  It is done in a relatively naive way, checking the ATGCNX
content of the sequence.  If it's above a user-configurable cutoff
value it's determined to be nucleotide, else polypeptide.

As a warning, this was written using existing functions quickly and
is not very memory efficient.  The current implementation isn't the
best choice if you have a file with millions of sequences.


Input: 

- A (multi-)FASTA file.  
- Optional composition cutoff (-c) of ATGCNX content for each sequence
  to be considered a nucleotide (else it's a polypeptide)  (1-100). If
  not passed, the default (80) is used.

Output:

The output can be individual FASTA files for the protein sequences
using the -p option, nucleotide sequences using the -n option, 
or both.

'''

import argparse

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Split multi-FASTA file into separate protein and nucleotide files')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-p', '--protein', type=str, required=False, help='Path to a tab-delimited file with coordinates' )
    parser.add_argument('-n', '--nucleotide', type=str, required=False, help='Tabdel file column with molecule identifiers' )
    parser.add_argument('-c', '--cutoff', type=str, required=False, default=80, help='Min percent (1-100) of ATGCNX content to be considered a nucleotide sequence' )
    args = parser.parse_args()

    pout = nout = None
    
    if args.protein is not None:
        pout = open(args.protein, 'wt')

    if args.nucleotide is not None:
        nout = open(args.nucleotide, 'wt')

    ## the user should have specified at least one
    if pout is None and nout is None:
        raise Exception("ERROR: you must specify either -p or -n options (else why are you running this script?")
    
    seqs = utils.fasta_dict_from_file(args.input)

    for seq_id in seqs:
        seq = seqs[seq_id]
        seqcomp = nucleotide_composition( seq['s'] )
        seq_wrapped = wrapped(seq['s'], every=60)

        if seqcomp >= args.cutoff:
            ## it's a nucleotide
            if nout is not None:
                nout.write(">{0} {1}\n{2}\n".format(seq_id, seq['h'], seq_wrapped ) )

        else:
            ## it's a protein
            if pout is not None:
                pout.write(">{0} {1}\n{2}\n".format(seq_id, seq['h'], seq_wrapped ) )


def wrapped(string, every=60):
    ''' this runs orders of magnitude faster than using the textwrap module '''
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))
                  

def nucleotide_composition( seq ):
    seq = seq.upper()

    residue_count = 0

    for base in "ATGCNX":
        residue_count += seq.count(base)

    return (residue_count / len(seq)) * 100

if __name__ == '__main__':
    main()







