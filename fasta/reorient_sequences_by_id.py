#!/usr/bin/env python3

'''
Description:

Use this when you have a multi-FASTA file and you want to reverse or reverse
complement only select sequences within that file.  For example, strand-specific
RNA-Seq read alignment showed that some assembled transcripts were in the wrong
orientation, and this allowed us to correct them.

Input: 

- A (multi-)FASTA file.  The IDs in the header (up to the first
whitespace) must match the identifiers in the ID file.
- An ID file with one sequence ID per line.
- An action, either 'reverse' or 'revcomp'

Output:

Output will be multi-FASTA data printed to STDOUT or a file if you
pass the -o option.

Those sequences with entries within the ID file will be processed, the rest
will just be printed back out.

'''

import argparse
import sys

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Reverse or reverse-complement selected sequences within a multi-FASTA')

    ## output file to be written
    parser.add_argument('-f', '--fasta_file', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-i', '--id_file', type=str, required=True, help='Path to file with IDs to process' )
    parser.add_argument('-a', '--action', type=str, required=True, choices=['reverse', 'revcomp'], help='What should be done to the sequences in the ID file' )
    parser.add_argument('-o', '--output_file', type=str, required=False, default=None, help='Optional Path to an output file to be created' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output_file is not None:
        fout = open(args.output_file, 'wt')
    
    seqs = utils.fasta_dict_from_file(args.fasta_file)

    ids = list()

    for line in open(args.id_file):
        line = line.rstrip()
        ids.append(line)

    for seq_id in seqs:
        seq = seqs[seq_id]

        if seq_id in ids:
            if args.action == 'reverse':
                seq['s'] = seq['s'][::-1]
            elif args.action == 'revcomp':
                seq['s'] = utils.reverse_complement(seq['s'])
        
        ## write this sequence, 60bp per line
        fout.write(">{0}\n".format(seq_id))
        for i in range(0, len(seq['s']), 60):
            fout.write(seq['s'][i : i + 60] + "\n")
                  

if __name__ == '__main__':
    main()







