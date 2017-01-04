#!/usr/bin/env python3

'''
Description:

Accepts a (multi-FASTA) file as input and reports the base/residue
composition of the sequences found within it.

Input: 

- A (multi-)FASTA file.  The IDs in the header (up to the first
whitespace) must be unique.


Output:

Output will be a tab-delimited report with the following columns:

     1. Base/residue
     2. Count of base/residue in column 1
     3. Percentage of base/residue relative to entire sequence content

The first line of output is a comment line showing the total residue
count in the file.

Example output:



This will all be printed to STDOUT or a file if you pass the -o option.

'''

import argparse
import sys
from collections import Counter

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Reports base/residue composition of a FASTA file')

    ## output file to be written
    parser.add_argument('-f', '--fasta_file', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-o', '--output_file', type=str, required=False, default=None, help='Optional Path to an output file to be created' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output_file is not None:
        fout = open(args.output_file, 'wt')
    
    seqs = utils.fasta_dict_from_file(args.fasta_file)
    
    total_residues = 0
    total_residue_counts = Counter()

    for seq_id in seqs:
        total_residues += len(seqs[seq_id]['s'])
        seq_residue_counts = Counter(seqs[seq_id]['s'])

        for residue in seq_residue_counts:
            total_residue_counts[residue] += seq_residue_counts[residue]
        
                  
    fout.write("# Total residues found: {0}\n".format(total_residues))

    for residue in total_residue_counts:
        residue_count = total_residue_counts[residue]
        residue_perc  = (residue_count / total_residues) * 100
        fout.write("{0}\t{1}\t{2}\n".format(residue, residue_count, residue_perc))

if __name__ == '__main__':
    main()







