#!/usr/bin/env python3

"""

Randomly extract/subsample a specified number of sequences from an
input multi-FASTA file.  If the number of sequences requested is
greater than the number of sequences in the file then all of the input
sequences will be printed.  Otherwise, the exact number of
requested sequences will be returned, in the same order that they
appear in the original file. Note that the input multi-FASTA will be
scanned twice: once to determine the total sequence count and once to
extract the randomly-chosen sequences.

Usage example:

./subsample_fasta.py -i reference_contigs.fna -n 100 > 100_random_contigs.fna

"""

import argparse
import os
from random import randint, seed
import sys

def main():
    parser = argparse.ArgumentParser( description='Randomly extract/subsample a specified number of sequences from a FASTA file.')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input FASTA file to be read' )
    parser.add_argument('-n', '--num_seqs', type=int, required=True, help='Number of sequences to extract.' )
    parser.add_argument('-s', '--seed', type=int, required=False, help='Seed for random number generation (optional).' )
    args = parser.parse_args()

    # 1st pass: count sequences
    all_seqs = []
    seqnum = 0
    for line in open(args.input_file):
        if line.startswith(">"):
            seqnum += 1
            all_seqs.append(seqnum)
    n_seqs = seqnum

    # select sequences for extraction
    selected_seqs = {}
    n_unselected = n_seqs
    if args.seed is not None:
        seed(args.seed)

    num_seqs = args.num_seqs
    if num_seqs > n_seqs:
        num_seqs = n_seqs
    for i in range(0, num_seqs):
        # choose one sequence from all_seqs and move it into selected_seqs
        chosen = randint(0, n_unselected-1)
        selected_seqs[all_seqs[chosen]] = 1
        # remove it from all_seqs
        del all_seqs[chosen]
        n_unselected -= 1

    # 2nd pass: extract selected sequences
    seqnum = 0
    for line in open(args.input_file):
        if line.startswith(">"):
            seqnum += 1
        if seqnum in selected_seqs:
            sys.stdout.write(line)

if __name__ == '__main__':
    main()

