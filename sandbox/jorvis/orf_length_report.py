#!/usr/bin/env python3

"""

Reports the longest ORF length in all frames.  This is only appropriate for bacterial
genes or eukaryotic mRNAs.

"""

import argparse
import os
from biocode import utils

def main():
    parser = argparse.ArgumentParser( description='Reports longest ORF length in all frames')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    args = parser.parse_args()

    seqs = utils.fasta_dict_from_file(args.input_file)

    for seq_id in seqs:
        winning_frame = None
        winning_frame_length = 0

        for frame in range(1, 4):
            offset = frame - 1
            seq = seqs[seq_id]['s'][offset:]
            polyseq = utils.translate(seq)

            longest_len = 0
            current_len = 0
            
            for base in polyseq:
                if base == '*':
                    if current_len > longest_len:
                        longest_len = current_len
                        current_len = 0
                else:
                    current_len += 1

            if current_len > longest_len:
                longest_len = current_len

            if longest_len > winning_frame_length:
                winning_frame = frame
                winning_frame_length = longest_len

        print("{0}\t{1}\t{2}".format(seq_id, winning_frame, winning_frame_length))


if __name__ == '__main__':
    main()







