#!/usr/bin/env python3

"""

Read BAM file, split each line into columns, build single NON-nested dict of transcript ids.
  
"""

import argparse
import os
import re
import sys


def main():
    parser = argparse.ArgumentParser( description='Groups transcripts by mapped read pairs')
    parser.add_argument('-i', '--input_sam_file', type=str, required=False, help='SAM file of read alignments back to transcripts' )
    args = parser.parse_args()

    single_read_pairings = 0
    same_read_pairings = 0
    pairings = dict()

    print("INFO: parsing SAM file and creating transcript pairings")
    fh = sys.stdin
    if args.input_sam_file != None:
        fh = open(args.input_sam_file)
    n_lines = 0

    for line in fh:
        n_lines = n_lines + 1
        if line[0] == '@': continue

        cols = line.split("\t")
        ref_read = cols[0]
        ref_transcript = cols[2]
        other_transcript = cols[6]

        if (n_lines % 1000000) == 0:
            print("n_lines=" + str(n_lines))

        # we don't care about the lines where the mate is unmapped or mapped to the same
        if other_transcript == '*':
            single_read_pairings += 1
            continue
        elif other_transcript == '=':
            same_read_pairings += 1
            continue

        transcripts = sorted([ref_transcript, other_transcript])
        key = transcripts[0] + ":" + transcripts[1]
        
        if key not in pairings:
            pairings[key] = list()

    print ("read " + str(n_lines) + " line(s) from file, " + str(len(pairings)) + " pairings")

if __name__ == '__main__':
    main()







