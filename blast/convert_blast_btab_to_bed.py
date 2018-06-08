#!/usr/bin/env python3

"""

The input here is supposed to be either BLAST output using the -m 8 or -m 9
options or also the tabular output from RAPSearch2.

Input:

1. Query 
2. Subject 
3. % identity
4. aln-len
5. mismatch
6. gap-openings
7. q.start
8. q.end
9. s.start
10. s.end
11. E-value (or log() for RAPSearch2)
12. bit-score

The output is a BED file with the following tabs:

<queryID> <start> <end> <matchID> <score> <strand>

"""

import argparse
import os


def main():
    parser = argparse.ArgumentParser( description='Converts BLAST and RAPSearch2 tabular output to BED')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        raise Exception("Input file passed {0} wasn't found".format(args.input_file))

    ofh = open(args.output_file, 'wt')

    for line in open(args.input_file):
        if line.startswith('#'): continue
        line = line.rstrip()
        cols = line.split("\t")
        if len(cols) < 12: continue

        strand = '+' if int(cols[6]) < int(cols[7]) else '-'
        ofh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(cols[0], int(cols[6]) - 1, cols[7], cols[1], cols[11], strand))
        

    

    ofh.close()


if __name__ == '__main__':
    main()







