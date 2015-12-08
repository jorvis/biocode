#!/usr/bin/env python3

import argparse
import nx
import string
import sys

# TODO - compute and print connected components (and do both in a *deterministic* fashion)

def add_edge(g, n1, n2):
    sys.stderr.write("add edge " + n1 + " - " + n2 + "\n")
    

def main():
    parser = argparse.ArgumentParser( description='Groups transcripts paired by find_paired_rnaseq_transcripts_by_read_alignment.py')
    parser.add_argument('-i', '--input_pair_file', type=str, required=False, help='Output file from find_paired_rnaseq_transcripts_by_read_alignment.py' )
    args = parser.parse_args()

    sys.stderr.write("INFO: parsing transcript pair file\n")
    
    fh = sys.stdin
    if args.input_pair_file != None:
        fh = open(args.input_pair_file)

    g = nx.Graph()

    # read pairs and convert them into a NetworkX-compatible graph
    for line in fh:
        (t1, t2) = line.rstrip().split(':')
        add_edge(g, t1,t2)

if __name__ == '__main__':
    main()
