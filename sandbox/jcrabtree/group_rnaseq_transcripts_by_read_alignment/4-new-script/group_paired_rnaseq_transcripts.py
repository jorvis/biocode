#!/usr/bin/env python3

"""

Second part of proposed replacement for group_rnaseq_transcripts_by_read_alignment.py
(https://github.com/jorvis/biocode/blob/master/sandbox/jorvis/group_rnaseq_transcripts_by_read_alignment.py) 
This script reads the transcript pairs produced by find_paired_rnaseq_transcripts_by_read_alignment.py 
and runs a connected components analysis to determine which transcripts should be grouped together.
The transcript groups are output in a deterministic order to facilitate later comparisons.

"""

import argparse
import networkx as nx
import pprint
import string
import sys

def main():
    parser = argparse.ArgumentParser( description='Groups transcripts paired by find_paired_rnaseq_transcripts_by_read_alignment.py')
    parser.add_argument('-i', '--input_pair_file', type=str, required=False, help='Output file from find_paired_rnaseq_transcripts_by_read_alignment.py. Reads from stdin if not specified.' )
    args = parser.parse_args()

    sys.stderr.write("INFO: parsing transcript pair file\n")

    # read from stdin if --input_pair_file not specified
    fh = sys.stdin
    if args.input_pair_file != None:
        fh = open(args.input_pair_file)

    # create new undirected graph object
    g = nx.Graph()
    n_lines = 0

    # read transcript pairs and convert them into a NetworkX graph
    for line in fh:
        (t1, t2) = line.rstrip().split(':')
        g.add_edge(t1,t2)
        n_lines += 1

    sys.stderr.write("INFO: read " + str(n_lines) + " from input pair file\n")
    sys.stderr.write("INFO: built graph with " + str(len(g)) + " node(s) and " + str(len(g.edges())) + " edge(s)\n")

    # compute connected components
    ccs = nx.connected_components(g)
    n_groups = 0

    # print the groups, making sure that they are sorted (both within and between groups) so that 
    # it's possible to directly diff the output group files
    groups = []
    for cc in ccs:
        snodes = sorted(cc)
        groups.append({"n_nodes": len(snodes), "group": " ".join(snodes)})

    for g in sorted(groups, key=lambda x: (x["n_nodes"], x["group"]), reverse=False):
        print(g["group"])
        n_groups += 1

    sys.stderr.write("INFO: printed " + str(n_groups) + " groups/connected components\n")

if __name__ == '__main__':
    main()
