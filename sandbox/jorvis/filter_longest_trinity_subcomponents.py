#!/usr/bin/env python3

'''
DESCRIPTION

Written in response to a request on the Trinity user's list, asking for a script to export the
longest isoform for each component.  Offered with the caveat so well explained by "Pseudonym":

http://seqanswers.com/forums/showthread.php?t=19129

  The output of Inchworm is the set of contigs. The output of Chrysalis is the set of components,
  plus the reads which are called as belonging to those components. The output of Butterfly is the
  called transcripts for each component.

  When Butterfly rebuilds a graph for each component and does cleanup, it sometimes finds that
  the resulting graph is disconnected. This is usually because the inital contigs discovered by
  Inchworm were not "real" contigs. Each connected component in the rebuilt graph is called a
  "subcomponent".

  So the "comp" is the component, "c" is the subcomponent, and "seq" is the extracted sequence
  from the subcomponent.

  Trinity does not reason in terms of genes, loci and alternative splicing events. It solves a
  graph problem, though of course the heuristics are tuned to the needs of biology. So while it's
  highly likely that all the isoforms of a given gene belong to the same subcomponent, you
  shouldn't assume that a subcomponent is a gene.

INPUT:

- A (multi-)FASTA file, usually the Trinity.fasta file from an assembly
- The IDs of each entry in the file should follow the standard convention:

  >comp0_c0_seq1 len=281 path=[1174:0-280]
  >comp2_c0_seq1 len=222 path=[2222:0-221]
  >comp2_c1_seq1 len=225 path=[41:0-224]
  >comp3_c0_seq1 len=237 path=[2021:0-236]
  >comp4_c0_seq1 len=326 path=[56:0-325]

'''

import argparse
import os
import re
import sys
import utils

def main():
    parser = argparse.ArgumentParser( description='Filters trinity output for longest subcomponents based on naming convention')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-o', '--output', type=str, required=False, help='Output file to be created.  Default = STDOUT' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output is not None:
        fout = open(args.output, 'wt')

    seqs = utils.fasta_dict_from_file(args.input)

    components = dict()

    for seq_id in seqs:
        m = re.search("(comp\d+)_", seq_id)
        if m:
            component_id = m.group(1)

            if component_id not in components or len(seqs[seq_id]['s']) > len(components[component_id]['s']):
                components[component_id] = seqs[seq_id]
                components[component_id]['longest_id'] = seq_id
        else:
            raise Exception("ERROR: This ID wasn't in the expected format of compN_cN_seqN: {0}".format(seq_id))

    for c_id in components:
        seq_wrapped = utils.wrapped_fasta(components[c_id]['s'], every=60)
        fout.write(">{0} {1}\n{2}\n".format(components[c_id]['longest_id'], components[c_id]['h'], seq_wrapped))


if __name__ == '__main__':
    main()







