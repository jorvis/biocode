#!/usr/bin/env python3

import argparse
import sys

from biocode import utils

'''
Description:

There are many genome files in GenBank with abnormally-long homopolymeric repeats of non-N
bases.  For example:

>gi|257136525|ref|NZ_GG699286.1| Xanthomonas campestris pv. vasculorum NCPPB702 genomic scaffold scf_7293_715, whole genome shotgun sequence

Examples of the homopolymeric stretches:

WARNING: Sequence ID gi|257136525|ref|NZ_GG699286.1| contains a homopolymer run (T) of length 45972
WARNING: Sequence ID gi|257136529|ref|NZ_GG699290.1| contains a homopolymer run (A) of length 131072
WARNING: Sequence ID gi|257136550|ref|NZ_GG699311.1| contains a homopolymer run (T) of length 51385
WARNING: Sequence ID gi|257136550|ref|NZ_GG699311.1| contains a homopolymer run (A) of length 262144
WARNING: Sequence ID gi|257136567|ref|NZ_GG699328.1| contains a homopolymer run (A) of length 61064

This causes errors with some tools (such as bowtie2), so this script can be used to replace
these with Ns, which are usually better-tolerated.

WARNINGS:

I didn't go with a memory-efficient implementation here, so this version of the script loads the
input genomes fully into memory during processing.

'''

def main():
    parser = argparse.ArgumentParser( description='Replaces long homopolymeric stretches with N characters')

    parser.add_argument('-i', '--input', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-o', '--output', type=str, required=False, help='Path to an output FASTA file to be created' )
    parser.add_argument('-hll', '--homopolymer_length_limit', type=int, required=True, help='Stretches of non-N residues longer than this will be replaced with Ns' )
    args = parser.parse_args()

    if args.output is None:
        out_fh = sys.stdout
    else:
        out_fh = open( args.output, 'wt' )

    sys.stderr.write("INFO: Parsing input FASTA\n")
    sys.stderr.flush()
    seqs = utils.fasta_dict_from_file(args.input)

    sys.stderr.write("INFO: Looking for homopolymeric runs > {0} bp\n".format(args.homopolymer_length_limit))
    sys.stderr.flush()
    for seq_id in seqs:
        seq = seqs[seq_id]
        current_seq = seq['s']
        current_homopolymer_base = None
        current_homopolymer_length = 0
        current_homopolymer_start_idx = 0
        base_index = 0

        for base in list(seq['s']):
            if base == current_homopolymer_base:
                current_homopolymer_length += 1
            else:
                if current_homopolymer_length > args.homopolymer_length_limit and current_homopolymer_base != 'N':
                    sys.stderr.write("WARNING: Replacing {3} bp of {2}s in Sequence ID {0} starting at position {1}\n".format(
                        seq_id, current_homopolymer_start_idx + 1, current_homopolymer_base, current_homopolymer_length))
                    sys.stderr.flush()

                    current_seq = "{0}{1}{2}".format(seq['s'][0:current_homopolymer_start_idx],
                                                     'N' * current_homopolymer_length,
                                                     seq['s'][base_index:])

                current_homopolymer_base = base
                current_homopolymer_length = 1
                current_homopolymer_start_idx = base_index

            base_index += 1

        ## check after the last row for any runs which terminate the sequence
        if current_homopolymer_length > args.homopolymer_length_limit and current_homopolymer_base != 'N':
             sys.stderr.write("WARNING: Replacing {3} bp of {2} bases in Sequence ID {0} starting at position {1}\n".format(
                 seq_id, current_homopolymer_start_idx, current_homopolymer_base, current_homopolymer_length))
             sys.stderr.flush()

             current_seq = "{0}{1}{2}".format(current_seq[0:current_homopolymer_start_idx],
                                              'N' * current_homopolymer_length,
                                              current_seq[base_index:])

        seqs[seq_id]['s'] = current_seq
        out_fh.write(">{0} {1}\n".format(seq_id, seqs[seq_id]['h']))
        out_fh.write(utils.wrapped_fasta(seqs[seq_id]['s']))
        out_fh.write("\n")
            

if __name__ == '__main__':
    main()







