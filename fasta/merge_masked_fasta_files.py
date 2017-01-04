#!/usr/bin/env python3

"""
If you parallelize tools like Repeatmasker to run on different libraries, you
need a way to merge the results in the end into a single masked FASTA file.  This
script wears that cape.

The amount of RAM used here is proportional to the sequence content of two of your
input files.  This is wonderful if you're passing ten, less so if you're only passing
two.

This script checks that the identifiers are shared between the files and that the
residue lengths for each with shared IDs are the same.

NOTE: This only currently supports hard-masked input ('N' characters)
"""

import argparse
import sys

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Merge masked FASTA files')

    ## output file to be written
    parser.add_argument('fasta_files', metavar='N', type=str, nargs='+', help='Pass one or more FASTA files')
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    ofh = sys.stdout
    if args.output_file is not None:
        ofh = open(args.output_file, 'wt')

    files = args.fasta_files

    # pull off a file and index it
    seqs = utils.fasta_dict_from_file(files.pop())

    # python strings are immutable, so we need to transform these into lists
    for seq_id in seqs:
        seqs[seq_id]['s'] = list(seqs[seq_id]['s'])

    for fasta_file in args.fasta_files:
        new_seqs = utils.fasta_dict_from_file(fasta_file)

        for seq_id in new_seqs:
            # make sure it exists in the source file
            if seq_id not in seqs:
                raise Exception("ERROR: Seq ID {0} was found in file {1} but not in the seed file".format(seq_id, fasta_file) )

            # they should also be the same length
            if len(seqs[seq_id]) != len(new_seqs[seq_id]):
                raise Exception("ERROR: Seq ID {0} was found in {1} and the seed file but had different lengths".format(seq_id, fasta_file))

            i = 0
            for base in new_seqs[seq_id]['s']:
                if base != seqs[seq_id]['s'][i]:
                    if base == 'N':
                        seqs[seq_id]['s'][i] = 'N'
                    elif seqs[seq_id]['s'][i] != 'N':
                        print("WARNING: Disagreement {0}-{1} at position {2}".format(base, seqs[seq_id]['s'][i], i) )

                i += 1

    # now done, print out the results
    for seq_id in seqs:
        ofh.write( ">{0} {1}\n{2}\n".format(seq_id, seqs[seq_id]['h'], utils.wrapped_fasta(''.join(seqs[seq_id]['s']))))


if __name__ == '__main__':
    main()







