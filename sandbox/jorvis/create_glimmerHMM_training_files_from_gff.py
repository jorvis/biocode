#!/usr/bin/env python3

"""
The GlimmerHMM training documentation says that two files are needed for training.  One
is a multi-fasta file of what are presumably transcript sequences (this is never stated,
so it could be CDS?) and the coordinates of the exons relative to the same sequences.

This script generates that given an input GFF3.  Currently the mRNA features are written
to the FASTA.

This script definitely assumes that features have been sorted.

"""

import argparse
import sys

from biocode import gff


def main():
    parser = argparse.ArgumentParser('Filter the genes of a GFF3 file by mRNA child IDs')

    ## output file to be written
    parser.add_argument('-i', '--input_gff', type=str, required=True, help='GFF file of source annotation' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Optional output file path (else STDOUT)' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output_file is not None:
        fout = open(args.output_file, 'wt')

    current_mRNA_id = None
    current_mol_id = None
    current_fragments = list()
    current_direction = None
        
    for line in open(args.input_gff):
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) != 9:
            continue

        # grab the ID and Parent columns if any
        id = gff.column_9_value(cols[8], 'ID')
        parent = gff.column_9_value(cols[8], 'Parent')
        mol_id = cols[0]
        type = cols[2]

        if type == 'mRNA':
            if current_mRNA_id is not None and id != current_mRNA_id:
                # purge the existing one first
                write_transcript(fout, current_mol_id, current_fragments, current_direction)
                current_fragments = list()
                
            current_mRNA_id = id
            current_mol_id = cols[0]
            current_direction = cols[6]
            
        elif type == 'exon':
            
            if cols[6] == '+':
                current_fragments.append({'start':cols[3], 'end':cols[4]})
            else:
                current_fragments.append({'start':cols[4], 'end':cols[3]})

    write_transcript(fout, current_mol_id, current_fragments, current_direction)


def write_transcript(fh, mol_id, fragments, direction):
    if direction == '-':
        fragments = reversed(fragments)
    elif direction != '+':
        raise Exception("ERROR: unrecognized direction: {0}".format(direction))

    for frag in fragments:
        fh.write("{0} {1} {2}\n".format(mol_id, frag['start'], frag['end']) )

    fh.write("\n")



if __name__ == '__main__':
    main()







