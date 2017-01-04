#!/usr/bin/env python3

"""
Expects output from the biocode/sam_bam/report_coverage_gaps.py script

Molecules will not be exported in the same order in which they were
found in the --fasta_file.
"""

import argparse

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Splits FASTA file based on reported coverage gaps')

    ## output file to be written
    parser.add_argument('-g', '--gaps_file', type=str, required=True, help='Path to an input gaps file to be read' )
    parser.add_argument('-f', '--fasta_file', type=str, required=True, help='Path to an input FASTA file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-mfl', '--min_fragment_length', type=int, required=False, help='Min length required for a fragment to be exported' )
    parser.add_argument('-mgl', '--min_gap_length', type=int, required=False, help='Ignore gaps reported under this min length' )
    args = parser.parse_args()

    fasta = utils.fasta_dict_from_file(args.fasta_file)

    # this is just to keep track of which we've exported
    molecules_split = dict()
    ofh = open(args.output_file, 'wt')

    last_molecule_id = None
    last_end_coordinate = None

    for line in open(args.gaps_file):
        cols = line.split('\t')
        mol_id, start, stop = cols[0], int(cols[1]), int(cols[2])

        # skip if this is too short
        if args.min_gap_length is not None and (stop - start + 1) < args.min_gap_length:
            #print("DEBUG: skipping short gap {0} : {1}-{2}".format(mol_id, start, stop))
            continue

        if last_molecule_id is None:
            # first entry, export only beginning of molecule to gap start
            export_fragment(ofh, fasta, mol_id, 1, start - 1, args.min_fragment_length, molecules_split)
            last_molecule_id = mol_id
            last_end_coordinate = stop
            
        elif mol_id != last_molecule_id:
            # new molecule, export end of last molecule
            last_molecule_length = len(fasta[last_molecule_id]['s'])
            export_fragment(ofh, fasta, last_molecule_id, last_end_coordinate + 1, last_molecule_length, args.min_fragment_length, molecules_split)
            # now export the beginning of this one unless the start is 1
            if start != 1:
                export_fragment(ofh, fasta, mol_id, 1, start - 1, args.min_fragment_length, molecules_split)
                
            last_molecule_id = mol_id
            last_end_coordinate = stop
        else:
            # same molecule as we just saw, but new entry for it
            # export end of last gap until beginning of this one
            export_fragment(ofh, fasta, mol_id, last_end_coordinate + 1, start - 1, args.min_fragment_length, molecules_split)
            last_molecule_id = mol_id
            last_end_coordinate = stop

    # Now export the full sequences of any which weren't split
    for id in fasta:
        if id not in molecules_split:
            ofh.write(">{0}\n{1}\n".format(id, utils.wrapped_fasta(fasta[id]['s'])))

def export_fragment(fh, fasta, id, start, stop, min_length, exported):
    # human, 1-based coordinates expected.
    adjusted_start = start - 1
    fragment = fasta[id]['s'][adjusted_start:stop]
    if min_length is None or len(fragment) >= min_length:
        #print("DEBUG: exporting {0} : {1}-{2}".format(id, start, stop))
        fh.write(">{0}__{1}-{2}\n{3}\n".format(id, start, stop, utils.wrapped_fasta(fragment)))
        exported[id] = 1
        



if __name__ == '__main__':
    main()







