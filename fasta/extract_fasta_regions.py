#!/usr/bin/env python3

'''
Description:

Accepts a tab-delimited input file that has at least molecule ID, 
start and stop coordinates.  You must specify the column numbers
in which each of these are found.  You can optionally pass a column
that contains identifiers which will be used in the header of each
sequence region extracted.

Input: 

- A (multi-)FASTA file.  The IDs in the header (up to the first
whitespace) must match the identifiers in the tab file.
- A tab file of sequence ranges to extract.  The order can be  user
specified, but these columns must be present:
    - molecule ID
    - start position
    - end position

Optionally, another column can be used to provide a label for the 
extracted range, will which appear in the FASTA header.  If not
specified one will be fabricated from the molecule + coordinates.

Column numbers provided are numbered starting at 1

The coordinates in the file are assumed be of the human-readable sort
and are internally transformed into 0-based, interbase.  A future
feature will be to add direct support for this coordinate system.

Output:

Output will be multi-FASTA data printed to STDOUT or a file if you
pass the -o option.

Examples:


'''

import argparse
import sys

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Extract regions from a multi-FASTA file')

    ## output file to be written
    parser.add_argument('-f', '--fasta_file', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-c', '--coords_file', type=str, required=True, help='Path to a tab-delimited file with coordinates' )
    parser.add_argument('-m', '--mol_col', type=int, required=True, help='Tabdel file column with molecule identifiers' )
    parser.add_argument('-x', '--start_coord_col', type=int, required=True, help='Tabdel file column with coordinate start positions' )
    parser.add_argument('-y', '--stop_coord_col', type=int, required=True, help='Tabdel file column with coordinate stop positions' )
    parser.add_argument('-n', '--name_col', type=int, required=False, default=None, help='Optional tabdel file column with name for exported fragment' )
    parser.add_argument('-o', '--output_file', type=str, required=False, default=None, help='Optional Path to an output file to be created' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output_file is not None:
        fout = open(args.output_file, 'wt')
    
    seqs = utils.fasta_dict_from_file(args.fasta_file)

    start_col = args.start_coord_col - 1
    stop_col  = args.stop_coord_col - 1
    mol_col   = args.mol_col - 1

    for line in open(args.coords_file):
        line = line.rstrip()
        cols = line.split('\t')

        if len(cols) < 3:
            continue

        (fmin, fmax, strand) = utils.humancoords_to_0interbase(int(cols[start_col]), int(cols[stop_col]))
        mol_id = cols[mol_col]

        if mol_id not in seqs:
            raise Exception("ERROR: molecule ID ({0}) not found in FASTA file".format(mol_id))

        seq = seqs[mol_id]['s'][fmin:fmax]

        seq_id = None
        if args.name_col is None:
            seq_id = "{0}___{1}.{2}.{3}".format( mol_id, fmin, fmax, strand  )
        else:
            seq_id = cols[int(args.name_col) - 1]

        if strand == -1:
            seq = utils.reverse_complement(seq)
        
        ## write this sequence, 60bp per line
        fout.write(">{0}\n".format(seq_id))
        for i in range(0, len(seq), 60):
            fout.write(seq[i : i + 60] + "\n")
                  

if __name__ == '__main__':
    main()







