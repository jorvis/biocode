#!/usr/bin/env python3.2

import argparse
import os
import re


def report_coverage_summary(refs, ofile_path):
    ofh = open(ofile_path, 'w')
    ofh.write("ref_id\tuncovered_bases\tref_length\n")
    
    for ref in refs:
        last_pos = 0
        c_uncovered = 0
        mol_length = refs[ref]['length']

        for range in refs[ref]['ranges']:
            ## this will only be true if the current range doesn't overlap the previous one
            if ( range[0] > last_pos ):
                c_uncovered += range[0] - last_pos

            ## now reset the last_pos, as long as this range is past it
            if (range[1] > last_pos):
                last_pos = range[1]

        ## now do the last position to the end of the molecule
        c_uncovered += mol_length - last_pos
        c_perc = ((mol_length - c_uncovered) / mol_length) * 100
        ofh.write("{0}\t{1}\t{2}\t{3:.2f}\n".format(ref, c_uncovered, mol_length, c_perc))


def main():
    ''' 
    Expectations:
    - The 4th line should have the column headers
    - show-coords is run with at least the -r, -c and -l options
    '''
    parser = argparse.ArgumentParser( description='Parses a nucmer coords file and reports a summary of the reference genome coverage')

    ## output file to be written
    parser.add_argument('-c', '--coords_file', type=str, required=True, help='Nucmer alignment coords file, run with the -r, -c and -l options' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-p', '--pct_id_cutoff', type=float, required=False, default=0.0, help='Only consider alignment ranges with a percent identity higher than this value (0-100)' )
    args = parser.parse_args()

    ## lookup for the indexes of needed columns (since show-coords options can make these change positions)
    column_idx = {}
    current_line_number = 0

    ## each element key is a reference molecule, with subkeys: ranges[], length
    references = {}
    current_ref_id = None
    
    for line in open(args.coords_file):
        line = line.strip()
        current_line_number += 1

        if ( current_line_number == 4 ):
            cols = re.split('\s{2,}', line)
            col_num = 0
            for col in cols:
                if col != '|' and col != '':
                    column_idx[col] = col_num

                col_num += 1
        elif current_line_number > 5:
            cols = line.split()
            this_ref_id = cols[ column_idx['| [TAGS]'] + 1]
            fmin = int(cols[ column_idx['[S1]'] ]) - 1
            fmax = int(cols[ column_idx['[E1]'] ])
            mol_length = int(cols[ column_idx['[LEN R]'] ])
            pct_id = float(cols[ column_idx['[% IDY]'] ])

            if current_ref_id == None or this_ref_id != current_ref_id:
                current_ref_id = this_ref_id
                references[current_ref_id] = { 'ranges': [], 'length': mol_length }

            if pct_id >= args.pct_id_cutoff:
                references[current_ref_id]['ranges'].append( [fmin, fmax] )

    report_coverage_summary(references, args.output_file)

if __name__ == '__main__':
    main()







