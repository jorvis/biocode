#!/usr/bin/env python3

'''
This script examines each feature in a GFF/GTF file and validates that
the coordinates of each are within range of the referenced molecule on
which that feature is located.

An example of where this might be necessary is if an underlying assembly
has been updated but the feature identifiers were not versioned (terrible.)

Follow the GFF3 specification!

Author:  Joshua Orvis
'''

import argparse

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Checks the CDS features against a genome sequence to report/correct phase columns.')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-g', '--genome_fasta', type=str, required=False, help='Optional.  You must specify this unless the FASTA sequences for the molecules are embedded in the GFF')
    args = parser.parse_args()

    ## features, structure like:
    # [ {line:$entire_row, mol_id:tp.assembly.2342.1, start:50000, end:55611}, ... ]
    features = list()
    molecules = dict()
    in_fasta_section = False
    current_fasta_id = None
    
    for line in open(args.input_file):
        if line.startswith('##FASTA'):
            in_fasta_section = True
            continue

        if in_fasta_section:
            m = re.search('>(\S+)\s*(.*)', line)
            if m:
                current_fasta_id = m.group(1)
            else:
                # python 2.6+ makes string concatenation amortized O(n)
                molecules[current_fasta_id] += str(line.rstrip())
            
        else:
            cols = line.split("\t")

            if line.startswith('#') or len(cols) != 9:
                continue

            features.append( {'line':line, 'mol_id':cols[0], 'start':int(cols[3]), 'end':int(cols[4])} )
            molecules[cols[0]] = ''
    
    # deal with the FASTA file if the user passed one
    if args.genome_fasta is not None:
        process_fasta(molecules, args.genome_fasta)

    error_count = 0

    for feat in features:
        mol_id = feat['mol_id']
        # make sure we know about this molecule and have a seq for it
        if len(molecules[mol_id]) == 0:
            print("ERROR: This line references molecule {0} which wasn't found in FASTA data\n{1}".format(mol_id, feat['line']))
            error_count += 1

        elif feat['end'] > len(molecules[mol_id]):
            print("ERROR: This line contains a feature with end coordinates > molecule length ({1})\n{0}".format(feat['line'],len(molecules[mol_id])))
            error_count += 1


    ########################
    # Because, grammar.
    if error_count == 1:
        print("1 feature was found outside of its molecule's range")
    else:
        print("{0} features were found outside of their molecule's range".format(error_count))
        
            
def process_fasta(mols, fasta_file):
    fasta_seqs = utils.fasta_dict_from_file(fasta_file)

    for mol_id in mols:
        # check if the FASTA file provides sequence for this
        if mol_id in fasta_seqs:
            mols[mol_id] = fasta_seqs[mol_id]['s']

        

if __name__ == '__main__':
    main()






