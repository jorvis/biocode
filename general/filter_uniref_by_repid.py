#!/usr/bin/env python3

"""
Filter a UniRef FASTA release by identifiers in the header.  For example, if you have
a list of IDs like 'A7IAG2_METB6', this corresponds to this UniRef entry:

>UniRef100_A7IAG2 Farnesyltranstransferase n=1 Tax=Methanoregula boonei (strain DSM 21) TaxID=456442 RepID=A7IAG2_METB6

The RepID key is what is searched to match the ID list.

"""

import argparse
import os
import re

def main():
    parser = argparse.ArgumentParser( description='Filter a Uniref FASTA file by taxonomy')

    ## output file to be written
    parser.add_argument('-i', '--input_fasta', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_fasta', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-idl', '--id_list', type=str, required=True, help='Plain text file with on ID per line' )
    args = parser.parse_args()

    keep = dict()
    for line in open(args.id_list):
        line = line.rstrip()
        keep[line] = True

    fout = open(args.output_fasta, 'wt')
    keep_entry = False

    for line in open(args.input_fasta):
        if line[0] == '>':
            
            m = re.search('RepID=(\S+)', line)
            if m:
                rep_id = m.group(1)

                if rep_id in keep:
                    keep_entry = True
                else:
                    keep_entry = False
            else:
                keep_entry = False

        if keep_entry:
            fout.write(line)
        
    fout.close()

if __name__ == '__main__':
    main()







