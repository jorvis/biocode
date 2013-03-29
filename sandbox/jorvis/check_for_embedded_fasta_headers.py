#!/usr/bin/env python3.2

import argparse
import os
import re
import sys

from biocodeutils import read_list_file

def main():
    parser = argparse.ArgumentParser( description='Looks for FASTA entries which seem to have been embedded within another')
    parser.add_argument('-l', '--input_list', type=str, required=False, help='A list file with the paths of all files' )
    parser.add_argument('-f', '--fasta_file', type=str, required=False, help='The single FASTA file to be processed' )
    args = parser.parse_args()

    if args.input_list is None and args.fasta_file is None:
        raise Exception("ERROR: You must specify either --input_list or --fasta_file to define input");

    if args.input_list is not None:
        for file in read_list_file( args.input_list ):
            fail_count = process_fasta_file( file )

            if ( fail_count > 0 ):
                print( "{0}\t{1}\n".format(file, fail_count) )

    if args.fasta_file is not None:
        fail_count = process_fasta_file( args.fasta_file )

        if ( fail_count > 0 ):
                print( "{0}\t{1}\n".format(file, fail_count) )
        


def process_fasta_file( file ):
    embedded_header_count = 0

    print("INFO: processing file: {0}".format(file), file=sys.stderr )

    ifh = open(file)
    
    for line in ifh:
        ## if the line contains a '>' but it's not at the beginning, flag the entry:
        if line.find('>') != -1 and not line.startswith('>'):
            embedded_header_count += 1
    
    
    return embedded_header_count


if __name__ == '__main__':
    main()







