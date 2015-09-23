#!/usr/bin/env python3

"""

--lab_name:
  Documentation: http://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation#protein_id

"""

import argparse
import biocodegff
import biocodetbl
import os


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-go', '--go_obo', type=str, required=False, help='GO terms will not be exported unless you pass the path to a GO OBO file')
    parser.add_argument('-ln', '--lab_name', type=str, required=True, help='Required by NCBI to identify the submitting group' )
    args = parser.parse_args()

    ## Do a full round-trip of the file and see how there are differences
    (assemblies, features) = biocodegff.get_gff3_features( args.input_file )
    
    ofh = open(args.output_file, 'wt')
    biocodetbl.print_tbl_from_assemblies(assemblies=assemblies, ofh=ofh, go_obo=args.go_obo, lab_name=args.lab_name)

if __name__ == '__main__':
    main()







