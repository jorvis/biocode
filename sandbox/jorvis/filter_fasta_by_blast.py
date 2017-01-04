#!/usr/bin/env python3.2

import argparse
import re

from biocode.utils import read_list_file


def process_blast_file( file, ids, args ):
    ''' Parses through a tab-delimited BLAST file and saves all IDs matching 
        the cutoff into the ids array
    '''
    passed_cutoff_count = 0
    
    for line in open(file, 'r'):
        line.rstrip()
        line_passed = 1
        keep_entry = 1

        cols = line.split("\t")
        if len(cols) == 12:
            qry_id = cols[0]
            ## have we seen this qry yet?
            if qry_id in ids:
                if ids[qry_id] == 1:
                    ## we already found a passing match for this query, don't override it
                    continue
            else:
                ids[qry_id] = 0
            
            ## check for an e-val cutoff
            if args.eval_cutoff is not None and float(cols[10]) < float(args.eval_cutoff):
                ids[qry_id] = 1
                passed_cutoff_count += 1

    return passed_cutoff_count;


def export_file( ipath, opath, mode, ids ):
    ofh = open(opath, 'w')
    keep_entry = 1
    c_entries_kept = 0
    c_total_entries = 0

    for line in open(ipath, 'r'):
        m = re.match( '>(\S+)', line )
        if m:
            acc = m.group(1);
            c_total_entries += 1

            if mode == 'include':
                if acc in ids and ids[acc] == 1:
                    keep_entry = 1
                    c_entries_kept += 1
                else:
                    keep_entry = 0
            elif mode == 'exclude':
                if acc in ids and ids[acc] == 1:
                    keep_entry = 0
                else:
                    keep_entry = 1
                    c_entries_kept += 1

        if keep_entry == 1:
            ofh.write( line )
    
    return [ c_total_entries, c_entries_kept ]



def main():
    parser = argparse.ArgumentParser( description='Parses BLAST results to include or exclude entries from a FASTA file')
    parser.add_argument('-l', '--blast_tab_list', type=str, required=True, help='A list file with the paths of all files in -m 8 output format from NCBI blast' )
    parser.add_argument('-f', '--fasta_file', type=str, required=True, help='The FASTA file to be filtered' )
    parser.add_argument('-m', '--mode', type=str, required=True, choices=['include','exclude'], help='Enter either "include" or "exclude" based on the mode of filtering you want based on the blast hits.')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-e', '--eval_cutoff', type=str, required=True, help='Filter entries with an E-value lower than this' )
    args = parser.parse_args()

    ## 3 different options for each accession in the fasta_file
    #   - have no key here - wasn't found in the BLAST matches at all
    #   - key with value of 0 - found in BLAST, but with scores not meeting cutoffs
    #   - key with value of 1 - found in BLAST, with scores meeting cutoffs
    ids_matched = {}
    c_ids_matching_cutoffs = 0

    for file in read_list_file( args.blast_tab_list ):
        c_ids_matching_cutoffs += process_blast_file(file, ids_matched, args)

    [c_total, c_kept] = export_file( args.fasta_file, args.output_file, args.mode, ids_matched )

    print("INFO: Exported {0} entries out of {1} total in the source file".format(c_kept, c_total) )

    


if __name__ == '__main__':
    main()







