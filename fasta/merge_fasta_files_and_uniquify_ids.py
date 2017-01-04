#!/usr/bin/env python3

"""
The use case for this was that I generated several different transcriptome
assemblies with TGICL, and each creates new IDs within the contigs file like
this:

   >CL1Contig1
   >CL1Contig2

But every file started with CL1Contig1.  So to merge all the contig files you
also have to have a way to uniquify the IDs.  This script does this by allowing
you to specify a base name and it generates IDs from that, appending a numeric
counter to the end, like:

   ./merge_fasta_files_and_uniquify_ids.py -l foo.list -b tgicl.contigs

It would generate IDs like:

   tgicl.contigs.1
   tgicl.contigs.2
   ...

Output will be written to STDOUT unless --output_file is passed.  Headers are retained,
if present, and there is no assumption that input IDs are unique even within individual
files.

"""

import argparse
import re
import sys

from biocode.utils import read_list_file


def main():
    parser = argparse.ArgumentParser( description='Looks for FASTA entries which seem to have been embedded within another')
    parser.add_argument('-l', '--input_list', type=str, required=False, help='A list file with the paths of all files' )
    parser.add_argument('-f', '--fasta_file', type=str, required=False, help='The single FASTA file to be processed' )
    parser.add_argument('-b', '--basename', type=str, required=True, help='Base name to be used for new FASTA files' )
    parser.add_argument('-o', '--output_file', type=str, required=False, default=None, help='Optional path to an output file to be created' )
    args = parser.parse_args()

    if args.input_list is None and args.fasta_file is None:
        raise Exception("ERROR: You must specify either --input_list or --fasta_file to define input");

    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output_file is not None:
        fout = open(args.output_file, 'wt')

    files_to_process = list()

    if args.input_list is not None:
        for file in read_list_file( args.input_list ):
            files_to_process.append(file)

    if args.fasta_file is not None:
        files_to_process.append(args.fasta_file)

    next_id_idx = 1

    for file in files_to_process:
        for line in open(file):
            if line[0] == '>':
                line = line.rstrip()
                next_id = "{0}.{1}".format(args.basename, next_id_idx)
                next_id_idx += 1
                
                m = re.match(">\S+ (.*)", line)
                if m:
                    fout.write(">{0} {1}\n".format(next_id, m.groups(1)))
                else:
                    fout.write(">{0}\n".format(next_id))
            else:
                fout.write(line)
        

if __name__ == '__main__':
    main()







