#!/usr/bin/env python3

"""
Join columns like the linux command 'join' does but via indexing of the two files
so we can avoid all the difficult sorting issues.

Also assumes that each input file has a header line which is skipped unless the
option --no-header is passed.

Warning: this approach requires the loading of file1 into RAM.

OUTPUT

The output will be all the columns from file1 followed by all those from file2, only
for those rows which joined.
"""

import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser( description='Hash-based implementation (limited) of linux join')

    parser.add_argument('input_file', nargs=2, help='<options> file1 file2')
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    parser.add_argument('-n', '--no_header', action='store_true', help="Pass if files don't heave header lines")
    parser.add_argument('-s', '--skip_na', action='store_true', help="Pass if it should skip rows whose index value is 'NA'")
    parser.add_argument('-1', '--file1_index', type=int, required=True, help='1-based column number of index in file 1' )
    parser.add_argument('-2', '--file2_index', type=int, required=True, help='1-based column number of index in file 2' )
    args = parser.parse_args()

    if args.output_file is None:
        ofh = sys.stdout
    else:
        ofh = open(args.output_file, 'wt')

    f1_idx = dict()
    f1_header = None
    line_num = 0
    
    for line in open(args.input_file[0]):
        line_num += 1
        line = line.rstrip()
        cols = line.split("\t")

        if not args.no_header:
            if line_num == 1:
                f1_header = cols
                continue

        f1_key = cols[args.file1_index - 1]

        if f1_key == 'NA' and args.skip_na:
            continue

        # Have we seen this index already?
        if f1_key in f1_idx:
            raise Exception("ERROR: Duplicate index ({}) in file1".format(f1_key))
        else:
            f1_idx[f1_key] = cols

    line_num = 0

    for line in open(args.input_file[1]):
        line_num += 1
        line = line.rstrip()
        cols = line.split("\t")

        if not args.no_header and line_num == 1:
            print("\t".join(f1_header + cols))
            continue

        f2_key = cols[args.file2_index - 1]

        if f2_key == 'NA' and args.skip_na:
            continue

        if f2_key in f1_idx:
            print("\t".join(f1_idx[f2_key] + cols), file=ofh)

    ofh.close()

if __name__ == '__main__':
    main()







