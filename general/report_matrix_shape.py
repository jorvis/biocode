#!/usr/bin/env python3

"""

Given an input tab-delimited file representing a data matrix (such as rnaSeq expression counts)
this script checks that all lines in the file have the same number of columns and, if so, reports
the total column and row counts.

"""

import argparse
import os


def main():
    parser = argparse.ArgumentParser( description='Check to make sure a tab-delimited matrix file has uniform shape')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to check' )
    args = parser.parse_args()

    line_num = 0
    current_count = None

    for line in open(args.input_file):
        line_num += 1
        line = line.rstrip();
        col_count = len(line.split('\t'))

        if current_count is None:
            current_count = col_count
            continue
        else:
            if current_count != col_count:
                raise Exception(
                    "ERROR: expected column count of {0} but on line {1} got a count of {2}".format(current_count, line_num, col_count))

    print("Success. Input file is uniform matrix with {0} rows and {1} columns".format(line_num, current_count))

if __name__ == '__main__':
    main()







