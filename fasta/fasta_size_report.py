#!/usr/bin/env python3

'''

Description:

    Accepts an input fasta file(s) and prints a synopsis file containing the ID, the size of the sequence, the original defline, and the originating file.


Usage Options:

    --file,-f
        Single or multiple fasta files of interest
        -f /path/to/file.fa /path/to/file.fa

    --list,-l
        Input list containing a set of fasta files
        -l /path/to/list.txt

    --dir,-d
        A directory containing .fsa, .fasta, .fa files
        -d /path/to/fa_files/

    --summary,-s    (Optional)
        If passed, also exports summary lines (each begins with a # symbol)

    --log,-l
        Log file

    --help,-h
        This help message

Output:
    Tab delimited list of fasta file info like following example. Optionally, if '-s' is set it also exports summary lines, each beginning with a comment.

    ID    SIZE    DEFLINE    FILE
    1     100     >1 blah    /path/to/A.fasta

Original Perl Author:

    Anu Ganapathy
    alenas@gmail.com

Python Translation:

    Chris
    RPINerd@gmail.com

'''

import argparse
import os


def main():
    parser = argparse.ArgumentParser(description='Prints a synopsis file containing the ID, the size of the sequence, the original defline, and the originating file')

    # Input file options
    parser.add_mutually_exclusive_group()
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )

    # Optional flags
    parser.add_argument('-l', '--log', required=False, help='Generate a log file if set', action="store_true")
    parser.add_argument('-s', '--summary', required=False, help='Generate summary lines in report', action="store_true")
    args = parser.parse_args()

    # Your code goes here


if __name__ == '__main__':
    main()







