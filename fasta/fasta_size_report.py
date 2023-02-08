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


def main(args):

    # Enable logging if option set
    # TODO if log, redirect print to file
    if args.log:
        log_file = open(args.log, "w")

    # Convert file argument into parsable list
    file_list = handle_file_input(args)

    # Process input fasta file(s)
    records = []
    for fa in file_list:

        # Check for valid file, but skip gracefully if missing
        if not os.path.exists(fa):
            print(f"The fasta file {fa} was not found! Skipping..")
            continue

        records.append(parse_file(fa))

    # Print metrics
    for record in records:
        print("\t".join(record))


def parse_file(fasta):

    record = []

    return record


def handle_file_input(args):

    file_list = []

    if args.dir:

        assert os.path.exists(args.dir), f"Directory not found at location {args.dir}"

        for file in os.listdir(args.dir):
            if file.endswith(".fa", ".fasta"):
                file_list.append(file)

    elif args.list:

        assert os.path.exists(args.list), f"List file {args.list} was not found!"

        with open(args.list, "r") as fa_list:
            for line in fa_list:
                file_list.append(line.strip())

    else:
        file_list.append(args.file.split(" "))

    return file_list


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prints a synopsis file containing the ID, the size of the sequence, the original defline, and the originating file')

    # Input file options
    input_files = parser.add_mutually_exclusive_group()
    input_files.add_argument('-f', '--file', type=str, required=True, help='Path to an input file to be read' )
    input_files.add_argument('-l', '--list', type=str, required=True, help='Path to a file containing a list of fastas' )
    input_files.add_argument('-d', '--dir', type=str, required=True, help='Path to a directory containing fasta files' )

    # Optional flags
    parser.add_argument('-l', '--log', required=False, help='Generate a log file if set')
    parser.add_argument('-s', '--summary', required=False, help='Generate summary lines in report', action="store_true")
    args = parser.parse_args()

    main(args)
