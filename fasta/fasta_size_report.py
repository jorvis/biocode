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
        If passed, also exports summary metrics

    --help,-h
        This help message

Output:
    Tab delimited list of fasta file info like following example. Optionally, if '-s' is set it also exports summary metrics.

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

    # Convert file argument into parsable list
    file_list = handle_file_input(args)

    # Process input fasta file(s)
    records = []
    for fa in file_list:

        # Check for valid file, but skip gracefully if missing
        file = f"{os.getcwd()}/{str(fa)}"
        if not os.path.exists(file):
            print(f"The fasta file {file} was not found! Skipping..")
            continue

        records.extend(parse_file(file))

    # Print metrics
    print("\t".join(["ID", "SIZE", "DEFLINE", "FILE"]))
    total_records = len(records)
    total_residues = 0
    for record in records:
        total_residues += int(record[1])
        print("\t".join(record))

    if args.summary:
        print(f"Summary:\n\tRecords:\t{total_records}\n\tResidues:\t{total_residues}")


# Given a single file, parse and calculate metrics
def parse_file(fasta):

    record = []
    with open(fasta, "r") as fa_file:
        short_file = fasta.split("/")[-1]
        length = 0
        defline = ""
        id = ""

        for line in fa_file:
            line = line.strip()

            if line.startswith(">"):

                if length > 0:
                    record.append([id, str(length), defline, short_file])
                    length = 0

                defline = str(line)

                id = line.lstrip(">").split(" ")[0]
            else:
                length += len(line)

        record.append([id, str(length), defline, short_file])

    return record


# Based on input type, builds and returns list of file names
def handle_file_input(args):

    file_list = []

    if args.dir:

        assert os.path.exists(args.dir), f"Directory not found at location {args.dir}"

        for file in os.listdir(args.dir):
            if file.endswith((".fa", ".fasta")):
                file_list.append(file)

    elif args.list:

        assert os.path.exists(args.list), f"List file {args.list} was not found!"

        with open(args.list, "r") as fa_list:
            for line in fa_list:
                file_list.append(line.strip())

    else:
        file_list.extend(args.file.split(" "))

    assert len(file_list) > 0, "At least one input file is required!"

    return file_list


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prints a synopsis file containing the ID, the size of the sequence, the original defline, and the originating file')

    # Input file options
    input_files = parser.add_mutually_exclusive_group(required=True)
    input_files.add_argument('-f', '--file', type=str, help='Path to an input file to be read' )
    input_files.add_argument('-l', '--list', type=str, help='Path to a file containing a list of fastas' )
    input_files.add_argument('-d', '--dir', type=str, help='Path to a directory containing fasta files' )

    # Optional flags
    parser.add_argument('-s', '--summary', required=False, help='Generate summary lines in report', action="store_true")

    args = parser.parse_args()
    main(args)
