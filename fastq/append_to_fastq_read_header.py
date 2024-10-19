#!/usr/bin/env python3

"""
Running with no options other than -i and -o, this script will transform headers like these:

    @SN7001163:78:C0YG5ACXX:6:1101:1129:2043 1:N:0:CCTAGGT

to headers like this:

    @SN7001163:78:C0YG5ACXX:6:1101:1129:2043/1

That is, it will see the '1' or '2' after the whitespace and transform the header to this instead (in FASTA):

    >SN7001163:78:C0YG5ACXX:6:1101:1129:2043/1

Otherwise, if you want a specific string to be appended to all headers, use the --append option.

Author: Joshua Orvis (jorvis AT gmail)

"""

import argparse
import gzip
import re
import sys


def parse_args() -> argparse.Namespace:
    """
    Simple argument parser for this script.

    :return: argparse.Namespace object containing the input options
    """

    parser = argparse.ArgumentParser(description="Append strings to the end of FASTQ read headers")

    ## output file to be written
    parser.add_argument(
        "-i", "--input_file", type=str, required=True, help="Path to an input file to be read (gzip supported)"
    )
    parser.add_argument("-o", "--output_file", type=str, required=False, help="Path to an output file to be created")
    parser.add_argument("-a", "--append", type=str, required=False, help="String to be appended to each read header")

    return parser.parse_args()


def append_header(line: str, append: str) -> str:
    """
    Append a custom string to the end of a FASTQ header line.

    :param line: FASTQ header line
    :param append: string to append

    :return: FASTQ header line with appended string
    """

    header = re.search("^(\@\S+)", line)
    if header:
        return f"{header.group(1)}{append}\n"
    else:
        raise Exception(f"ERROR: FASTQ header line found that didn't match expected format: {line}")


def reformat_header(line: str) -> str:
    """
    Reformat a FASTQ header line to the standard format.

    :param line: FASTQ header line

    :return: reformatted FASTQ header line
    """

    header = re.search("^(\@.+?) (\d)", line)
    if header:
        return f"{header.group(1)}/{header.group(2)}\n"
    else:
        raise Exception(f"ERROR: FASTQ header line found that didn't match expected format: {line}")


def main(args: argparse.Namespace) -> None:

    if args.output_file is None:
        ofh = sys.stdout
    else:
        ofh = open(args.output_file, "wt")

    line_count = 0

    if args.input_file.endswith(".gz"):
        fh = gzip.open(args.input_file, "rb")
        is_compressed = True
    else:
        fh = open(args.input_file, "rU")
        is_compressed = False

    for line in fh:
        if is_compressed:
            line = line.decode()

        line_count += 1

        if line_count % 4 == 1:
            header = append_header(line, args.append) if args.append else reformat_header(line)
            ofh.write(header)
        else:
            ofh.write(line)

    fh.close()
    ofh.close()


if __name__ == "__main__":

    args = parse_args()
    main(args)
