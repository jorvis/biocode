#!/usr/bin/env python3

"""This script takes two fasta files and interleaves them

Usage:
    interleave-fasta.py fasta_file1 fasta_file2

Author: Sébastien Boisvert
Py3 conversion: RPINerd
"""

import argparse
import os


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Interleave two fasta files")
    parser.add_argument("-1", "--fasta1", type=str, help="First fasta file")
    parser.add_argument("-2", "--fasta2", type=str, help="Second fasta file")
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:

    # Verify the files exist
    if not os.path.exists(args.fasta1):
        raise FileNotFoundError(f"Error: file {args.fasta1} does not exist\n")
    if not os.path.exists(args.fasta2):
        raise FileNotFoundError(f"Error: file {args.fasta2} does not exist\n")

    with open(args.fasta1, "r") as f1, open(args.fasta2, "r") as f2:
        while True:
            line = f1.readline()
            if line.strip() == "":
                break

            print(line.strip())
            print(f1.readline().strip())
            print(f2.readline().strip())
            print(f2.readline().strip())


if __name__ == "__main__":
    args = parse_args()
    main(args)
