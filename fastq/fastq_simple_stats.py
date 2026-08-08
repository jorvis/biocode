#!/usr/bin/env python3

"""
Reads one or more FASTQ files and reports some basic statistics.  Example:

$ ./fastq_simple_stats.py my.fastq

File: my.fastq
    Total sequence entries: 277971801
    Total bases: 28353123702
    Avg sequence length: 102.0
    GC percentage: 44.5

If you pass multiple files via *.fastq or something similar this will report stats
on these as an aggregate unless you pass the --individual flag, which will cause
reporting on each file individually.

Also handles files if they are GZIP compressed.

Contact: jorvis@gmail.com
"""

import argparse
import codecs
import gzip
import sys
from typing import TextIO


def parse_args() -> argparse.Namespace:
    """
    Simple argument parser for this script.

    :return: argparse.Namespace object containing the input options
    """

    parser = argparse.ArgumentParser(description="Provides simple quantitative statistics for a given FASTQ file")

    parser.add_argument(
        "input_files", metavar="N", type=str, nargs="+", help="Path to one or more input files, separated by spaces"
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=False,
        help="Optional path to an output file to be created, else prints on STDOUT",
    )
    parser.add_argument("--individual", action="store_true", help="Report stats on each file individually")
    parser.add_argument(
        "-p",
        "--progress_interval",
        type=int,
        required=False,
        help="Pass an integer to show progress ever N entries on STDERR",
    )

    return parser.parse_args()


def report_stats(
    fout: TextIO | codecs.StreamWriter, input_file: str, entry_count: int, total_bases: int, gc_count: int
) -> None:

    try:
        avg_entry_length = total_bases / entry_count
    except ZeroDivisionError:
        avg_entry_length = 0.0
    except Exception as e:
        raise Exception(f"Error calculating average entry length: {e}")

    fout.write(f"\nFile: {input_file}\n")
    fout.write(f"\tTotal sequence entries: {entry_count:,}\n")
    fout.write(f"\tTotal bases: {total_bases:,}\n")
    fout.write(f"\tAvg sequence length: {avg_entry_length:,.1f}\n")

    try:
        fout.write(f"\tGC percentage: {((gc_count / total_bases) * 100):.1f}\n")
    except ZeroDivisionError:
        fout.write("\tGC percentage: N/A\n")
    except Exception as e:
        raise Exception(f"Error calculating GC percentage: {e}")


def main(args: argparse.Namespace) -> None:

    # Open the output file, if specified
    if args.output_file is None:
        fout = codecs.getwriter("utf8")(sys.stdout.buffer)
    else:
        fout = open(args.output_file, "w")

    # Values that will be reported
    entry_count: int = 0
    total_bases: int = 0
    gc_count: int = 0

    for input_file in args.input_files:
        print(f"Processing file: {input_file}", file=sys.stderr)
        line_number = 0

        if input_file.endswith(".gz"):
            fh = gzip.open(input_file, "rb")
            is_compressed = True
        else:
            fh = open(input_file, "r")
            is_compressed = False

        if args.individual:
            entry_count = 0
            total_bases = 0
            gc_count = 0

        for line in fh:
            if is_compressed:
                line = line.decode()

            line = line.rstrip()
            line_number += 1

            # Every 4th line is the start of a sequence entry (check still)
            if line_number % 4 == 1:
                if line[0] == "@":
                    entry_count += 1

                    if args.progress_interval:
                        if not entry_count % args.progress_interval:
                            print(f"{entry_count:,} entries processed", file=sys.stderr, end="\r")
                else:
                    raise Exception(f"Line {line_number} expected to start with '@': {line}")
            elif line_number % 4 == 2:
                total_bases += len(line)
                gc_count += line.count("G")
                gc_count += line.count("C")

        if args.individual:
            report_stats(fout, input_file, entry_count, total_bases, gc_count)

    if not args.individual:
        if len(args.input_files) == 1:
            report_stats(fout, args.input_files[0], entry_count, total_bases, gc_count)
        else:
            report_stats(fout, "multiple", entry_count, total_bases, gc_count)


if __name__ == "__main__":
    args = parse_args()
    main(args)
