#!/usr/bin/env python3

"""
A memory-efficient script for validating mates in paired-end FASTQ files.

This actually came to be because of a help request thread in /r/bioinformatics:

https://www.reddit.com/r/bioinformatics/comments/owed06/a_tool_to_verify_paired_end_files/

Tested against headers conventions like this:

@SN7001163:145:C36LWACXX:4:1101:1178:2099 1:N:0:TTAGGCA
@SN7001163:145:C36LWACXX:4:1101:1178:2099 2:N:0:TTAGGCA

@61JCNAAXX100503:5:100:10001:18267/1
@61JCNAAXX100503:5:100:10001:18267/2

Assumptions:

Input reads can be in raw FASTQ format or compressed with a .gz extension.  If
passing paired-end data, all pairwise entries must be in the same order within
their files (this is very common.)

Author:
Joshua Orvis
"""

import argparse
import gzip


def parse_args() -> argparse.Namespace:
    """
    Simple argument parser for this script.

    :return: argparse.Namespace object containing the input options
    """

    parser = argparse.ArgumentParser(description="Mate validation in paired-end FASTQ files")
    parser.add_argument(
        "-l", "--left", type=str, required=True, help="Comma-separated list of left reads in a paired set"
    )
    parser.add_argument(
        "-r", "--right", type=str, required=True, help="Comma-separated list of right reads in a paired set"
    )

    return parser.parse_args()


def compare_reads(ifh, ifh_right) -> None:
    """
    Parse through two read files and compare the read headers.

    :param ifh: file handle for read 1
    :param ifh_right: file handle for read 2

    :return: None
    """

    # Step through each file and compare the read headers
    try:

        line_type = 1
        first_line = True
        read_style = None

        for r1, r2 in zip(ifh, ifh_right, strict=True):

            # For the first entry, determine the read pair style (i.e. /1 or 1:N:0:TTAGGCA)
            if first_line:
                r1 = r1.decode().strip()
                r2 = r2.decode().strip()
                if r1[-2] == "/" and r2[-2] == "/":
                    read_style = 1
                elif r1.split(" ")[1].startswith("1") and r2.split(" ")[1].startswith("2"):
                    read_style = 2
                else:
                    raise Exception(f"ERROR: Read style not recognized! ({r1}, {r2})")
                first_line = False

            # We are only interested in the ID lines
            if line_type == 1:
                r1 = r1.decode().strip()
                r2 = r2.decode().strip()

                # Compare the read headers
                if read_style == 1:
                    if r1[:-2] != r2[:-2]:
                        raise Exception(f"ERROR: Headers do not match! ({r1}, {r2})")
                elif read_style == 2:
                    if r1.split(" ")[0] != r2.split(" ")[0]:
                        raise Exception(f"ERROR: Headers do not match! ({r1}, {r2})")

            # Reset the line type if we reach the end of the read
            if line_type == 4:
                line_type = 1

            line_type += 1

    except ValueError:
        # If the files are not the same length, we will get a ValueError
        raise Exception("ERROR: Files are not the same length!")


def main(args: argparse.Namespace) -> None:
    """
    Main function for the script.

    :param args: argparse.Namespace object containing the input options
    :return: None
    """

    left_files = args.left.split(",")
    right_files = args.right.split(",")

    # Validate that there are equal numbers of R1 and R2 files
    if len(left_files) != len(right_files):
        raise Exception(
            f"ERROR: Count of left files ({len(left_files)}) must be the same as right files ({len(right_files)})"
        )

    # Build a list of tuples with the read pairs and their gzip status
    pairs = []
    for lfile, rfile in zip(left_files, right_files):
        if lfile.endswith(".gz") and rfile.endswith(".gz"):
            gz_status = True
        elif not lfile.endswith(".gz") and not rfile.endswith(".gz"):
            gz_status = False
        else:
            raise Exception(f"ERROR: gzip status must be the same for both files! ({lfile}, {rfile})")

        pairs.append((lfile, rfile, gz_status))

    # Iterate over the files and compare the reads
    for lfile, rfile, gz_status in pairs:
        if gz_status:
            with gzip.open(lfile, "rb") as ifh, gzip.open(rfile, "rb") as ifh_right:
                compare_reads(ifh, ifh_right)
        else:
            with open(lfile, "r") as ifh, open(rfile, "r") as ifh_right:
                compare_reads(ifh, ifh_right)


if __name__ == "__main__":
    args = parse_args()
    main(args)
