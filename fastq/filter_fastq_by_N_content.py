#!/usr/bin/env python3

"""

Assumptions:
- The sequence residue lines are in upper case.  I've never seen FASTQ files otherwise,
  but if they aren't this will report no reads being filtered.  A check would slow the
  script too much.

Author: Joshua Orvis (jorvis AT gmail)

"""

import argparse
from typing import TextIO


def parse_arguments() -> argparse.Namespace:
    """Simple argument parser for the script."""
    parser = argparse.ArgumentParser(description="Filter FASTQ file by N content")

    parser.add_argument("-l", "--left", type=str, required=False, help="FASTQ: Left (R1) mate pair file")
    parser.add_a("-r", "--right", type=str, required=False, help="FASTQ: Right (R2) mate pair file")
    parser.add_argument("-s", "--singletons", type=str, required=False, help="FASTQ: singleton reads file")
    parser.add_argument("-lo", "--left_out", type=str, required=False, help="Output file of Left reads kept")
    parser.add_argument("-ro", "--right_out", type=str, required=False, help="Output file of Right reads kept")
    parser.add_argument("-so", "--singletons_out", type=str, required=False, help="Output file of singleton reads kept")
    parser.add_argument(
        "-p",
        "--percent_n_cutoff",
        type=int,
        options=range(1, 101),
        required=True,
        help="1-100: Percentage of Ns in read >= this will cause removal",
    )
    parser.add_argument(
        "-or", "--output_report", type=str, required=False, help="Optional report of read filtering statistics"
    )

    return parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> int:
    """
    Validate the arguments passed to the script. For cleaner code, the input options are mapped to a binary state.

    One of the following sets must be passed:
    - [left, right]
    - [left, right, singletons]
    - [singletons]
    """

    left = 0 if args.left is not None else 1
    right = 0 if args.right is not None else 2
    singletons = 0 if args.singletons is not None else 4
    arg_state = left + right + singletons

    # Valid states are 3 (1+2), 7 (1+2+4), and 4
    match arg_state:
        case 0:
            raise Exception("ERROR: You must pass at least one of --left, --right, or --singletons")
        case 1, 5:
            raise Exception("ERROR: If you specify --left you must also specify --right")
        case 2, 6:
            raise Exception("ERROR: If you specify --right you must also specify --left")
        case 3:
            assert args.left_out is not None, "ERROR: If you specify --left you must also specify --left_out"
            assert args.right_out is not None, "ERROR: If you specify --right you must also specify --right_out"
            assert (
                args.singletons_out is not None
            ), """ERROR: If you pass --left and --right you must still pass --singletons_out.
                This is so orphaned mate pairs can still be exported."""
            return arg_state
        case 4:
            assert (
                args.singletons_out is not None
            ), "ERROR: If you pass --singletons you must also pass --singletons_out"
            if args.left_out is not None or args.right_out is not None:
                print("--left and --right arguments not passed, --left_out and --right_out options will be ignored")
            return arg_state
        case 7:
            assert args.left_out is not None, "ERROR: If you specify --left you must also specify --left_out"
            assert args.right_out is not None, "ERROR: If you specify --right you must also specify --right_out"
            assert (
                args.singletons_out is not None
            ), "ERROR: If you pass --singletons you must also pass --singletons_out"
            return arg_state
        case _:
            raise Exception("ERROR: Invalid argument state")


def left_right_processing(
    counts: dict,
    pct_cutoff: int,
    l_in_fh: TextIO,
    r_in_fh: TextIO,
    l_out_fh: TextIO,
    r_out_fh: TextIO,
    s_out_fh: TextIO,
    line_limit=None,
) -> dict:
    """Process left and right reads."""

    l_keep = False
    r_keep = False
    last_lheader = None
    last_rheader = None
    line_count = 0

    for l_line, r_line in zip(l_in_fh, r_in_fh):
        line_count += 1

        if line_limit is not None and line_count > line_limit:
            print("WARNING: Program exited after {0} lines for debugging".format(line_limit))
            break

        if line_count % 4 == 1:
            last_lheader = l_line
            last_rheader = r_line
            l_keep = False
            r_keep = False

        elif line_count % 4 == 2:
            counts["total_reads"] += 2
            l_pct_n = (l_line.count("N") / (len(l_line) - 1)) * 100
            r_pct_n = (r_line.count("N") / (len(r_line) - 1)) * 100

            # Left bad, right good
            if l_pct_n >= pct_cutoff and r_pct_n < pct_cutoff:
                l_keep = False
                r_keep = s_out_fh
                r_keep.write(last_rheader)
                r_keep.write(r_line)
                counts["left_only_discarded"] += 1
                counts["right_only_kept"] += 1

            # Left good, right bad
            elif l_pct_n < pct_cutoff and r_pct_n >= pct_cutoff:
                r_keep = False
                l_keep = s_out_fh
                l_keep.write(last_lheader)
                l_keep.write(l_line)
                counts["left_only_kept"] += 1
                counts["right_only_discarded"] += 1

            # Both good
            elif l_pct_n < pct_cutoff and r_pct_n < pct_cutoff:
                l_keep = l_out_fh
                r_keep = r_out_fh
                l_keep.write(last_lheader)
                l_keep.write(l_line)
                r_keep.write(last_rheader)
                r_keep.write(r_line)
                counts["pairs_kept"] += 1

            # Both bad
            else:
                counts["pairs_discarded"] += 1

        # Handle the third/fourth lines
        else:
            if r_keep:
                r_keep.write(r_line)

            if l_keep:
                l_keep.write(l_line)

    return counts


def singleton_processing(counts: dict, pct_cutoff: int, s_in_fh: TextIO, s_out_fh: TextIO) -> dict:
    """Process singleton reads."""

    if args.singletons is not None:
        line_count = 0
        last_header = None
        keep = False

        for line in s_in_fh:
            line_count += 1

            if line_count % 4 == 1:
                last_header = line
                keep = False

            elif line_count % 4 == 2:
                counts["total_reads"] += 1
                pct_n = (line.count("N") / (len(line) - 1)) * 100

                if pct_n < pct_cutoff:
                    keep = True
                    s_out_fh.write(last_header)
                    s_out_fh.write(line)
                    counts["singletons_kept"] += 1
                else:
                    counts["singletons_discarded"] += 1
            else:
                if keep:
                    s_out_fh.write(line)

    return counts


def main(args) -> None:

    # Debugging on large files only
    line_limit = None  # noqa

    # Determine the set of inputs and outputs
    run_state = validate_arguments(args)

    counts = {
        "lines": 0,
        "records": 0,
        "total_reads": 0,
        "pairs_kept": 0,
        "pairs_discarded": 0,
        "left_only_kept": 0,
        "left_only_discarded": 0,
        "right_only_kept": 0,
        "right_only_discarded": 0,
        "singletons_kept": 0,
        "singletons_discarded": 0,
    }

    # Left and right provided
    if run_state == 3 or run_state == 7:
        with (
            open(args.left, "r") as l_in_fh,
            open(args.right, "r") as r_in_fh,
            open(args.left_out, "wt") as l_out_fh,
            open(args.right_out, "wt") as r_out_fh,
            open(args.singletons_out, "wt") as s_out_fh,
        ):
            counts = left_right_processing(
                counts, args.percent_n_cutoff, l_in_fh, r_in_fh, l_out_fh, r_out_fh, s_out_fh
            )

    # Singletons provided
    if run_state == 4 or run_state == 7:
        with open(args.singletons, "r") as s_in_fh, open(args.singletons_out, "wt") as s_out_fh:
            counts = singleton_processing(counts, args.percent_n_cutoff, s_in_fh, s_out_fh)

    # Output the report
    if args.output_report is not None:
        with open(args.output_report, "wt") as report_fh:
            report_fh.write("Total input reads: {0}\n".format(counts["total_reads"]))
            report_fh.write("Read pairs kept: {0}\n".format(counts["pairs_kept"]))
            report_fh.write("Read pairs discarded: {0}\n".format(counts["pairs_discarded"]))
            report_fh.write("Left only reads kept: {0}\n".format(counts["left_only_kept"]))
            report_fh.write("Left only reads discarded: {0}\n".format(counts["left_only_discarded"]))
            report_fh.write("Right only reads kept: {0}\n".format(counts["right_only_kept"]))
            report_fh.write("Right only reads discarded: {0}\n".format(counts["right_only_discarded"]))
            report_fh.write("Singleton reads kept: {0}\n".format(counts["singletons_kept"]))
            report_fh.write("Singleton reads discarded: {0}\n".format(counts["singletons_discarded"]))


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
