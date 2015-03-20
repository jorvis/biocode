#!/usr/bin/env python3

"""

Assumptions:
- The sequence residue lines are in upper case.  I've never seen FASTQ files otherwise,
  but if they aren't this will report no reads being filtered.  A check would slow the
  script too much.

Author: Joshua Orvis (jorvis AT gmail)

"""

import argparse
import re
import sys

def main():
    parser = argparse.ArgumentParser( description='Filter FASTQ file by N content')

    ## output file to be written
    parser.add_argument('-l', '--left', type=str, required=False, help='FASTQ: Left (R1) mate pair file' )
    parser.add_argument('-r', '--right', type=str, required=False, help='FASTQ: Right (R2) mate pair file' )
    parser.add_argument('-s', '--singletons', type=str, required=False, help='FASTQ: singleton reads file' )
    parser.add_argument('-lo', '--left_out', type=str, required=False, help='Output file of Left reads kept' )
    parser.add_argument('-ro', '--right_out', type=str, required=False, help='Output file of Right reads kept')
    parser.add_argument('-so', '--singletons_out', type=str, required=False, help='Output file of singleton reads kept')
    parser.add_argument('-p', '--percent_n_cutoff', type=int, required=True, help='1-100: Percentage of Ns in read >= this will cause removal')
    parser.add_argument('-or', '--output_report', type=str, required=False, help='Optional report of read filtering statistics')

    args = parser.parse_args()

    l_in_fh = None
    r_in_fh = None
    s_in_fh = None

    l_out_fh = None
    r_out_fh = None
    s_out_fh = None

    # Let the option checking and file handle bonanza begin!

    # One of these must be passed: [left, right], [left, right, singletons], [singletons]
    if args.left is not None and args.right is not None:
        # singletons is optional if L/R is passed
        if args.singletons is not None:
            s_in_fh = open(args.singletons)

        if args.left_out is None:
            raise Exception("ERROR: If you specify --left you must also specify --left_out")
        else:
            l_out_fh = open(args.left_out, 'wt')

        if args.right_out is None:
            raise Exception("ERROR: If you specify --right you must also specify --right_out")
        else:
            r_out_fh = open(args.right_out, 'wt')

        if args.singletons_out is None:
            raise Exception("ERROR: The --singletons_out option must be passed if --left and --right are specified, even " \
                            "if --singletons isn't.  This is so that orphaned mate pairs can still be exported.")
        else:
            s_out_fh = open(args.singletons_out, 'wt')
        
    elif args.singletons is not None:
        s_in_fh = open(args.singletons)

        if args.singletons_out is None:
            raise Exception("ERROR: If you pass --singletons you must also pass --singletons_out")
        else:
            s_out_fh = open(args.singletons_out, 'wt')
        
    else:
        raise Exception("ERROR.  One of these must be passed: [left, right], [left, right, singletons], [singletons]")
        
    
    ##############################################################
    #  now the actual functioning code
    ##############################################################
    line_count = 0
    record_count = 0
    last_lheader = None
    last_rheader = None

    #counts = {'left_kept':0, 'left_discarded':0, 'right_kept':0, 'right_discarded':0, 'singlets_kept':0, 'singlets_discarded':0}
    counts = {'total_reads':0, 'pairs_kept':0, 'pairs_discarded':0, \
              'left_only_kept':0, 'left_only_discarded':0, \
              'right_only_kept':0, 'right_only_discarded':0, \
              'singletons_kept':0, 'singletons_discarded':0
    }

    # for debugging purposes on large files only
    line_limit = None

    if args.left is not None and args.right is not None:
        l_keep = False
        r_keep = False
        
        # Read left and right files at once.
        with open(args.left) as l_in_fh, open(args.right) as r_in_fh:
            for l_line, r_line in zip(l_in_fh, r_in_fh):
                line_count += 1

                if line_limit is not None and line_count > line_limit:
                    print("WARNING: Program exited after {0} lines for debugging".format(line_limit))
                    break

                if line_count % 4 == 1:
                    record_count += 1
                    last_lheader = l_line
                    last_rheader = r_line
                    l_keep = False
                    r_keep = False
                    
                elif line_count % 4 == 2:
                    counts['total_reads'] += 2
                    l_pct_n = (l_line.count('N') / (len(l_line) - 1)) * 100
                    r_pct_n = (r_line.count('N') / (len(r_line) - 1)) * 100

                    #print("LN%: {0}\tRN%{1}: ".format(l_pct_n, r_pct_n))

                    # left bad, right good
                    if l_pct_n >= args.percent_n_cutoff and r_pct_n < args.percent_n_cutoff:
                        l_keep = False
                        r_keep = s_out_fh
                        r_keep.write(last_rheader)
                        r_keep.write(r_line)
                        counts['left_only_discarded'] += 1
                        counts['right_only_kept'] += 1

                    # left good, right bad
                    elif l_pct_n < args.percent_n_cutoff and r_pct_n >= args.percent_n_cutoff:
                        r_keep = False
                        l_keep = s_out_fh
                        l_keep.write(last_lheader)
                        l_keep.write(l_line)
                        counts['left_only_kept'] += 1
                        counts['right_only_discarded'] += 1

                    # both good
                    elif l_pct_n < args.percent_n_cutoff and r_pct_n < args.percent_n_cutoff:
                        l_keep = l_out_fh
                        r_keep = r_out_fh
                        l_keep.write(last_lheader)
                        l_keep.write(l_line)
                        r_keep.write(last_rheader)
                        r_keep.write(r_line)
                        counts['pairs_kept'] += 1

                    # both bad
                    else:
                        counts['pairs_discarded'] += 1
                else:
                    # handle the third/fourth lines
                    if r_keep:
                        r_keep.write(r_line)

                    if l_keep:
                        l_keep.write(l_line)
                

    if args.singletons is not None:
        # s_in_fh
        # s_out_fh
        line_count = 0
        record_count = 0
        last_header = None
        keep = False
        
        for line in s_in_fh:
            line_count += 1

            if line_count % 4 == 1:
                record_count += 1
                last_header = line
                keep = False

            elif line_count % 4 == 2:
                counts['total_reads'] += 1
                pct_n = (line.count('N') / (len(line) - 1)) * 100

                if pct_n < args.percent_n_cutoff:
                    keep = True
                    s_out_fh.write(last_header)
                    s_out_fh.write(line)
                    counts['singletons_kept'] += 1
                else:
                    counts['singletons_discarded'] += 1
            else:
                if keep:
                    s_out_fh.write(line)
        
    
    if args.output_report is not None:
        report_fh = open(args.output_report, 'wt')

        report_fh.write("Total input reads: {0}\n".format(counts['total_reads']))
        report_fh.write("Read pairs kept: {0}\n".format(counts['pairs_kept']))
        report_fh.write("Read pairs discarded: {0}\n".format(counts['pairs_discarded']))
        report_fh.write("Left only reads kept: {0}\n".format(counts['left_only_kept']))
        report_fh.write("Left only reads discarded: {0}\n".format(counts['left_only_discarded']))
        report_fh.write("Right only reads kept: {0}\n".format(counts['right_only_kept']))
        report_fh.write("Right only reads discarded: {0}\n".format(counts['right_only_discarded']))
        report_fh.write("Singleton reads kept: {0}\n".format(counts['singletons_kept']))
        report_fh.write("Singleton reads discarded: {0}\n".format(counts['singletons_discarded']))


if __name__ == '__main__':
    main()







