#!/usr/bin/env python3

"""
Initially written to deal with a specific issue where a user has two different
fastq files 'left.fq' and 'right.fq'.  In each the directionality (/1 or /2) has been
stripped from the ID in the read headers (but it will still work if they're present.)
Also, the reads pair members are not in pairwise order in the two files as they should be.

This script takes these two files, adds the direction back, then orders them.  Any
singletons are written out to a separate file.

The ID hashing requirements of the left read file dictate the amount of memory
required to run this, which could be considerable.  Attempts are made to minimize
this by iterating both left and right files at once, deleting keys as pairs are
found.

WARNING: For speed, this script completely trusts that the first record starts
on the first line, that there are no comments or whitespace lines, and that each
record is 4 lines in length.

WARNING AGAIN: Also, the script currently assumes that the record header lines
contain on the @ symbol followed by the identifier, with no annotations afterward,
like this:

    @SN180:361:C2G2WACXX:7:1101:20475:1998
    or
    @SN180:361:C2G2WACXX:7:1101:20475:1998/1
    or
    @SN180:361:C2G2WACXX:7:1101:20475:1998/2

"""

import argparse
import fileinput
import re

# key = id
# value = list where [0]=seq [1]=label [2]=qual
total_reads_read = 0
left_reads = dict()
left_reads_written = 0
right_reads = dict()
right_reads_written = 0
singletons_written = 0

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-il', '--input_left', type=str, required=True, help='Input left reads file' )
    parser.add_argument('-ir', '--input_right', type=str, required=True, help='Input right reads file' )
    parser.add_argument('-ol', '--output_left', type=str, required=True, help='Output left reads file' )
    parser.add_argument('-or', '--output_right', type=str, required=True, help='Output right reads file' )
    parser.add_argument('-os', '--output_singleton', type=str, required=True, help='Output singletons reads file' )
    args = parser.parse_args()

    line_count = 0
    current_left_header = None
    current_left_cols = list()

    current_right_header = None
    current_left_cols = list()

    global total_reads_read
    global left_reads
    global right_reads

    ## open input/output files
    lifh = open(args.input_left)
    rifh = open(args.input_right)
    lofh = open(args.output_left, 'wt')
    rofh = open(args.output_right, 'wt')
    sofh = open(args.output_singleton, 'wt')

    for left_line, right_line in zip(lifh, rifh):
        line_count += 1

        # id line
        if line_count % 4 == 1:
            total_reads_read += 2
            # process the existing entries
            if current_left_header is not None:
                process_current_entries( current_left_header, current_left_cols, current_right_header, \
                                         current_right_cols, lofh, rofh, )

            current_left_header = left_line.lstrip('@').rstrip()
            current_left_cols = list()

            m = re.match("(.+)\/1", current_left_header)
            if m:
                current_left_header = m.group(1)

            current_right_header = right_line.lstrip('@').rstrip()
            current_right_cols = list()

            m = re.match("(.+)\/1", current_right_header)
            if m:
                current_right_header = m.group(1)
        else:
            current_left_cols.append(left_line)
            current_right_cols.append(right_line)

    if current_left_header is not None:
        process_current_entries( current_left_header, current_left_cols, current_right_header, \
                                 current_right_cols, lofh, rofh)

    # Second-pass after parsing is complete
    for lid in left_reads:
        if lid in right_reads:
            export_reads(lid, left_reads[lid], lid, right_reads[lid], lofh, rofh)
            del right_reads[lid]
        else:
            export_singleton(sofh, lid, left_reads[lid], '/1')

    # if we get here they must all be singletons
    for rid in right_reads:
        export_singleton(sofh, rid, right_reads[rid], '/2')

    print("DEBUG: total reads read: {0}".format(total_reads_read) )
    print("DEBUG: left reads written: {0}".format(left_reads_written) )
    print("DEBUG: right reads written: {0}".format(right_reads_written) )
    print("DEBUG: singletons written: {0}".format(singletons_written) )

def export_reads(lheader, lcols, rheader, rcols, lofh, rofh):
    global left_reads_written
    global right_reads_written

    lofh.write("@{0}/1\n".format(lheader))
    lofh.write("".join(lcols))
    left_reads_written += 1
    
    rofh.write("@{0}/2\n".format(rheader))
    rofh.write("".join(rcols))
    right_reads_written += 1

def export_singleton(ofh, id, cols, direction):
    global singletons_written

    ofh.write("@{0}{1}\n".format(id, direction))
    ofh.write("".join(cols))
    singletons_written += 1


def process_current_entries(lheader, lcols, rheader, rcols, lofh, rofh):
    global left_reads
    global right_reads

    left_reads[lheader] = lcols
    right_reads[rheader] = rcols

    if lheader in right_reads:
        export_reads(lheader, lcols, lheader, right_reads[lheader], lofh, rofh)
        del left_reads[lheader]
        del right_reads[lheader]
        
    elif rheader in left_reads:
        export_reads(rheader, left_reads[rheader], rheader, rcols, lofh, rofh)
        del left_reads[rheader]
        del right_reads[rheader]
        



if __name__ == '__main__':
    main()







