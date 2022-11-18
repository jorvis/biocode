#!/usr/bin/env python
from __future__ import division

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
import os
import re

def main():
    parser = argparse.ArgumentParser( description='Mate validation in paired-end FASTQ files')

    ## output file to be written
    parser.add_argument('-l', '--left', type=str, required=True, help='Comma-separated list of left reads in an paired set' )
    parser.add_argument('-r', '--right', type=str, required=True, help='Comma-separated list of right reads in an paired set' )
    args = parser.parse_args()

    check_input_options(args)
    left_files = args.left.split(',')
    right_files = args.right.split(',')

    lindex = 0
    
    for lpath in left_files:
        if lpath.endswith('.gz'):
            ifh = gzip.open(lpath, 'rb')
            left_is_compressed = True

            if args.right is not None:
                ifh_right = gzip.open(right_files[lindex], 'rb')
        else:
            ifh = open(lpath, 'rU')
            left_is_compressed = False
            if args.right is not None:
                ifh_right = open(right_files[lindex], 'rU')

        while True:
            line1 = ifh.readline().replace("\n", '')
            print("line1: ({0})".format(line1))
            rline1 = ifh_right.readline()

            # Check if either file ended early
            if not line1:
                if not rline1:
                    # expected, matching end of both files
                    break
                else:
                    raise Exception("ERROR: Left file ended before right file")

            if not rline1:
                raise Exception("ERROR: Right file ended before left file")

            # We don't actually need these for now.  Could just read to advance the cursor
            line2 = ifh.readline()
            line3 = ifh.readline()
            line4 = ifh.readline()
            rline2 = ifh_right.readline()
            rline3 = ifh_right.readline()
            rline4 = ifh_right.readline()

            print(line1)

            if len(line1):
                print("line1 is ({0})".format(line1))
                lread_base = line1.split()[0]
            else:
                lread_base = False
            
            if len(rline1):
                rread_base = rline1.split()[0]
            else:
                rread_base = False

            ## did we reach the natural ending?
            if not lread_base and not rread_base:
                break
            
            ## Let's first check the easy case that they match:
            if lread_base == rread_base:
                pass
            else:
                if lread_base[-2:] == '/1' and rread_base[-2:] == '/2':
                    # compare without the mate info then
                    if lread_base[:-2] != rread_base[:-2]:
                        raise Exception("Mismatch found: left:{0} right:{1}".format(lread_base, rread_base))
                else:
                    raise Exception("Mismatch found: left:{0} right:{1}".format(lread_base, rread_base))

            ## Now those with file names with ids ending like /1 and /2 to denote mate pairs
            print("Comparing {0} with {1}".format(lread_base, rread_base))

    lindex += 1


def check_input_options(args):
    left_files = args.left.split(',')
    right_files = args.right.split(',')
        
    if len(left_files) != len(right_files):
        raise Exception("ERROR: Count of left files ({0}) must be the same as right files ({1})".format(
            len(left_files), len(right_files)))
        
if __name__ == '__main__':
    main()







