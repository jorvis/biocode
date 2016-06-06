#!/usr/bin/env python3

"""

Put some general, high-level documentation here

"""

import argparse
import os
import re

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here' )

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    args = parser.parse_args()

    f_reads_correctly_mapped = 0
    f_reads_incorrectly_mapped = 0
    r_reads_correctly_mapped = 0
    r_reads_incorrectly_mapped = 0

    for line in open(args.input_file):
        # lines are like: comp0_c0_seq1   1-T:6   1-F:0   2-T:0   2-F:5
        m = re.search('(.+)\t1-T:(\d+)\t1-F:(\d+)\t2-T:(\d+)\t2-F:(\d+)', line)
        if m:
            f_reads_correctly_mapped += int(m.group(2))
            f_reads_incorrectly_mapped += int(m.group(3))
            r_reads_correctly_mapped += int(m.group(5))
            r_reads_incorrectly_mapped += int(m.group(4))

    correct_count = f_reads_correctly_mapped + r_reads_correctly_mapped
    incorrect_count = f_reads_incorrectly_mapped + r_reads_incorrectly_mapped
    correct_perc = (correct_count / (correct_count + incorrect_count) ) * 100
            
    print("Forward reads correctly mapped: {0}".format(f_reads_correctly_mapped))
    print("Forward reads incorrectly mapped: {0}".format(f_reads_incorrectly_mapped))
    print("Reverse reads correctly mapped: {0}".format(r_reads_correctly_mapped))
    print("Reverse reads incorrectly mapped: {0}".format(r_reads_incorrectly_mapped))
    print("Overall correct percentage: {0}".format(correct_perc))

if __name__ == '__main__':
    main()







