#!/usr/bin/env python3

"""

Read BAM file, split each line into columns, build single NON-nested dict of transcript ids.
Also parses CIGAR line to compute coordinates and parses the read name AND stores each 
alignment in the list() within the transcript id-keyed dict.
  
"""

import argparse
import os
import re
import sys


def main():
    parser = argparse.ArgumentParser( description='Groups transcripts by mapped read pairs')
    parser.add_argument('-i', '--input_sam_file', type=str, required=False, help='SAM file of read alignments back to transcripts' )
    args = parser.parse_args()

    single_read_pairings = 0
    same_read_pairings = 0
    pairings = dict()

    print("INFO: parsing SAM file and creating transcript pairings")
    fh = sys.stdin
    if args.input_sam_file != None:
        fh = open(args.input_sam_file)
    n_lines = 0

    for line in fh:
        n_lines = n_lines + 1
        if line[0] == '@': continue

        cols = line.split("\t")
        ref_read = cols[0]
        ref_transcript = cols[2]
        other_transcript = cols[6]

        if (n_lines % 1000000) == 0:
            print("n_lines=" + str(n_lines))

        # we don't care about the lines where the mate is unmapped or mapped to the same
        if other_transcript == '*':
            single_read_pairings += 1
            continue
        elif other_transcript == '=':
            same_read_pairings += 1
            continue

        m = re.match("(.+)__[12]$", ref_read)
        if m:
            ref_read_base = m.group(1)
        else:
            m = re.match("(.+)\/[12]$", ref_read)
            if m:
                ref_read_base = m.group(1)
            else:
                ref_read_base = ref_read

        transcripts = sorted([ref_transcript, other_transcript])
        cigar = cols[5]

        # query
        qstart = 1
        m = re.match("^(\d+)[SH]", cigar)
        if m:
            qstart += int(m.group(1))

        qlen = 0
        for m in re.finditer("(\d+)[M=XI]", cigar):
            qlen += int(m.group(1))

        qend = qstart + qlen - 1

        # reference
        rstart = int(cols[3])
        rlen = 0
        for m in re.finditer("(\d+)[M=XDN]", cigar):
            rlen += int(m.group(1))

        rend = rstart + rlen - 1
        
        if transcripts[0] not in pairings:
            pairings[transcripts[0]] = dict()
            
        if transcripts[1] not in pairings[transcripts[0]]:
            pairings[transcripts[0]][transcripts[1]] = list()

        if ref_read_base not in pairings[transcripts[0]][transcripts[1]]:
            pairings[transcripts[0]][transcripts[1]].append( {
                'id': ref_read_base, 'qstart':qstart, 'qend':qend, 'rstart':rstart, 'rend':rend
            } )

    print ("read " + str(n_lines) + " line(s) from file, " + str(len(pairings)) + " pairings")
    print ("single read pairings=" + str(single_read_pairings) + " same_read_pairings=" + str(same_read_pairings))

if __name__ == '__main__':
    main()







