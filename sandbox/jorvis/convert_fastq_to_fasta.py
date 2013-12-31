#!/usr/bin/env python3

import argparse
#import os
import re

"""
Timing information on a 1.3GB input file:

 0m0.035s - v1. Skeleton with output file creation
 0m1.849s - v2. v1 with loop through input file (pass only)
 0m4.875s - v3. v2 with conditional to check line number
 0m7.617s - v4. v3 writing FASTA output but empty header
 

"""


def main():
    #parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    #parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    #args = parser.parse_args()

    input_file = '/home/jorvis/Dropbox/jhu/metagenomics_410.734/lessons/week03/SRS012902/SRS012902.denovo_duplicates_marked.trimmed.singleton.fastq'
    output_file = '/home/jorvis/Dropbox/jhu/metagenomics_410.734/lessons/week03/SRS012902/SRS012902.denovo_duplicates_marked.trimmed.singleton.fasta'

    ofh = open(output_file,'w')
    line_count = 0
    record_count = 0
    last_header = None

    for line in open(input_file, 'rU'):
        line_count += 1

        if line_count % 4 == 0:
            record_count += 1
            m = re.search('^\@(.+?)\s*', line)
            
        elif line_count % 4 == 2:
            ofh.write(">\n" + line)
        

    ofh.close()


if __name__ == '__main__':
    main()







