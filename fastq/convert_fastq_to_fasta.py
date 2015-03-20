#!/usr/bin/env python3

"""
This is a minimal, speed-optimized version of a FASTQ->FASTA conversion script.  It lacks a
lot of options other scripts might have, but I got tired of runtimes taking hours with
billion-read FASTQ files.

Example timing information on a 5.4GB input file (with 21.7 million records)**:

    - 55.3 seconds with default settings
    - 2 minutes,  2 seconds when run with --detect_direction option
    - 2 minutes,  8 seconds when run with the --width=60 option
    - 3 minutes, 17 seconds when run with both --detect_direction and --width=60 options

** Timing was run on a laptop with an Intel(R) Core(TM) i7-3740QM CPU @ 2.70GHz, 
utilizing only one core.

OUTPUT
------

Running with no options other than -i and -o, this script will transform headers like these:

   @SN7001163:78:C0YG5ACXX:6:1101:1129:2043 1:N:0:CCTAGGT
   @61JCNAAXX100503:5:100:10001:18267/2

to this in the FASTA output:   
   
   >SN7001163:78:C0YG5ACXX:6:1101:1129:2043 1:N:0:CCTAGGT
   >61JCNAAXX100503:5:100:10001:18267/2

If you have have headers like this:

   @SN7001163:78:C0YG5ACXX:6:1101:1129:2043 1:N:0:CCTAGGT

The optional --detect_direction option will see the '1' or '2' after the whitespace and transform
the header to this instead (in FASTA):

   >SN7001163:78:C0YG5ACXX:6:1101:1129:2043/1

Note that this can increase the runtime up to 3x (from my testing), since a regex is used.

The FASTQ format dictates that the sequence residues of each entry be all contained within a
single line, while the convention for FASTA format is a fixed-with for all residue lines (usually
60bp) where longer ones are wrapped to the next line.  By default this script will simply copy
the residue lines over from the FASTQ to the FASTA file, but you can use the --width option
to specify that the lines in the output FASTA should be wrapped at some width.  This comes with
a performance penalty.

Author: Joshua Orvis (jorvis AT gmail)

"""

import argparse
import re
import sys

def main():
    parser = argparse.ArgumentParser( description='Convert FASTQ to FASTA format')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    parser.add_argument('-d', '--detect_direction', action='store_true', help='Pass this flag to auto-detect the mate pair direction.  See full documentation for more info' )
    parser.add_argument('-w', '--width', type=int, required=False, help='Defines the width (or number of bases) on each FASTA sequence line' )
    args = parser.parse_args()

    if args.output_file is None:
        ofh = sys.stdout
    else:
        ofh = open( args.output_file, 'wt' )
    
    line_count = 0
    record_count = 0
    last_header = None

    for line in open(args.input_file, 'rU'):
        line_count += 1

        if line_count % 4 == 1:
            record_count += 1

            if args.detect_direction:
                m = re.search('^\@(.+?) (\d)', line)
                if m:
                    last_header = "{0}/{1}\n".format(m.group(1), m.group(2))
                else:
                    raise Exception("ERROR: FASTQ header line found that didn't match expected format: {0}".format(line))
            else:
                line = line.lstrip('@')
                last_header = line
                
        elif line_count % 4 == 2:
            if args.width:
                ofh.write(">{0}{1}\n".format(last_header, wrapped(line, args.width)))
            else:
                ofh.write(">{0}{1}".format(last_header, line))
        
    ofh.close()

    print("{0} records written to the output FASTA file".format(record_count))


def wrapped(string, every=60):
    string = string.rstrip()
    ''' this runs orders of magnitude faster than using the textwrap module '''
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))


if __name__ == '__main__':
    main()







