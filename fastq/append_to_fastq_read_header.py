#!/usr/bin/env python3

"""
Running with no options other than -i and -o, this script will transform headers like these:

   @SN7001163:78:C0YG5ACXX:6:1101:1129:2043 1:N:0:CCTAGGT

to headers like this:

   @SN7001163:78:C0YG5ACXX:6:1101:1129:2043/1

That is, it will see the '1' or '2' after the whitespace and transform the header to this instead (in FASTA):

   >SN7001163:78:C0YG5ACXX:6:1101:1129:2043/1

Otherwise, if you want a specific string to be appended to all headers, use the --append option.

Author: Joshua Orvis (jorvis AT gmail)

"""

import argparse
import gzip
import re
import sys

def main():
    parser = argparse.ArgumentParser( description='Append strings to the end of FASTQ read headers')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read (gzip supported)' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    parser.add_argument('-a', '--append', type=str, required=False, help='String to be appended to each read header' )
    args = parser.parse_args()

    if args.output_file is None:
        ofh = sys.stdout
    else:
        ofh = open( args.output_file, 'wt' )

    if args.append is None:
        detect_direction = True
    else:
        detect_direction = False
        
    line_count = 0

    if args.input_file.endswith('.gz'):
        fh = gzip.open( args.input_file, 'rb')
        is_compressed = True
    else:
        fh = open(args.input_file, 'rU')
        is_compressed = False

    for line in fh:
        if is_compressed:
            line = line.decode()
        
        line_count += 1

        if line_count % 4 == 1:

            if detect_direction:
                m = re.search('^(\@.+?) (\d)', line)
                if m:
                    ofh.write("{0}/{1}\n".format(m.group(1), m.group(2)))
                else:
                    raise Exception("ERROR: FASTQ header line found that didn't match expected format: {0}".format(line))
            else:
                m = re.search('^(\@\S+)', line)
                if m:
                    ofh.write("{0}{1}\n".format(m.group(1), args.append))
                else:
                    raise Exception("ERROR: FASTQ header line found that didn't match expected format: {0}".format(line))

            
        else:
            ofh.write(line)
        
    ofh.close()




if __name__ == '__main__':
    main()







