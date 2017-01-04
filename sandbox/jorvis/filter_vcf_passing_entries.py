#!/usr/bin/env python3.2

'''
Removes entries in a VCF file which don't pass the following conditions:

    - The 5th column must not be '.'
    - The 7th column must be 'PASS'

Comment lines are retained.

Will NOT currently work with merged sample files.
'''


import argparse
import gzip
import re


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-v', '--vcf_file', type=str, required=True, help='Input VCF file' )
    args = parser.parse_args()

    file = args.vcf_file

    fh = None
    if file.endswith(".gz"):
        fh = gzip.GzipFile(file, "r")
    else:
        fh = open(file)

    total_seq_length = 0
    sample_id = None
        
    for line in fh:
        #line = line.decode().rstrip()
        line = line.rstrip()
        cols = line.split("\t")

        if line.startswith("##contig="):
            print(line)
            m = re.match( ".+length\=(\d+)\>", line )
            if m:
                total_seq_length += int(m.group(1))
                continue
            else:
                raise Exception("ERROR: found a contig line with no length defined: {0}".format(line) )

        elif line.startswith("#CHROM"):
            sample_id = cols[9]
            print(line)
            continue
        
        elif line.startswith("#"):
            print(line)
            continue
        
        cols = line.split("\t")
        
        if len(cols) < 7 or cols[4] == '.' or cols[6] != 'PASS':
            continue

        phred_score = get_scaled_phred_quality_score(cols[8], cols[9])

        print(line)
    

    

def get_scaled_phred_quality_score(keys, values):
    keylist = keys.split(':')
    vallist = values.split(':')

    if keylist[3] == 'GQ':
        return float(vallist[3])
    else:
        raise Exception('ERROR: expected GQ in 4th column of this string: {0}'.format(keys))
    



if __name__ == '__main__':
    main()












    
