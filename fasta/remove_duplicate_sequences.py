#!/usr/bin/env python3

'''
Description:

Reads a multi-FASTA file sequence and removes duplicates based on MD5sums
of the sequence strings.  Of course, this means there's a possibility for
hash collisions but for the purposes it was written for this was
acceptable.  Use at your own discretion.

Input: 

- A (multi-)FASTA file.  

Output:

Output will be multi-FASTA data printed to STDOUT or a file if you
pass the -o option.  The output sequences are written with 60 characters
per line.

'''

import argparse
import hashlib
import sys

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Read a multi-FASTA file sequence and remove duplicates (by MD5 hash)')

    ## output file to be written
    parser.add_argument('-i', '--fasta_file', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-o', '--output_file', type=str, required=False, default=None, help='Optional Path to an output file to be created' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output_file is not None:
        fout = open(args.output_file, 'wt')
    
    seqs = utils.fasta_dict_from_file(args.fasta_file)
    found = dict()

    m = hashlib.md5()

    for seq_id in seqs:
        seq = seqs[seq_id]['s']
        m.update(seq.encode())
        md5sum = m.hexdigest()
        
        ## write this sequence, 60bp per line
        if md5sum not in found:
            fout.write(">{0} {1}\n".format(seq_id, seqs[seq_id]['h']))
            for i in range(0, len(seq), 60):
                fout.write(seq[i : i + 60] + "\n")
            found[md5sum] = 1

if __name__ == '__main__':
    main()







