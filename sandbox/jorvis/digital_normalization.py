#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict

"""

You probably don't want to run this.

I've been using two different digital normalization tools for a while now and am only
writing this rudimentary implementation to get a better feel for the algorithm works.

If you need digital normalization code go here:

https://github.com/ged-lab/khmer

And yes, I know this could be written more succinctly to the point of obfuscation, but
I'm using it to teach the method.  It works.

"""

def main():
    parser = argparse.ArgumentParser( description='Rudimentary digital normalization implementation')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input FASTQ file' )
    parser.add_argument('-k', '--kmer_size', type=int, required=False, default=20, help='K-mer size to use (default=20)' )
    parser.add_argument('-C', '--coverage_target', type=int, required=False, default=20, help='Default coverage target over which reads are discarded' )
    args = parser.parse_args()

    global_kmer_counts = defaultdict(int)
    current_line_number = 1
    reads_kept = 0
    reads_total = 0
    
    for seq in open(args.input_file):
        if current_line_number % 4 == 2:
            seq = seq.rstrip()
            reads_total += 1
            read_counts = defaultdict(int)

            for i in range(len(seq) - args.kmer_size):
                kmer = seq[i:i+args.kmer_size]
                if kmer in global_kmer_counts:
                    read_counts[kmer] += global_kmer_counts[kmer]
                else:
                    read_counts[kmer] += 1

            m = median(read_counts.values())

            if m < args.coverage_target:
                reads_kept += 1
                print("Keeping read {2}/{3} (med: {1})".format(seq, m, reads_kept, reads_total))
                for kmer in read_counts:
                    global_kmer_counts[kmer] += read_counts[kmer]

        current_line_number += 1

    print("Completed.  {0} reads out of {1} would have been kept.".format(reads_kept, reads_total) )
    
def median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if not length % 2:
        return (sorts[int(length / 2)] + sorts[int(length / 2 - 1)]) / 2.0
    return sorts[int(length / 2)]


if __name__ == '__main__':
    main()







