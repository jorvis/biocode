#!/usr/bin/env python3

"""

Trinity transcripts are named like this:

TRINITY_DN16084_c1_g1_i5
TRINITY_DN16084_c1_g1_i6
TRINITY_DN16084_c1_g1_i7
TRINITY_DN2006_c3_g2_i1
TRINITY_DN2006_c3_g2_i2
TRINITY_DN2006_c3_g2_i3



"""

import argparse
import os
import re

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')
    parser.add_argument('-ia', '--input_abundance', type=str, required=True, help='Input abundance estimation file' )
    parser.add_argument('-tpm', '--tpm_cutoff', type=float, required=True, help='TPM cutoff' )
    #parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    abundance_format = detect_abundance_format(args.input_abundance)

    if abundance_format == 'kallisto':
        process_kallisto(args)
    else:
        raise Exception(
            "ERROR: Couldn't detect abundance file tool source. Currently only Kallisto is supported")

def detect_abundance_format(path):
    with open(path) as f:
        first_line = f.readline()

        if first_line == 'target_id\tlength\teff_length\test_counts\ttpm\n':
            return 'kallisto'
        else:
            return None

def process_kallisto(args):
    abund = read_kallisto_abundances(args.input_abundance)

    abund_group_scores = dict()

    for id in abund:
        m = re.match("(TRINITY_.+_c\d+_g\d+)_i\d+", id)
        if not m:
            raise Exception("The following trinity contig didn't match naming expectations.  See the doc.  id:{0}".format(id))

        base = m.group(1)
        if base in abund_group_scores:
            abund_group_scores[base]['score'] += abund[id]
            abund_group_scores[base]['count'] += 1
        else:
            abund_group_scores[base] = {'score': abund[id], 'count': 1}

    passing_group_count = 0
    passing_transcript_count = 0
            
    for base in abund_group_scores:
        if abund_group_scores[base]['score'] >= args.tpm_cutoff:
            passing_group_count += 1
            passing_transcript_count += abund_group_scores[base]['count']

    print("Aggregate TPM cutoff: {0}".format(args.tpm_cutoff))
    print("Input transcript count: {0}".format(len(abund)))
    print("Groups passing cutoff: {0}/{1}".format(passing_group_count, len(abund_group_scores)))
    print("Transcripts in groups passing cutoff: {0}/{1}".format(passing_transcript_count, len(abund)))

def read_kallisto_abundances(path):
    # key is sequence ID, value is TPM
    abund = dict()
    
    with open(path) as f:
        first_line = f.readline()

        for line in f:
            cols = line.split()
            abund[cols[0]] = float(cols[4].rstrip())

    return abund


if __name__ == '__main__':
    main()







