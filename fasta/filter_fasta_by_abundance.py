#!/usr/bin/env python3

"""

Use case example:  You've assembled a transcriptome using Trinity and ran an abundance
estimator like Kallisto.  Now you want to filter the transcript FASTA file to only
include those sequences above a certain TPM cutoff value.

Currently only Kallisto abundance input is supported (TPM), but that's easily expanded if
you create an issue here with example files: https://github.com/jorvis/biocode/issues

The FASTA file is not read into memory, allowing this to be run on large input files
without issue.

Sequences with score higher than or equal to the --tpm_cutoff will be exported to the
output file.

"""

import argparse
import os
import re


def main():
    parser = argparse.ArgumentParser( description='Filter FASTA file using abundance estimation output')

    parser.add_argument('-if', '--input_fasta', type=str, required=True, help='Input FASTA file' )
    parser.add_argument('-ia', '--input_abundance', type=str, required=True, help='Input abundance estimation file' )
    parser.add_argument('-tc', '--tpm_cutoff', type=float, required=True, help='TPM cutoff to use' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
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

    current_header_line = None
    keep_seq = False

    ofh = open(args.output_file, 'wt')

    with open(args.input_fasta) as f:
        for line in f:
            line = line.rstrip()
            
            if line.startswith('>'):
                m = re.match(">(\S+)", line)

                if m:
                    seq_id = m.group(1)
                else:
                    raise Exception("Failed to parse sequence ID from FASTA header line: {0}".format(line))

                if seq_id not in abund:
                    raise Exception("Failed to find a score for sequence:{0} within the abundance file".format(seq_id))

                keep_seq = True if abund[seq_id] >= args.tpm_cutoff else False
                
            if keep_seq:
                ofh.write("{0}\n".format(line))

    ofh.close() 

                    
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







