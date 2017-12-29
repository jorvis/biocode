#!/usr/bin/env python3

"""

We now often deal with BAM files which are 1TB in size or more.  When doing evaluations
of different utilities, we want mpileup output for a portion of the entries within
the huge BAM file, and generating mpileup for the entire thing can take a week.

This script accepts a list of IDs and uses the index from a BAM file to extract only
the mpileup output for that set of IDs.  If you've aligned billions of reads to a set
of transcripts, for example, this would accept a list of transcript IDs for which
you want coverage output.

NOTE:  By default, some read filtering options are turned off via the samtools mpileup
options "-A" and "-Q 0".  You can add more or set these to an empty string with
the --samtools_options parameter.

"""

import argparse
import os
import datetime
import subprocess


def main():
    parser = argparse.ArgumentParser( description='Generates an mpileup file for a subset of molecules in a large BAM file')

    ## output file to be written
    parser.add_argument('-b', '--bam_file', type=str, required=True, help='Path to an input BAM file to be read (actually, the index base name/path)' )
    parser.add_argument('-l', '--list_file', type=str, required=True, help='Path to an input list file to be read' )
    parser.add_argument('-f', '--fasta_file', type=str, required=True, help='FASTA file corresponding to the IDs in the list file' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-sp', '--samtools_path', type=str, required=False, default='samtools', help='Use a specific path of samtools, else PATH will be used' )
    parser.add_argument('-so', '--samtools_options', type=str, required=False, default='-A -Q 0', help='Options to pass to samtools mpileup' )
    args = parser.parse_args()

    # make sure output file starts out empty
    fh = open(args.output_file, 'wt')
    fh.close()

    for mol_id in open(args.list_file):
        mol_id = mol_id.rstrip()
        if len(mol_id) < 1: continue

        cmd = "{0} mpileup {4} -f {1} {2} -r '{3}' >> {5}".format(args.samtools_path, args.fasta_file, args.bam_file, mol_id, args.samtools_options, args.output_file)
        run_command(cmd)
        
        
def run_command(cmd):
    print("INFO: Running command: {0}".format(cmd), flush=True)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        raise Exception("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))


if __name__ == '__main__':
    main()







