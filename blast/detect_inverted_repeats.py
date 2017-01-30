#!/usr/bin/env python3

"""

The product of some transcriptome assemblies had sequences which contained large 
inverted repeats within some transcripts.  We wanted a quick utility to detect
these and generate a simple report on them.

The input is just a FASTA file.  This script will iterate over them, writing each to
a temporary file, then calling blastn (which must be in your PATH) and parsing the
results.  This has only been tested on transcriptome files, not genomic.

The current output will look like this:

INVERSION of 1225 bp in TR61583|c7_g2_i3: 1     1225    2765    1541
INVERSION of 112 bp in TR100233|c9_g5_i2: 1     112     112     1
INVERSION of 1391 bp in TR109248|c7_g2_i5: 3143 4533    1435    45
INVERSION of 1391 bp in TR109248|c7_g2_i5: 45   1435    4533    3143
INVERSION of 1343 bp in TR109248|c7_g2_i5: 1618 2960    2960    1618
INVERSION of 1071 bp in TR64784|c1_g3_i8: 1159  2229    1071    1
INVERSION of 1071 bp in TR64784|c1_g3_i8: 1     1071    2229    1159

Messages will be printed on STDOUT after every 100 entries in the input file
are processed.

"""

import argparse
import datetime
import os
import subprocess

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Use BLAST to identify internal inverted repeats')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-n', '--min_repeat_size', type=int, required=True, help='Minimum size of a repeat to consider' )
    parser.add_argument('-pid', '--percent_identity', type=float, required=False, default=98.0, help='Percent identity cutoff' )
    args = parser.parse_args()

    # parse FASTA input, storing into a dict keyed by ID
    seqs = utils.fasta_dict_from_file(args.input_file)

    ofh = open(args.output_file, 'wt')

    seqs_processed = 0
    print_interval = 100

    for id in seqs:
        # Write a FASTA of just this sequence
        fasta_name = "{0}.temp.input.fasta".format(os.getpid())
        blast_name = "{0}.temp.blast.out".format(os.getpid())
        fasta_fh = open(fasta_name, 'wt')
        fasta_fh.write(">{0}\n{1}".format(id, seqs[id]['s']))
        fasta_fh.close()

        # Perform the blast using bl2seq
        #cmd = "bl2seq -i {0} -j {0} -p blastn -e 1e-10 -D 1 -o {1} -W {2}".format(fasta_name, blast_name, args.min_repeat_size)
        cmd = "blastn -query {0} -subject {0} -outfmt 6 -out {1} -word_size {2} -perc_identity {3}".format(fasta_name, blast_name, args.min_repeat_size, args.percent_identity)
        run_command(cmd)

        # Parse the result file to look for inverted repeats
        for line in open(blast_name):
            if line.startswith('#'):
                continue

            cols = line.split()
            qstart, qend, sstart, send = int(cols[6]), int(cols[7]), int(cols[8]), int(cols[9])

            if qstart < qend:
                q_orientation = 'F'
                match_len = qend - qstart + 1
            else:
                q_orientation = 'R'
                match_len = qstart - qend + 1

            if sstart < send:
                s_orientation = 'F'
            else:
                s_orientation = 'R'

            if s_orientation != q_orientation and match_len >= args.min_repeat_size:
                ofh.write("INVERSION of {5} bp in {4}: {0}\t{1}\t{2}\t{3}\n".format(qstart, qend, sstart, send, cols[0], match_len))

            if s_orientation == q_orientation and match_len >= args.min_repeat_size:
                if (qstart >= sstart and qstart <= send) or (qend >= sstart and qend <= send):
                    pass
                else:
                    ofh.write("DIRECT REPEAT of {5} bp in {4}: {0}\t{1}\t{2}\t{3}\n".format(qstart, qend, sstart, send, cols[0], match_len))
                    #ofh.write("# ^^ {0}".format(line))

        seqs_processed += 1

        if seqs_processed % print_interval == 0:
            print("INFO: processed {0} input sequences".format(seqs_processed))

    ofh.close()
    
def run_command(cmd):
    #print("INFO: Running command: {0}",format(cmd), flush=True)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        print("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))
        #raise Exception("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))


if __name__ == '__main__':
    main()







