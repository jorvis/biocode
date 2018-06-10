#!/usr/bin/env python3

"""

The input here is the HTAB output from HMMer searches.  If the --map_to_gff_coords option
is passed, the query IDs will be searched within the passed GFF3 file so that the 
coordinates of each match can be transformed onto genomic space (appropriate for 
changing protein HMM results into files usable by genome browsers.)

Input:

    col  perl-col   description
    1      [0]      HMM accession
    2      [1]      Date search was run (if available), otherwise date of htab parse
    3      [2]      Length of the HMM (not populated if -s is used)
    4      [3]      Search program
    5      [4]      Database file path
    6      [5]      Sequence accession
    7      [6]      Alignment start position on HMM match - hmm-f
    8      [7]      Alignment end position on HMM match - hmm-t
    9      [8]      Alignment start position on sequence - seq-f
    10     [9]      Alignment end position on sequence - seq-t
    11     [10]     frame (only populated if --frames search is run on nucleotide sequence)
    12     [11]     Domain score
    13     [12]     Total score
    14     [13]     Index of domain hit
    15     [14]
    16     [15]     HMM description (may be truncated by hmmsearch or hmmpfam if -s is used)
    17     [16]     Sequence description (may be truncated by hmmsearch or hmmpfam)
    18     [17]     Total score trusted cutoff (not populated if -s is used)
    19     [18]     Total score noise cutoff (not populated if -s is used)
    20     [19]     Expect value for total hit
    21     [20]     Expect value for domain hit
    22     [21]     Domain score trusted cutoff (egad..hmm2.trusted_cutoff2) (not populated if -s is used)
    23     [22]     Domain score noise cutoff (egad..hmm2.noise_cutoff2) (not populated if -s is used)
    24     [23]     Total score gathering threshold (not populated if -s is used)
    25     [24]     Domain score gathering threshold (not populated if -s is used)


The output is a BED file with the following tabs:

<queryID> <start> <end> <matchID> <score> <strand>

What goes in the 'score' field depends on the options you pass.  In BED format the value must be
an integer between 0 and 1000, which doesn't let itself well to being either Total Score or E-value.

By default the Total Score is used but maxes at 1000.  If you pass --score_by=evalue we have the
same issue of needing to scale it into something useful between 0 and 1000.  Here I've chosen to
use 0 - log(e-value) limiting the range within that span.

"""

import argparse
import math
import os
import re

from biocode import gff

def main():
    parser = argparse.ArgumentParser( description='Converts HMMer tabular output to BED')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-m', '--map_to_gff_coords', type=str, required=False, help='GFF3 file to use to map coordinates onto genomic space' )
    parser.add_argument('-s', '--score_by', type=str, required=False, default='totalscore', help='Populates the score column by totalscore or evalue' )
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        raise Exception("Input file passed {0} wasn't found".format(args.input_file))

    ofh = open(args.output_file, 'wt')

    if args.map_to_gff_coords is not None:
        (assemblies, features) = gff.get_gff3_features(args.map_to_gff_coords)

    line_num = 0

    for line in open(args.input_file):
        line_num += 1

        if line.startswith('#'): 
            continue

        line = line.rstrip()
        cols = line.split("\t")
        if len(cols) < 24: continue

        if args.map_to_gff_coords is None:
            start = cols[8]
            stop = cols[9]
            mol_id = cols[5]
            strand = '+'
        else:
            if cols[5] in features:
                feat_loc = features[cols[5]].location()
                strand = '+' if feat_loc.strand == 1 else '-'

                if strand == '+':
                    start = feat_loc.fmin + ((int(cols[8]) - 1) * 3)
                    stop  = feat_loc.fmin + (int(cols[9]) * 3)
                else:
                    start = feat_loc.fmin + ((int(cols[8]) - 1) * 3)
                    stop  = feat_loc.fmin + ((int(cols[9]) + 1) * 3)

                mol_id = feat_loc.on.id
            else:
                raise Exception("ERROR: Failed to find feature {0} in GFF3 feature set".format(cols[0]))

        if args.score_by == 'totalscore':
            adj_score = int(float(cols[12]))
        else:
            e_val = float(cols[19])
            
            if e_val == 0:
                adj_score = 1000
            else:
                try:
                    adj_score = int(0 - math.log(e_val))
                except ValueError:
                    raise Exception("ERROR doing math on this log value: ({0}) on line {1}".format(e_val, line_num))

        if adj_score < 0:
            adj_score = 0
        elif adj_score > 1000:
            adj_score = 1000
            
        ofh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(mol_id, start, stop, cols[0], adj_score, strand))

    ofh.close()


if __name__ == '__main__':
    main()







