#!/usr/bin/env python3

"""

The input here is supposed to be either BLAST output using the -m 8 or -m 9
options or also the tabular output from RAPSearch2.  If the --map_to_gff_coords option
is passed, the query IDs will be searched within the passed GFF3 file so that the 
coordinates of each match can be transformed onto genomic space (appropriate for 
changing protein BLASTS into files usable by genome browsers.)

Input:

1. Query 
2. Subject 
3. % identity
4. aln-len
5. mismatch
6. gap-openings
7. q.start
8. q.end
9. s.start
10. s.end
11. E-value (or log() for RAPSearch2)
12. bit-score

The output is a BED file with the following tabs:

<queryID> <start> <end> <matchID> <score> <strand>

What goes in the 'score' field depends on the options you pass.  In BED format the value must be
an integer between 0 and 1000, which doesn't let itself well to being either bitscore or E-value.

By default the bitscore is used but maxes at 1000.  If you pass --score_by=evalue we have the
same issue of needing to scale it into something useful between 0 and 1000.  Here I've chosen to
use 0 - log(e-value) limiting the range within that span.

"""

import argparse
import math
import os
import re

from biocode import gff

def main():
    parser = argparse.ArgumentParser( description='Converts BLAST and RAPSearch2 tabular output to BED')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-m', '--map_to_gff_coords', type=str, required=False, help='GFF3 file to use to map coordinates onto genomic space' )
    parser.add_argument('-s', '--score_by', type=str, required=False, default='bitscore', help='Populates the score column by bitscore or evalue' )
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        raise Exception("Input file passed {0} wasn't found".format(args.input_file))

    ofh = open(args.output_file, 'wt')

    if args.map_to_gff_coords is not None:
        (assemblies, features) = gff.get_gff3_features(args.map_to_gff_coords)

    rapsearch_detected = False

    for line in open(args.input_file):
        if line.startswith('#'): 
            m = re.search('RAPSearch', line)
            if m:
                rapsearch_detected = True

            continue

        line = line.rstrip()
        cols = line.split("\t")
        if len(cols) < 12: continue

        if args.map_to_gff_coords is None:
            start = cols[6]
            stop = cols[7]
            mol_id = cols[0]
            strand = '+'
        else:
            if cols[0] in features:
                feat_loc = features[cols[0]].location()
                strand = '+' if feat_loc.strand == 1 else '-'
                
                protein_length = int(cols[7]) - int(cols[6]) + 1

                start = feat_loc.fmin + ((int(cols[6]) - 1) * 3)
                stop = start + ((protein_length + 1) * 3)
                    
                mol_id = feat_loc.on.id
            else:
                raise Exception("ERROR: Failed to find feature {0} in GFF3 feature set".format(cols[0]))

        if args.score_by == 'bitscore':
            adj_score = int(float(cols[11]))
        else:
            e_val = float(cols[10])
            if rapsearch_detected == True:
                adj_score = int(0 - e_val)
            else:
                adj_score = int(0 - math.log(e_val))

        if adj_score < 0:
            adj_score = 0
        elif adj_score > 1000:
            adj_score = 1000
            
        ofh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(mol_id, start, stop, cols[1], adj_score, strand))

    ofh.close()


if __name__ == '__main__':
    main()







