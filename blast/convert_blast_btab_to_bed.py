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

"""

import argparse
import os

from biocode import gff

def main():
    parser = argparse.ArgumentParser( description='Converts BLAST and RAPSearch2 tabular output to BED')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-m', '--map_to_gff_coords', type=str, required=False, help='GFF3 file to use to map coordinates onto genomic space' )
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        raise Exception("Input file passed {0} wasn't found".format(args.input_file))

    ofh = open(args.output_file, 'wt')

    if args.map_to_gff_coords is not None:
        (assemblies, features) = gff.get_gff3_features(args.map_to_gff_coords)

    for line in open(args.input_file):
        if line.startswith('#'): continue
        line = line.rstrip()
        cols = line.split("\t")
        if len(cols) < 12: continue

        strand = '+' if int(cols[6]) < int(cols[7]) else '-'

        if args.map_to_gff_coords is None:
            start = int(cols[6]) - 1
            stop = cols[7]
            mol_id = cols[0]
        else:
            if cols[0] in features:
                feat_loc = features[cols[0]].location()
                start = feat_loc.fmin + 1
                stop = feat_loc.fmax
                mol_id = feat_loc.on.id
            else:
                raise Exception("ERROR: Failed to find feature {0} in GFF3 feature set".format(cols[0]))

        ofh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(mol_id, start, stop, cols[1], int(float(cols[11])), strand))

    ofh.close()


if __name__ == '__main__':
    main()







