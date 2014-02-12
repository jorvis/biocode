#!/usr/bin/env python3

import argparse
import os
import biothings
import biocodegff


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-o1a', '--organism1_annotation', type=str, required=True, help='Annotation GFF for organism 1' )
    parser.add_argument('-o1p', '--organism1_protein_alignments', type=str, required=True, help='Path to AAT GFF3 (match/match_part)' )
    #parser.add_argument('-o2a', '--organism2_annotation', type=str, required=True, help='Annotation GFF for organism 2' )
    args = parser.parse_args()

    print("INFO: Parsing organism1 annotation")
    (assemblies, features) = biocodegff.get_gff3_features( args.organism1_annotation )
    # keys are assembly IDs, value for each is a list of matches on them
    aat_matches = dict()
    aat_match_count = 0

    ## IDs of features in organism 1 which overlap AAT
    o1_with_aat = list()

    print("INFO: Parsing organism1 protein alignments")
    for line in open(args.organism1_protein_alignments):
        cols = line.split("\t")

        if len(cols) != 9 or cols[2] != 'nucleotide_to_protein_match':
            continue

        assembly_id = cols[0]

        # skip this match if there were not predicted genes on the same assembly
        if assembly_id not in assemblies:
            continue
        
        fmin = int(cols[3]) - 1
        fmax = int(cols[4])
        strand = cols[6]
        match_id = biocodegff.column_9_value(cols[8], 'ID').replace('"', '')

        match = biothings.Match( id=match_id, subclass='nucleotide_to_protein_match', length=fmax - fmin )
        match.locate_on( target=assemblies[assembly_id], fmin=fmin, fmax=fmax, strand=strand )

        if assembly_id not in aat_matches:
            aat_matches[assembly_id] = list()

        aat_matches[assembly_id].append(match)
        aat_match_count += 1

    print("INFO: Parsed {0} protein alignment chains".format(aat_match_count))

    print("INFO: Comparing organism1's mRNAs with AAT match coordinates")
    for assembly_id in assemblies:
        if assembly_id not in aat_matches:
            continue
        
        assembly = assemblies[assembly_id]

        for gene in assembly.genes():
            for mRNA in gene.mRNAs():
                for aat_match in aat_matches[assembly_id]:
                    if mRNA.overlaps_with(aat_match):
                        o1_with_aat.append(mRNA.id)
                        break   # only need to see if one matched

    print("INFO: Found {0} mRNAs in organism 1 with overlapping AAT coordinates".format(len(o1_with_aat)))
            


if __name__ == '__main__':
    main()







