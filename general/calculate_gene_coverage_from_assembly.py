#!/usr/bin/env python3

"""

Original use case description:

Specifically, if a gene contains a highly variable region, this will cause the alignment 
to break. So, even though the gene is there completely, because there is no alignment 
containing a fraction of the gene, that fraction is assumed to be missing in the new assembly.

Here is the fix I propose. If the matches to the start and end of a gene in the reference 
genome are both in a single contig in the new assembly, that gene should be considered to 
be present in its entirety, even if some parts are not in an alignment.

This is a demo to test this.

Input nucmer coords format:
[S1] start of the alignment region in the reference sequence 
[E1] end of the alignment region in the reference sequence 
[S2] start of the alignment region in the query sequence 
[E2] end of the alignment region in the query sequence 
[LEN 1] length of the alignment region in the reference sequence 
[LEN 2] length of the alignment region in the query sequence 
[% IDY] percent identity of the alignment 
[LEN R] length of the reference sequence 
[LEN Q] length of the query sequence 

Output:

The output is a 2-column tab-delimited file where the first column is the gene ID
and the second is it's coverage class.  Here are the possible values for those:

single_contig_full_coverage - Both sides of the gene are covered by the same 
aligned contig.

multi_contig_full_coverage - Both sides of the are covered, but by different
aligned contigs.

fmax_coverage - Only the right-most side of the gene (unstranded) has coverage.

fmin_coverage - Only the left-most side of the gene (unstranded) has coverage.

no_coverage - Neither side of this gene overlaps an aligned contig.

"""

import argparse
import os

from biocode import gff, things

DEBUG_COORDINATES = False

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')
    parser.add_argument('-c', '--coords_file', type=str, required=True, \
                            help='Path to a nucmer coords file with non-overlapping results (requires -l -r -T options of show-coords)' )
    parser.add_argument('-a', '--annotation_file', type=str, required=True, help='Path to a sorted GFF3 annotation file' )
    parser.add_argument('-r', '--reference_fasta', type=str, required=True, help='Path to the reference file used with nucmer' )
 
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.annotation_file)

    # matches[ref_id][qry_id]
    matches = parse_nucmer_matches(assemblies, args.coords_file)

    for assembly_id in assemblies:
        assembly = assemblies[assembly_id]

        for gene in assembly.genes():
            start_covered = False
            end_covered = False
            full_single_contig_coverage = False

            if DEBUG_COORDINATES:
                print("gene:{0}\tat coordinates {1}-{2}".format(gene.id, gene.location().fmin, gene.location().fmax))
            else:
                print(gene.id, end='\t')

            # check that the assembly has any matches at all
            if assembly_id in matches:
                for contig_id in matches[assembly_id]:
                    if DEBUG_COORDINATES: print("\tchecking matches from contig: {0}".format(contig_id))
                    start_covered_by_this_contig = False
                    end_covered_by_this_contig = False

                    for match in matches[assembly_id][contig_id].parts:
                        if DEBUG_COORDINATES: print("\t\t{0}-{1}".format(match.location().fmin, match.location().fmax))
                        if gene.overlaps_with(match):
                            if DEBUG_COORDINATES: print("\t\t\toverlaps")
                            # is it completely contained
                            if gene.contained_within(match):
                                if DEBUG_COORDINATES: print("\t\t\tcontained_within")
                                start_covered_by_this_contig = True
                                end_covered_by_this_contig = True

                            # ok, what about the start side?
                            elif match.overlaps_min_side_of(gene):
                                start_covered = True
                                start_covered_by_this_contig = True

                            elif match.overlaps_max_side_of(gene):
                                end_covered = True
                                end_covered_by_this_contig = True

                    # no more comparisons need to happen. report the type of this one
                    if start_covered_by_this_contig and end_covered_by_this_contig:
                        full_single_contig_coverage = True
                        break

            if full_single_contig_coverage:
                print("single_contig_full_coverage")
            elif start_covered and end_covered:
                print("multi_contig_full_coverage")
            elif start_covered:
                print("fmin_coverage")
            elif end_covered:
                print("fmax_coverage")
            else:
                print("no_coverage")
    
def parse_nucmer_matches(assemblies, match_file):
    # [ref_id][qry_id] = Match
    matches = dict()
    
    with open(match_file) as fh:
        in_alignment_region = False
        next_match_num = 1
        
        for line in fh:
            if line.startswith('[S1]'):
                in_alignment_region = True
                continue

            if in_alignment_region:
                line = line.rstrip()
                cols = line.split()
                ref_id = cols[9]
                qry_id = cols[10]

                if int(cols[0]) < int(cols[1]):
                    ref_start = int(cols[0]) - 1
                    ref_end = int(cols[1])
                    ref_strand = 1
                else:
                    ref_start = int(cols[1]) - 1
                    ref_end = int(cols[0])
                    ref_strand = -1

                if ref_id not in matches:
                    matches[ref_id] = dict()

                if qry_id not in matches[ref_id]:
                    match = things.Match(id="match_{0}".format(next_match_num), subclass='match')
                    next_match_num += 1
                    matches[ref_id][qry_id] = match

                # check assumptions
                if ref_start > ref_end:
                    raise Exception("Error, unexpected ref_start ({0}) > ref_end ({1})".format(ref_start, ref_end))

                # each line is a MatchPart
                match_part = things.MatchPart(parent=matches[ref_id][qry_id], length=cols[8])
                match_part.locate_on(target=assemblies[ref_id], fmin=ref_start, fmax=ref_end, strand=ref_strand)
                matches[ref_id][qry_id].add_part(match_part)

    return matches

if __name__ == '__main__':
    main()







