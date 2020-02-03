#!/usr/bin/env python3

"""

The biocode script calculate_coverage_of_reference.py uses a base-level approach at
reporting on the coverage of genes in a reference.  This doesn't work for some
use cases, as reported by J. Da Silva:

Specifically, if a gene contains a highly variable region, this will cause the alignment 
to break. So, even though the gene is there completely, because there is no alignment 
containing a fraction of the gene, that fraction is assumed to be missing in the new 
assembly.

Here is the fix I propose. If the matches to the start and end of a gene in the reference 
genome are both in a single contig in the new assembly, that gene should be considered 
to be present in its entirety, even if some parts are not in an alignment.

---------------------------------------------

This script uses this simpler approach.  Ouput classes for each gene:

unmapped - Neither end of the gene was found within any of the newly-mapped contigs.
split_mapping - Both ends were mapped but on separate contigs in the new assembly
fully_mapped - Both ends were mapped on the same contig in the new assembly
3_prime_mapped_only - You've got this.
5_prime_mapped_only - This one too.

Example input:

[S1]    [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [TAGS]
1       255     75      328     255     254     92.40   2540030 11810   tp.assembly.567468735.1 NODE_38_length_11810_cov_73.2673
3       481     1583618 1583156 479     463     91.74   2540030 1583618 tp.assembly.567468735.1 NODE_1_length_1583618_cov_74.0429
1732    2389    1759    2413    658     655     90.90   2540030 3677    tp.assembly.567468735.1 NODE_60_length_3677_cov_38.595
2210    9290    1253    8326    7081    7074    94.18   2540030 11810   tp.assembly.567468735.1 NODE_38_length_11810_cov_73.2673
10092   10694   8552    9154    603     603     92.74   2540030 11810   tp.assembly.567468735.1 NODE_38_length_11810_cov_73.2673

Validation:

Gene:FFF274896AFEDB774BE80BF4026117B0 on tp.assembly.567497685.1 at 548145-549821
	matching NODE_12_length_224693_cov_77.5294 at 547008-548379
	matching NODE_2_length_910331_cov_77.0225 at 548104-548472

tp.assembly.567497685.1	.	gene	548146	549821	.	+	.	ID=FFF274896AFEDB774BE80BF4026117B0
547008	548379	22890	21508	1372	1383	90.46	570487	224693	tp.assembly.567497685.1	NODE_12_length_224693_cov_77.5294
548104	548472	15028	14645	369	384	82.73	570487	910331	tp.assembly.567497685.1	NODE_2_length_910331_cov_77.0225
"""

import argparse
import os
from biocode import utils, gff, things
import sys

def main():
    parser = argparse.ArgumentParser( description='Gene coverage calculator')

    parser.add_argument('-c', '--coords', type=str, required=True, help='Coords file, probably from a nucmer alignment' )
    parser.add_argument('-a', '--annotation_file', type=str, required=True, help='GFF3 annotation of input' )
    parser.add_argument('-r', '--reference_fasta', type=str, required=True, help='Reference FASTA from GFF3 coordinates' )
    parser.add_argument('-p', '--protein_encoding_only', action='store_true', help='Pass if you only want protein-coding genes included' )
    parser.add_argument('-d', '--debugging', action='store_true', help='If passed, prints debugging output on STDERR' )
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.annotation_file)

    matches = load_matches(args.coords, assemblies)

    for assembly_id in assemblies:
        assembly = assemblies[assembly_id]

        # some contigs have no mappings
        if assembly_id not in matches: continue
        match = matches[assembly_id]
        
        for gene in assembly.genes():
            if args.protein_encoding_only:
                if len(gene.mRNAs()) == 0:
                    if args.debugging:
                        print("DEBUG: Skipping non-encoding gene {0}".format(gene.id), file=sys.stderr)
                    continue
            
            fmin_matches = list()
            fmax_matches = list()
            gene_loc = gene.location_on(assembly)

            if args.debugging:
                print("DEBUG: Gene:{0} on {3} at {1}-{2}".format(gene.id, gene_loc.fmin, gene_loc.fmax, gene_loc.on.id), file=sys.stderr)

            for mp in match.parts:
                match_loc = mp.location_on(assembly)
                
                if mp.overlaps_min_side_of(thing=gene):
                    if args.debugging:
                        print("\tmatching {0} at {1}-{2}".format(mp.id, match_loc.fmin, match_loc.fmax), file=sys.stderr)
                        
                    fmin_matches.append(mp.id)

                if mp.overlaps_max_side_of(thing=gene):
                    if args.debugging:
                        print("\tmatching {0} at {1}-{2}".format(mp.id, match_loc.fmin, match_loc.fmax), file=sys.stderr)
                        
                    fmax_matches.append(mp.id)

            # everything starts as unmapped
            match_class = 'unmapped'

            # if there were matches to both ends the class gets upgraded to at least this
            if len(fmin_matches) and len(fmax_matches):
                match_class = 'split_mapping'

                # see if there was a shared match between the two ends
                if shared_match(fmin_matches, fmax_matches):
                    match_class = 'fully_mapped'

            if match_class == 'unmapped':
                # let's see if at least one of the ends is mapped
                if len(fmin_matches):
                    if gene_loc.strand == 1:
                        match_class = '5_prime_mapped_only'
                    else:
                        match_class = '3_prime_mapped_only'

                if len(fmax_matches):
                    if gene_loc.strand == 1:
                        match_class = '3_prime_mapped_only'
                    else:
                        match_class = '5_prime_mapped_only'

            print("{0}\t{1}\t{2}".format(assembly.id, gene.locus_tag, match_class))
                        
            
def shared_match(a, b): 
    a_set = set(a) 
    b_set = set(b) 
    if len(a_set.intersection(b_set)) > 0: 
        return(True)  
    return(False)                    
            
                
def load_matches(infile, assemblies):
    matches = dict()

    for line in open(infile):
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) < 9: continue
        # this matches header lines
        if line[0] == '[': continue

        ref_asmbl_id = cols[9]
        new_asmbl_id = cols[10]
        ref_start = int(cols[0])
        ref_end   = int(cols[1])
        
        if ref_start > ref_end:
            raise Exception("Unexpected start {0} > {1} end".format(ref_start, ref_end))

        if ref_asmbl_id not in matches:
            #matches[ref_asmbl_id] = list()
            matches[ref_asmbl_id] = things.Match()

        #matches[ref_asmbl_id].append({'start': ref_start, 'end': ref_end, 'mol_id': new_asmbl_id})
        mp = things.MatchPart(id=new_asmbl_id, parent=matches[ref_asmbl_id])
        mp.locate_on(target=assemblies[ref_asmbl_id], fmin=ref_start, fmax=ref_end)
        matches[ref_asmbl_id].add_part(mp)

    return matches
        
if __name__ == '__main__':
    main()







