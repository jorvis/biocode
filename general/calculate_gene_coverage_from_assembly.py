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

------------------------------------------------------------------------

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

------------------------------------------------------------------------

Validation example:

tp.assembly.567485385.1 TpMuguga_04g00003       split_mapping   65

DEBUG: Gene:2A9AF11CDD8F0F735E26394CF71B9EF6 on tp.assembly.567485385.1 at 8499-10209
        matching NODE_33_length_23380_cov_77.6652 at 6391-9449
        matching NODE_68_length_2687_cov_71.4871 at 10034-12716

[S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[LEN R]	[LEN Q]	[TAGS]
6391	9449	11945	15004	3059	3060	97.45	17691	23380	tp.assembly.567485385.1	NODE_33_length_23380_cov_77.6652
10034	12716	36	2687	2683	2652	94.13	17691	2687	tp.assembly.567485385.1	NODE_68_length_2687_cov_71.4871


                                                        10034              12716
                                                        |------------------|  contig:NODE_68_length_2687_cov_71.4871

      6391                                  9449
      |-------------------------------------|  contig:NODE_33_length_23380_cov_77.6652

                                 8499                     10209
5000                             |------------------------| gene: TpMuguga_04g00003                   15000
|-----------------------------------------------------------------------------------------------------| tp.assembly.567485385.1

Coverage blocks:
10209 - 10034 = 175
9449 - 8499   = 950
              -----
               1125bp of 1710 gene length covered
               65.7894737%

------------------------------------------------------------------------

Overlapping contig validation example (coverage maxed at 100%)

DEBUG: Gene:3CCE302F0B1982C331747490BC7E3B1A on tp.assembly.567478335.1 at 211885-212339
	matching NODE_117_length_594_cov_61.1456 at 211425-212019
	matching NODE_89_length_1565_cov_80.7712 at 211849-212321
	matching NODE_3_length_840898_cov_75.6474 at 212227-306140

           Error: coverage exceeded 100%, which shouldn't happen.
           Gene length: 454
           Fmin max length: 436
           Fmax max length: 112
           Got best coverage value > 100 (120)

                                     212227                                                                                     306140   
                                     |-------------------------------------------------------------------------------------------...  contig:NODE_3_length_840898_cov_75.6474

           211849                    212321
           |-------------------------|  contig:NODE_89_length_1565_cov_80.7712

         211425                     212019
         |--------------------------|  contig:NODE_117_length_594_cov_61.1456
         
            211885                    212339
210000      |-------------------------|  gene:TpMuguga_02g02195                                                       214000
|---------------------------------------------------------------------------------------------------------------------| tp.assembly.567478335.1

212339 - 212227 = 112
212321 - 211885 = 436

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
                        
                    fmin_matches.append(mp)

                if mp.overlaps_max_side_of(thing=gene):
                    if args.debugging:
                        print("\tmatching {0} at {1}-{2}".format(mp.id, match_loc.fmin, match_loc.fmax), file=sys.stderr)
                        
                    fmax_matches.append(mp)

            # everything starts as unmapped
            match_class = 'unmapped'
            coverage = 0

            # if there were matches to both ends the class gets upgraded to at least this
            if len(fmin_matches) and len(fmax_matches):
                match_class = 'split_mapping'

                # see if there was a shared match between the two ends
                if shared_match(fmin_matches, fmax_matches):
                    match_class = 'fully_mapped'
                    coverage = 100
                else:
                    coverage = calculate_coverage(gene_loc, fmin_matches, fmax_matches)

            if match_class == 'unmapped':
                # let's see if at least one of the ends is mapped
                if len(fmin_matches):
                    coverage = calculate_coverage(gene_loc, fmin_matches, fmax_matches)
                    if gene_loc.strand == 1:
                        match_class = '5_prime_mapped_only'
                    else:
                        match_class = '3_prime_mapped_only'

                if len(fmax_matches):
                    coverage = calculate_coverage(gene_loc, fmin_matches, fmax_matches)
                    if gene_loc.strand == 1:
                        match_class = '3_prime_mapped_only'
                    else:
                        match_class = '5_prime_mapped_only'

            print("{0}\t{1}\t{2}\t{3}".format(assembly.id, gene.locus_tag, match_class, coverage))
                        

def calculate_coverage(geneloc, fmins, fmaxes):
    """
    This is to report the percent coverage of a gene on its ends due to
    partial matches to two different assembly contigs. Since each end 
    could map to more than on contig, only the longest match on each end
    is considered.
    """
    gene_length = geneloc.fmax - geneloc.fmin
    best_coverage = 0
    fmin_max_length = 0
    fmax_max_length = 0
    
    for fmin_match in fmins:
        if len(fmin_match.locations) > 1:
            raise Exception("Unexpected.  Got match with more than one location")

        fmin_loc = fmin_match.locations[0]
        fmin_match_length = fmin_loc.fmax - geneloc.fmin

        if fmin_match_length > fmin_max_length:
            fmin_max_length = fmin_match_length

    for fmax_match in fmaxes:
        if len(fmax_match.locations) > 1:
            raise Exception("Unexpected.  Got match with more than one location")

        fmax_loc = fmax_match.locations[0]
        fmax_match_length = geneloc.fmax - fmax_loc.fmin

        if fmax_match_length > fmax_max_length:
            fmax_max_length = fmax_match_length

    best_coverage = int(((fmin_max_length + fmax_max_length) / gene_length) * 100)

    if best_coverage > 100:
        msg = """
           Error: coverage exceeded 100%, indicating missed contig assembly.  Reducing to 100%.
           Gene length: {0}
           Fmin max length: {2}
           Fmax max length: {3}
           Got best coverage value > 100 ({1})
        """.format(gene_length, best_coverage, fmin_max_length, fmax_max_length)
        best_coverage = 100
        print(msg, file=sys.stderr)
    
    return best_coverage
                
            
def shared_match(a, b): 
    for qry_match in a:
        for sbj_match in b:
            if qry_match.id == sbj_match.id:
                return True
    
    return False
            
                
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







