#!/usr/bin/env python3

import argparse
import os
import re
import biothings
import biocodegff
import biocodeutils

"""
Example with a match extending far past a gene
DEBUG: g8713.t1:(3529) overlaps (size:2893) nucleotide_to_protein_match.158742:(6245), match target id:jgi|Copci1|10482|CC1G_12482T0, length:1764
	mRNA % cov: 81.97789742136582
	target % cov: 164.00226757369614

NODE_19891_length_98951_cov_9.347667	AUGUSTUS	mRNA	38709	42237	.	+	.	ID=g8713.t1;Parent=g8713

NODE_19891_length_98951_cov_9.347667    nap     nucleotide_to_protein_match     35357   41601   .       +       .       ID=nucleotide_to_protein_match.158742;Target=jgi|Copci1|10482|CC1G_12482T0 317 504
NODE_19891_length_98951_cov_9.347667    nap     match_part      35357   35498   .       +       .       ID=match_part.493339;Parent=nucleotide_to_protein_match.158742;Target=jgi|Copci1|10482|CC1G_12482T0 317 364a
NODE_19891_length_98951_cov_9.347667    nap     match_part      36352   36741   .       +       .       ID=match_part.493340;Parent=nucleotide_to_protein_match.158742;Target=jgi|Copci1|10482|CC1G_12482T0 365 490
NODE_19891_length_98951_cov_9.347667    nap     match_part      41560   41601   .       +       .       ID=match_part.493341;Parent=nucleotide_to_protein_match.158742;Target=jgi|Copci1|10482|CC1G_12482T0 491 504
"""

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    parser.add_argument('-a', '--organism1_annotation', type=str, required=True, help='Annotation GFF for organism 1' )
    parser.add_argument('-p', '--organism1_aat_alignments', type=str, required=True, help='Path to AAT GFF3 (match/match_part)' )
    parser.add_argument('-aatdb', '--aat_fasta_db', type=str, required=True, help='Path to FASTA database that was used in AAT' )
    parser.add_argument('-b', '--organism1_blast_alignments', type=str, required=True, help='Path to BLASTp btab file vs.organism 2 proteins' )
    parser.add_argument('-be', '--blast_eval_cutoff', type=float, required=False, default=1e-5, help='BLAST e-value cutoff' )
    parser.add_argument('-bpi', '--blast_percent_identity_cutoff', type=float, required=False, default=0, help='BLAST %identity cutoff' )
    parser.add_argument('-ppc', '--aat_percent_coverage_cutoff', type=float, required=False, default=0, help='% coverage of the query protein by the AAT match' )
    parser.add_argument('-o', '--output_id_list', type=str, required=False, help='List of IDs from organism1 that passed' )
    args = parser.parse_args()

    ## if the output file wasn't passed build one from the other parameters
    if args.output_id_list is None:
        args.output_id_list = "training_ids.be_{0}.bpi_{1}.ppc_{2}.list".format(args.blast_eval_cutoff, args.blast_percent_identity_cutoff, args.aat_percent_coverage_cutoff)

    print("INFO: Parsing organism1 annotation")
    (assemblies, features) = biocodegff.get_gff3_features( args.organism1_annotation )

    print("INFO: Parsing AAT FASTA database")
    aat_seqs = biocodeutils.fasta_dict_from_file( args.aat_fasta_db )
    
    # keys are assembly IDs, value for each is a list of matches on them
    aat_matches = dict()
    aat_match_count = 0
    current_match = None

    ## IDs of features in organism 1 which overlap AAT
    o1_with_aat = list()
    o1_with_o2 = list()

    print("INFO: Parsing organism1 protein alignments")
    for line in open(args.organism1_aat_alignments):
        cols = line.split("\t")

        if line.startswith('#') or len(cols) != 9:
            continue

        assembly_id = cols[0]

        # skip this match if there were not predicted genes on the same assembly
        if assembly_id not in assemblies:
            continue

        if assembly_id not in aat_matches:
            aat_matches[assembly_id] = list()
        
        fmin = int(cols[3]) - 1
        fmax = int(cols[4])
        strand = cols[6]
        feature_id = biocodegff.column_9_value(cols[8], 'ID').replace('"', '')
        target = biocodegff.column_9_value(cols[8], 'Target')
        m = re.search("^(\S+)", target)
        if m:
            target = m.group(1)

        if cols[2] == 'nucleotide_to_protein_match':
            if current_match is not None:
                aat_matches[assembly_id].append(current_match)
                aat_match_count += 1
            
            current_match = biothings.Match( id=feature_id, target_id=target, subclass='nucleotide_to_protein_match', length=fmax - fmin )
            current_match.locate_on( target=assemblies[assembly_id], fmin=fmin, fmax=fmax, strand=strand )

        elif cols[2] == 'match_part':
            parent_id = biocodegff.column_9_value(cols[8], 'Parent').replace('"', '')
            match_part = biothings.MatchPart( id=feature_id, parent=parent_id, length=fmax - fmin )
            current_match.add_part(match_part)


    print("INFO: Parsed {0} protein alignment chains".format(aat_match_count))

    print("INFO: Comparing organism1's mRNAs with AAT match coordinates")
    for assembly_id in assemblies:
        if assembly_id not in aat_matches:
            continue
        
        assembly = assemblies[assembly_id]

        for gene in assembly.genes():
            for mRNA in gene.mRNAs():

                for aat_match in aat_matches[assembly_id]:
                    overlap_size = mRNA.overlap_size_with(aat_match)

                    if overlap_size is not None:
                        #print("DEBUG: {0}:({1}) overlaps (size:{2}) {3}:({4})".format(mRNA.id, mRNA.length, overlap_size, aat_match.id, aat_match.length) )
                        # this shouldn't be possible, but check just in case
                        if overlap_size > mRNA.length:
                            raise Exception("ERROR: overlap size ({0}) > mRNA length ({1})".format(overlap_size, mRNA.length))

                        if aat_match.target_id not in aat_seqs:
                            raise Exception("ERROR: Found match with target ID ({0}) but didn't find a FASTA entry for it via -aatdb".format(aat_match.target_id))

                        # 
                        # this is a protein length, so x3
                        match_target_length = len(aat_seqs[aat_match.target_id]['s']) * 3
                        
                        mRNA_percent_coverage = (overlap_size / mRNA.length) * 100
                        target_percent_coverage = (overlap_size / match_target_length) * 100

                        #if mRNA_percent_coverage >= args.aat_percent_coverage_cutoff and target_percent_coverage >= args.aat_percent_coverage_cutoff:
                        if mRNA_percent_coverage == 100 and target_percent_coverage == 100:
                            o1_with_aat.append(mRNA.id)
                            #print("DEBUG: {0}:({1}) overlaps (size:{2}) {3}:({4}), match target id:{5}, length:{6}".format(mRNA.id, mRNA.length, overlap_size, aat_match.id, aat_match.length, aat_match.target_id, match_target_length) )
                            #print("\tmRNA % cov: {0}".format(mRNA_percent_coverage))
                            #print("\ttarget % cov: {0}".format(target_percent_coverage))
                            break   # only need to see if one matched

    print("INFO: Found {0} mRNAs in org1 with overlapping fungi AAT coordinates".format(len(o1_with_aat)))

    # key=org1_transcript_id, value=org2_transcript_id
    top_blast_hits = dict()

    print("INFO: parsing BLAST results vs. org2")
    for line in open(args.organism1_blast_alignments):
        cols = line.split("\t")

        if float(cols[19]) > args.blast_eval_cutoff:
            continue

        if float(cols[10]) < args.blast_percent_identity_cutoff:
            continue
        
        # if we survived until here, this one's good.
        top_blast_hits[cols[0]] = cols[5]

    print("INFO: Comparing overlap between AAT-matched proteins and BLAST ones")
    for o1_mRNA_id in o1_with_aat:
        if o1_mRNA_id in top_blast_hits:
            o1_with_o2.append(o1_mRNA_id)

    print("INFO: Found {0} mRNAs in org1 with overlapping AAT coordinates and BLAST hit to org2".format(len(o1_with_o2)))

    id_list_fh = open(args.output_id_list, 'wt')
    for mRNA_id in o1_with_o2:
        id_list_fh.write("{0}\n".format(mRNA_id))


if __name__ == '__main__':
    main()







