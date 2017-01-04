#!/usr/bin/env python3

"""
When doing gene finding with tools like Augustus it's common to take a file of curated
genes/transcripts and split these into a set for training and another for evaluation.

This script accepts an input GFF3 file and then generates two files, one of genes to
use as training and another as evaluation.  Selection of these is (computationally) random.

If the --retain_composition option is used, the script will attempt to make the exon count
distribution remain consisent between the full file, training set and evaluation set rather
than relying on pure random selection to do this.

Run this with the -h option to view other parameters and defaults.

Author: Joshua Orvis
"""

import argparse
import random
from collections import defaultdict

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Split an annotation GFF3 into training and evaluation sets')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-ot', '--output_training_file', type=str, required=True, help='GFF3 file to be created with the training genes' )
    parser.add_argument('-oe', '--output_evaluation_file', type=str, required=True, help='GFF3 file to be created with the evaluation genes' )
    parser.add_argument('-ts', '--training_set_size', type=int, required=False, default=200, help='Number of transcripts to select for training' )
    parser.add_argument('-es', '--evaluation_set_size', type=int, required=False, default=100, help='Number of transcripts to select for evaluation' )
    parser.add_argument('-me', '--max_exon_count', type=int, required=False, help='Skips any mRNAs with more exons than this' )
    parser.add_argument('--retain_composition', dest='retain_composition',action='store_true')
    parser.add_argument('--no_retain_composition', dest='retain_composition',action='store_false')
    parser.set_defaults(retain_composition=False)
    args = parser.parse_args()

    if args.retain_composition is True:
        raise Exception("ERROR: --retain_composition option not yet implemented")

    (assemblies, features) = gff.get_gff3_features(args.input_file)

    # key: exon count, value = list of mRNA objects with that count
    # which of these gets used depends on whether --retain_composition is passed
    mRNAs_by_exon_count = defaultdict(lambda: list())
    mRNAs = list()
    mRNA_count = 0

    for asm_id in assemblies:
        for gene in assemblies[asm_id].genes():
            for mRNA in gene.mRNAs():
                exon_count = mRNA.exon_count()

                if args.max_exon_count is None or exon_count <= args.max_exon_count:
                    mRNA_count += 1
                    
                    if args.retain_composition is True:
                        mRNAs_by_exon_count[exon_count].append(mRNA)
                    else:
                        mRNAs.append(mRNA)

    # if you feel like printing a profile
    #for exon_count in mRNAs_by_exon_count:
    #    print("DEBUG: exons:{0}\tcount:{1}".format( exon_count, len(mRNAs_by_exon_count[exon_count]) ) )

    # sanity check on the number of available mRNAs
    if (args.training_set_size + args.evaluation_set_size) > mRNA_count:
        raise Exception("ERROR: acceptable mRNA count ({0}) is less than combined training_set_size ({1}) and evaluation_set_size ({2}) options".format(mRNA_count, args.training_set_size, args.evaluation_set_size) )

    training_mRNAs = list()
    evaluation_mRNAs = list()
    
    if args.retain_composition is True:
        print("DEBUG: retaining composition")
        pass
    else:
        training_mRNAs = random.sample( mRNAs, args.training_set_size )
        unselected_mRNAs = list(set(mRNAs) & set(set(mRNAs) ^ set(training_mRNAs)))
        evaluation_mRNAs = random.sample( unselected_mRNAs, args.evaluation_set_size )

    export_mRNAs_to_file(training_mRNAs, args.output_training_file)
    export_mRNAs_to_file(evaluation_mRNAs, args.output_evaluation_file)
    
    
def export_mRNAs_to_file( mRNAs, f ):
    fh = open(f, 'wt')
    fh.write("##gff-version 3\n")

    mRNA_ids = dict()
    genes_to_print = list()

    for mRNA in mRNAs:
        mRNA_ids[mRNA.id] = 1
        if mRNA.parent not in genes_to_print:
            genes_to_print.append(mRNA.parent)

    for gene in sorted(genes_to_print):
        mRNAs_to_keep = []
        original_mRNAs = []

        for mRNA in gene.mRNAs():
            original_mRNAs.append(mRNA)
            
            if mRNA.id in mRNA_ids:
                mRNAs_to_keep.append(mRNA)

        gene.children['mRNA'] = mRNAs_to_keep
        gene.print_as(fh=fh, source='SOURCE', format='gff3')
        gene.children['mRNA'] = original_mRNAs

    
if __name__ == '__main__':
    main()







