#!/usr/bin/env python3

"""
This script will take two generic gff3 files that contain gene structure
annotation and compare them. When a model from file A has identical coordinates
to a model contained in file B, this is considered an identity match.

Every part of the gene model is compared.  For each gene, the following
coordinates must match:  gene, mRNA, CDS, exon

USAGE:  report_basic_gff_model_agreement.py --ref fileA --qry fileB --output_base comparison

This script produces four output files, with a prefix defined by --output_base:

- $prefix.matches: a tab-delimited file with ids for models with identity matches,
  where column 1 is the id from gff3 file A and column 2 is the
  id from gff3 file B
- $prefix.summary: a summary text file with counts of agreement and disagreement between
  the two files (total identical models, total models in A not in B, total
  models in B not in A); 
- a gff3 file containing only the subset of models from gff3 file A
  that have identity matches to gff3 file B;
- a gff3 file containing only the subset of models from gff3 file B
  that have identity matches to gff3 file A.

Author: Joshua Orvis
"""

import argparse

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Basic comparison of two GFF3 files')

    ## output file to be written
    parser.add_argument('-r', '--ref', type=str, required=True, help='Path to the reference GFF3 file' )
    parser.add_argument('-q', '--qry', type=str, required=True, help='Path to the query GFF3 file' )
    parser.add_argument('-o', '--output_base', type=str, required=True, help='Base name/path of the output files to be created' )
    args = parser.parse_args()

    (assemblies, ref_features) = gff.get_gff3_features(args.ref)
    ref_genes = get_genes_from_dict(ref_features)
    
    (assemblies, qry_features) = gff.get_gff3_features(args.qry, assemblies=assemblies)
    qry_genes = get_genes_from_dict(qry_features)

    ref_matches_found = dict()
    qry_matches_found = dict()

    for ref_gene in ref_genes:
        for qry_gene in qry_genes:
            if ref_gene.has_same_coordinates_as( thing=qry_gene ) and \
               ref_gene.shares_exon_structure_with( thing=qry_gene ) and \
               ref_gene.shares_CDS_structure_with( thing=qry_gene ):

                ref_matches_found[ref_gene.id] = qry_gene.id
                qry_matches_found[qry_gene.id] = ref_gene.id

    # open our output files
    out_matches = open("{0}.matches".format(args.output_base), 'wt')
    out_summary = open("{0}.summary".format(args.output_base), 'wt')

    print("INFO: {0}/{1} reference genes had a match to a qry gene".format( len(ref_matches_found), len(ref_genes) ))
    print("INFO: {0}/{1} qry genes had a match to a reference gene".format( len(qry_matches_found), len(qry_genes) ))

    for ref_gene_id in ref_matches_found:
        out_matches.write("{0}\t{1}\n".format(ref_gene_id, ref_matches_found[ref_gene_id]))

    out_summary.write("Reference\t{0}\n".format(args.ref) )
    out_summary.write("Query\t{0}\n".format(args.ref) )
    out_summary.write("Total identical models (with respect to reference)\t{0}\n".format(len(ref_matches_found)))
    out_summary.write("Models in REF not in QRY\t{0}\n".format( len(ref_genes) - len(ref_matches_found) ))
    out_summary.write("Models in QRY not in REF\t{0}\n".format( len(qry_genes) - len(qry_matches_found) ))


## This should be replaced with a more general filter_by_type added to biocode utils
def get_genes_from_dict( features ):
    genes = list()

    for feat_id in features:
        feat = features[feat_id]
        if feat.__class__.__name__ == 'Gene':
            genes.append(feat)

    return genes

if __name__ == '__main__':
    main()







