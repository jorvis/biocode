#!/usr/bin/env python3

"""

check_overlaps.py: Given two GFF files (reference,query), print to STDOUT a list of all of the genes that overlap 1:1 by CDS. This script can be used to create a correspondence table for assigning locus tags. 

Assumptions:
- One isoform per gene
- Biocode-friendly GFF formats


"""

import argparse

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-r', '--reference', type=str, required=True, help='Reference GFF3' )
    parser.add_argument('-q', '--query', type=str, required=True, help='Query GFF3' )
    args = parser.parse_args()

    print("INFO: parsing reference features\n")
    (assemblies, ref_features) = gff.get_gff3_features(args.reference)

    print("INFO: parsing query features\n")
    (assemblies, qry_features) = gff.get_gff3_features(args.query, assemblies=assemblies)

    ref_genes = get_genes_from_dict( ref_features )
    qry_genes = get_genes_from_dict( qry_features )
    ref_gene_one_qry_overlap = {}
    qry_gene_one_ref_overlap = {}
#Find all of the query genes that overlap the reference gene at all
    for ref_gene in sorted(ref_genes):
        ref_loc = ref_gene.location()
        num_qry_overlaps = {} #keep track of the number of query RNAs that have at least one CDS that overlaps
        qry_to_ref = {} #keep track of the query:ref relationships for printing out later
        for qry_gene in sorted(qry_genes):
            qry_loc = qry_gene.location()
            if ref_gene.overlaps_with( qry_gene ):
                for ref_RNA in ref_gene.RNAs():
                    for qry_RNA in qry_gene.RNAs():
                        for ref_CDS in ref_RNA.CDSs(): 
                            ref_CDS_loc = ref_CDS.location()
                            for qry_CDS in qry_RNA.CDSs():
                                qry_CDS_loc = qry_CDS.location()
                                if ref_CDS.overlaps_with( qry_CDS ) and qry_CDS_loc.strand is ref_CDS_loc.strand: #Does the ref CDS overlap the query CDS?
                                    num_qry_overlaps[qry_gene.id] = 1 #If so, add the qry_gene.id to the list of overlaps
                                    qry_to_ref[qry_gene.id] = ref_gene.id #Also, keep track of the query:ref relationships
        #Store all of the reference genes that overlap only a single query gene
        for qry_gene in num_qry_overlaps:
            if len(num_qry_overlaps) == 1:
                ref_gene_one_qry_overlap[qry_to_ref[qry_gene]] = qry_gene
                #print(str(qry_gene) + "\t" + str(qry_to_ref[qry_gene]))$
    #Now do the same thing finding all ref overlaps for each query gene
    for qry_gene in sorted(qry_genes):
        qry_loc = qry_gene.location()
        num_ref_overlaps = {}
        ref_to_qry = {}
        for ref_gene in sorted(ref_genes):
            ref_loc = ref_gene.location()
            if qry_gene.overlaps_with( ref_gene ):
                for qry_RNA in qry_gene.RNAs():
                    for ref_RNA in ref_gene.RNAs():
                        for qry_CDS in qry_RNA.CDSs():
                            qry_CDS_loc = qry_CDS.location()
                            for ref_CDS in ref_RNA.CDSs():
                                ref_CDS_loc = ref_CDS.location()
                                if qry_CDS.overlaps_with( ref_CDS ) and qry_CDS_loc.strand is ref_CDS_loc.strand:
                                    num_ref_overlaps[ref_gene.id] = 1
                                    ref_to_qry[ref_gene.id] = qry_gene.id
        #Store all of the wry genes that overlap only a single reference gene
        for ref_gene in num_ref_overlaps:
            if len(num_ref_overlaps) == 1:
                qry_gene_one_ref_overlap[ref_to_qry[ref_gene]] = ref_gene
#Find all of the reference genes with only one query overlap and vice versa and print them out
    for qry_gene_id in qry_gene_one_ref_overlap:
       for ref_gene_id in ref_gene_one_qry_overlap:
            if qry_gene_id is ref_gene_one_qry_overlap[ref_gene_id] and ref_gene_id is qry_gene_one_ref_overlap[qry_gene_id]:
                print(qry_gene_id + "\t" + ref_gene_id)


def get_genes_from_dict( features ):
    genes = list()

    for feat_id in features:
        feat = features[feat_id]
        if feat.__class__.__name__ == 'Gene':
            genes.append(feat)

    return genes



if __name__ == '__main__':
    main()


