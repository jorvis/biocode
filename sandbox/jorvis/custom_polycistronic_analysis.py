#!/usr/bin/env python3.2

import argparse

from biocode import things


## tracked in JIRA-247 - Initially written for a specific purpose but will abstract out if useful.
# currently only cares about gene and mRNA features

def main():
    parser = argparse.ArgumentParser( description='Script for reporting of possible polycistronic genes transcripts based on a reference annotation and RNA-seq transcript assemblies')

    ## output file to be written
    parser.add_argument('-r', '--reference_file', type=str, required=True, help='GFF3 file of a reference annotation' )
    parser.add_argument('-q', '--query_file', type=str, required=True, help='GFF3 file with alternative annotation (such as an RNA-seq assemby)' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    (ref_assemblies, ref_genes) = parse_gff3( args.reference_file )
    (qry_assemblies, qry_genes) = parse_gff3( args.query_file )

    for ref_gene_id in ref_genes:
        ref_gene = ref_genes[ref_gene_id]

        for qry_gene_id in qry_genes:
            qry_gene = qry_genes[qry_gene_id]

            overlap = ref_gene.overlap_with(gene=qry_gene)

            if overlap:
                print("DEBUG: {0} and {1} appear to overlap".format(ref_gene_id, qry_gene_id) )


    

def parse_gff3(gff3_file):

    assemblies = dict()
    genes = dict()
    
    for line in open(gff3_file):
        cols = line.split("\t")

        if len(cols) != 9:
            continue

        mol_id = cols[0]
        
        ## initialize this assembly if we haven't seen it yet
        if mol_id not in assemblies:
            assemblies[mol_id] = things.Assembly(id=mol_id)

        current_assembly = assemblies[mol_id]
        rfmin = int(cols[3]) - 1
        rfmax = int(cols[4])
        rstrand = None
        feat_id = column_9_value(cols[8], 'ID')
        #print("Processing feature: ({0})".format(feat_id))

        if cols[6] == '-':
            strand = -1
        elif cols[6] == '+':
            strand = 1
        else:
            strand = 0

        if cols[2] == 'gene':
            gene = things.Gene(id=feat_id)
            gene.locate_on(assembly=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            genes[feat_id] = gene
        
        elif cols[2] == 'mRNA':
            mRNA = things.mRNA(id=feat_id)
            mRNA.locate_on(assembly=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            parent_id = column_9_value(cols[8], 'Parent')
            genes[parent_id].add_mRNA( mRNA )

    return (assemblies, genes)



if __name__ == '__main__':
    main()




