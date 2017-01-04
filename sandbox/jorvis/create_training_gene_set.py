#!/usr/bin/env python3

'''
This script was born from the need to auto-generate a set of training genes given
a selection of pre-generated evidence from tools like GeneMark-ES, CEGMA, AAT and
transcript evidence aligned with GMAP.

Initially written to support the annotation of 40+ mucormycosis-causing genomes,
this may not be applicable to other species, but I've tried to make it configurable
and document what it does.  Read on, and use at your own discretion.

INPUT



GENE SELECTION (in order of decreasing priority)
74 genes with GeneMark-ES, CEGMA and AAT agreement
28 CEGMA genes had no GeneMark-ES match but did have an AAT one
137 genes were shared perfectly between Genemark-ES and CEGMA
664 Genemark-ES genes had an exact AAT match


OUTPUT

Also writes IDs of each comparison in the following files:
  - gmes_aat.shared.ids
  - gmes_cegma.shared.ids
  - gmes_aat_cegma.shared.ids
  - cegma_aat_nogmes.shared.ids
  - training.all.ids


Author:  Joshua Orvis
'''

import argparse

from biocode import gff


def main():

    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-g', '--genemark', type=str, required=True, help='Path to the results from GeneMark-ES' )
    parser.add_argument('-c', '--cegma', type=str, required=True, help='Path to the results from CEGMA, converted to GFF3' )
    parser.add_argument('-a', '--aat', type=str, required=True, help='Path to the results from AAT, converted to GFF3' )
    parser.add_argument('-e', '--expression', type=str, required=False, help='Any expression data aligned using GMAP (in gff3_gene mode)' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-m', '--max_genes', type=int, required=False, help='Limits gene IDs exported to the top N by strongest evidence class' )
    args = parser.parse_args()

    print("INFO: parsing Genemark-ES data")
    (assemblies, gm_es_features) = gff.get_gff3_features(args.genemark)
    gm_es_genes = get_genes_from_dict(gm_es_features)
    print("\tINFO: Got {0} Genemark-ES genes".format(len(gm_es_genes)))

    print("INFO: parsing CEGMA data")
    (assemblies, cegma_features) = gff.get_gff3_features(args.cegma, assemblies=assemblies)
    cegma_genes = get_genes_from_dict(cegma_features)
    print("\tINFO: Got {0} CEGMA genes".format(len(cegma_genes)))

    print("INFO: parsing AAT results")
    (assemblies, aat_features) = gff.get_gff3_features(args.aat, assemblies=assemblies)
    aat_genes = get_genes_from_dict(aat_features)
    print("\tINFO: Got {0} AAT 'genes'".format(len(aat_genes)))

    expression_genes = list()
    if args.expression is not None:
        print("INFO: parsing expression results")
        (assemblies, expression_features) = gff.get_gff3_features(args.expression, assemblies=assemblies)
        expression_genes = get_genes_from_dict(expression_features)
        print("\tINFO: Got {0} expression 'genes'".format(len(expression_genes)))

    genemark_cegma_shared_genes = list()
    gmes_cegma_fh = open('gmes_cegma.shared.ids', 'wt')

    for gm_es_gene in gm_es_genes:
        for cegma_gene in cegma_genes:
            if gm_es_gene.has_same_coordinates_as( thing=cegma_gene ):
                if gm_es_gene.shares_exon_structure_with( thing=cegma_gene ) == True:
                    genemark_cegma_shared_genes.append(gm_es_gene)
                    gmes_cegma_fh.write("{0}\n".format(gm_es_gene.id))
                    break

    print("\n{0} genes were shared perfectly between Genemark-ES and CEGMA".format(len(genemark_cegma_shared_genes)) )

    #############################################################################

    genemark_cegma_expression_shared_genes = list()
    gmes_cegma_exp_fh = open('gmes_cegma_exp.shared.ids', 'wt')

    for gm_es_gene in genemark_cegma_shared_genes:
        for exp_gene in expression_genes:
            if gm_es_gene.shares_CDS_structure_with( exp_gene ):
                genemark_cegma_expression_shared_genes.append(gm_es_gene)
                break

    print("{0} genes were shared perfectly between Genemark-ES and CEGMA and expression data".format(len(genemark_cegma_expression_shared_genes)) )

    #############################################################################

    genemark_aat_shared_genes = list()
    gmes_aat_fh = open('gmes_aat.shared.ids', 'wt')

    for gm_es_gene in gm_es_genes:
        for aat_gene in aat_genes:
            if gm_es_gene.shares_exon_structure_with( thing=aat_gene, stop_tolerant=True ) == True:
            #if gm_es_gene.shares_exon_structure_with( thing=aat_gene ) == True:
                genemark_aat_shared_genes.append(gm_es_gene)
                gmes_aat_fh.write("{0}\n".format(gm_es_gene.id))
                break

    print("{0} Genemark-ES genes had an exact AAT match".format(len(genemark_aat_shared_genes)) )    

    ##############################################################################
    cegma_matching_gm_es = list()
    genemark_aat_cegma_shared_genes = list()
    gmes_aat_cegma_fh = open('gmes_aat_cegma.shared.ids', 'wt')
    
    for cegma_gene in cegma_genes:
        match_found = False

        for gm_es_gene in gm_es_genes:
            if cegma_gene.has_same_coordinates_as( thing=gm_es_gene ):
                if cegma_gene.shares_exon_structure_with( thing=gm_es_gene ) == True:
                    match_found = True

                    if gm_es_gene in genemark_aat_shared_genes and gm_es_gene not in genemark_aat_cegma_shared_genes:
                        genemark_aat_cegma_shared_genes.append(gm_es_gene)
                        gmes_aat_cegma_fh.write("{0}\n".format(gm_es_gene.id))
                        
                    break

        if match_found == True:
            cegma_matching_gm_es.append(cegma_gene)

    
    print("{0} genes with GeneMark-ES, CEGMA and AAT agreement".format(len(genemark_aat_cegma_shared_genes)) )
    training_fh = open('training_gene.ids', 'wt')
    
    for gene in genemark_aat_cegma_shared_genes:
        training_fh.write("{0}\n".format(gene.id) )

    ##############################################################################
    cegma_with_aat_not_gm_es = list()
    cegma_aat_nogmes_fh = open('cegma_aat_nogmes.shared.ids', 'wt')
    
    for cegma_gene in cegma_genes:
        if cegma_gene in cegma_matching_gm_es:
            continue

        for aat_gene in aat_genes:
            #if cegma_gene.shares_exon_structure_with( thing=aat_gene, stop_tolerant=True ) == True:
            if cegma_gene.shares_exon_structure_with( thing=aat_gene ) == True:
                cegma_with_aat_not_gm_es.append(cegma_gene)
                cegma_aat_nogmes_fh.write("{0}\n".format(cegma_gene.id))
                break
            
    print("{0} CEGMA genes had no GeneMark-ES match but did have an AAT one".format(len(cegma_with_aat_not_gm_es)) )


    ##############################################################################
    ## now to assemble the results
    training_ids = list()

    # 0. Start with genes shared between GeneMark-ES, CEGMA and expression evidence
    recruit_training_genes( training_ids, genemark_cegma_expression_shared_genes, args.max_genes )
    print("DEBUG: {0} genes after recruitment of GeneMark-ES, CEGMA and expression data".format(len(training_ids)))
    
    # 1. Pull in the genes with shared evidence across GeneMark-ES, CEGMA and AAT
    recruit_training_genes( training_ids, genemark_aat_cegma_shared_genes, args.max_genes )
    print("DEBUG: {0} genes after recruitment of GeneMark-ES, CEGMA and AAT".format(len(training_ids)))

    # 2. Next include those genes 
    recruit_training_genes( training_ids, cegma_with_aat_not_gm_es, args.max_genes )
    print("DEBUG: {0} genes after recruitment of CEGMA + AAT without GM-ES".format(len(training_ids)))

    recruit_training_genes( training_ids, genemark_cegma_shared_genes, args.max_genes )
    print("DEBUG: {0} genes after recruitment of GeneMark-ES + CEGMA".format(len(training_ids)))

    recruit_training_genes( training_ids, genemark_aat_shared_genes, args.max_genes )
    print("DEBUG: {0} genes after recruitment of GeneMark-ES + AAT".format(len(training_ids)))

    output_list_fh = open(args.output_file, 'wt')
    for training_id in training_ids:
        output_list_fh.write("{0}\n".format(training_id))


        
def recruit_training_genes( ids, genes, max_c ):
    for gene in genes:
        if max_c is not None and len(ids) >= max_c:
            break
        else:
            if gene.id not in ids:
                ids.append(gene.id)

def get_genes_from_dict( features ):
    genes = list()

    for feat_id in features:
        feat = features[feat_id]
        if feat.__class__.__name__ == 'Gene':
            genes.append(feat)

    return genes
    

if __name__ == '__main__':
    main()
