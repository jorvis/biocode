#!/usr/bin/env python3

'''
TYPE 1 - BEST. Has all of the following characteristics:
    1. GM-ES prediction identical to a CEGMA prediction.
    2. Identical to an AAT alignment of at least one C. muris or C. parvum homologue.
    3. Intron-exon boundaries of model supported by assembled transcripts, with start
       & stop coordinates contained within the first and last exons of the transcript,
       respectively. Non-CDS portion does not have to be identical, i.e. allowing for UTRs.

    Example: http://bit.ly/1irtqwZ
    Example: http://bit.ly/1emZ39O # not the best example, because while technically true,
             the terminal exon of the Cufflinks transcript it ridiculously long.

Type 2 - BEST. Has all of the following characteristics:
    1. GM-ES prediction identical to a CEGMA prediction.
    2. NOT identical to an AAT alignment of at least one C. muris or C. parvum homologue.
    3. Intron-exon boundaries of model supported by assembled transcripts, with start &
       stop coordinates contained within the first and last exons of the transcript,
       respectively. Non-CDS portion does not have to be identical, i.e. allowing for UTRs.

    Example:
    http://bit.ly/1g43oCO

TYPE 2 - BETTER. Has all of the following characteristics:
    1. No match to CEGMA model. By no match, I mean either no CEGMA model exists OR one
       does, but it does not have same coordinates.
    2. GM-ES model identical to an AAT alignment of at least one C. muris or C. parvum
       homologue.
    3. Intron-exon boundaries of model supported by assembled transcripts, with start &
       stop coordinates contained within the first and last exons of the transcript,
       respectively.
    
    Example:
    http://bit.ly/1e8CsTG

TYPE 3 - STILL BETTER. Has all of the following characteristics:
    1. CEGMA model (a case where we are using CEGMA model rather than GM-ES model)
    2. Does not match GM-ES model (imperfect match or no GMES model)
    3. CEGMA model identical to an AAT alignment of at least one C. muris or C. parvum
       homologue.
    4. Intron-exon boundaries of CEGMA model supported by assembled transcripts, with
       start & stop coordinates contained within the first and last exons of the transcript,
       respectively. 

    Example:
    http://bit.ly/Q0ZG39
    http://bit.ly/1gJ6Pe6

Author:  Joshua Orvis
'''

from biocode import gff


def main():

    gm_es_file = 'genemark_hmm.gff3'
    cegma_file = 'output.cegma.gff3'
    #aat_file = 'bail_training_genes.aat.1500maxintron.80percid.gff3'
    aat_file = 'aat.bail_hominis_filtered_training.gff3'
    #aat_file = 'aat.merged.gff3'
    

    print("INFO: parsing Genemark-ES data")
    (assemblies, gm_es_features) = gff.get_gff3_features(gm_es_file)
    gm_es_genes = get_genes_from_dict(gm_es_features)
    print("\tINFO: Got {0} Genemark-ES genes".format(len(gm_es_genes)))

    print("INFO: parsing CEGMA data")
    (assemblies, cegma_features) = gff.get_gff3_features(cegma_file, assemblies=assemblies)
    cegma_genes = get_genes_from_dict(cegma_features)
    print("\tINFO: Got {0} CEGMA genes".format(len(cegma_genes)))

    print("INFO: parsing AAT results")
    (assemblies, aat_muris_features) = gff.get_gff3_features(aat_file, assemblies=assemblies)
    aat_genes = get_genes_from_dict(aat_muris_features)
    print("\tINFO: Got {0} AAT 'genes'".format(len(aat_genes)))

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

    genemark_aat_shared_genes = list()
    gmes_aat_fh = open('gmes_aat.shared.ids', 'wt')

    for gm_es_gene in gm_es_genes:
        for aat_gene in aat_genes:
            if gm_es_gene.shares_exon_structure_with( thing=aat_gene, stop_tolerant=True ) == True:
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
            if cegma_gene.shares_exon_structure_with( thing=aat_gene, stop_tolerant=True ) == True:
                cegma_with_aat_not_gm_es.append(cegma_gene)
                cegma_aat_nogmes_fh.write("{0}\n".format(cegma_gene.id))
                break
            
    print("{0} CEGMA genes had no GeneMark-ES match but did have an AAT one".format(len(cegma_with_aat_not_gm_es)) )


    

        

def get_genes_from_dict( features ):
    genes = list()

    for feat_id in features:
        feat = features[feat_id]
        if feat.__class__.__name__ == 'Gene':
            genes.append(feat)

    return genes
    

if __name__ == '__main__':
    main()
