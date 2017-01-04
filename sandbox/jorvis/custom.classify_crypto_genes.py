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
    transcript_file = 'expression_data.gff3'
    aat_muris_file = 'cmuris.aat.gff3'
    aat_parvum_file = 'cparvum.aat.gff3'

    html_out_file = 'gene_classifications.html'
    html_out = open(html_out_file, 'wt')

    # for each gene, how many flanking bases should be shown on either side?
    flanking_bases = 1000

    type1_best = list()
    type2_best = list()
    type2_better = list()
    type3_still_better = list()

    print("INFO: parsing Genemark-ES data")
    (assemblies, gm_es_features) = gff.get_gff3_features(gm_es_file)
    gm_es_genes = get_genes_from_dict(gm_es_features)
    print("\tINFO: Got {0} Genemark-ES genes".format(len(gm_es_genes)))

    print("INFO: parsing CEGMA data")
    (assemblies, cegma_features) = gff.get_gff3_features(cegma_file, assemblies=assemblies)
    cegma_genes = get_genes_from_dict(cegma_features)
    print("\tINFO: Got {0} CEGMA genes".format(len(cegma_genes)))

    print("INFO: parsing expression data (Trinity, Cufflinks, GMAP cDNAs)")
    (assemblies, transcript_features) = gff.get_gff3_features(transcript_file, assemblies=assemblies)
    transcript_genes = get_genes_from_dict(transcript_features)
    print("\tINFO: Got {0} expression 'genes'".format(len(transcript_genes)))

    print("INFO: parsing AAT results (C. muris)")
    (assemblies, aat_muris_features) = gff.get_gff3_features(aat_muris_file, assemblies=assemblies)
    aat_muris_genes = get_genes_from_dict(aat_muris_features)
    print("\tINFO: Got {0} AAT (C. muris) 'genes'".format(len(aat_muris_genes)))

    print("INFO: parsing AAT results (C. parvum)")
    (assemblies, aat_parvum_features) = gff.get_gff3_features(aat_parvum_file, assemblies=assemblies)
    aat_parvum_genes = get_genes_from_dict(aat_parvum_features)
    print("\tINFO: Got {0} AAT (C. parvum) 'genes'".format(len(aat_parvum_genes)))
    
    #biocodeutils.add_assembly_fasta(assemblies, args.masked_fasta)

    genemark_cegma_shared_genes = list()

    for gm_es_gene in gm_es_genes:
        for cegma_gene in cegma_genes:
            if gm_es_gene.has_same_coordinates_as( thing=cegma_gene ):
                if gm_es_gene.shares_exon_structure_with( thing=cegma_gene ) == True:
                    genemark_cegma_shared_genes.append(gm_es_gene)
                    break

    print("\n{0} genes were shared perfectly between Genemark-ES and CEGMA".format(len(genemark_cegma_shared_genes)) )
    #############################################################################

    genemark_aat_shared_genes = list()

    for gm_es_gene in gm_es_genes:
        for aat_gene in aat_muris_genes:
            if gm_es_gene.shares_exon_structure_with( thing=aat_gene ) == True:
                genemark_aat_shared_genes.append(gm_es_gene)
                break

        if gm_es_gene not in genemark_aat_shared_genes:
            for aat_gene in aat_parvum_genes:
                if gm_es_gene.shares_exon_structure_with( thing=aat_gene ) == True:
                    genemark_aat_shared_genes.append(gm_es_gene)
                    break

    print("{0} Genemark-ES genes had an exact AAT match".format(len(genemark_aat_shared_genes)) )    

    ##############################################################################
    cegma_not_matching_gm_es = list()
    
    for cegma_gene in cegma_genes:
        match_found = False

        for gm_es_gene in gm_es_genes:
            if cegma_gene.has_same_coordinates_as( thing=gm_es_gene ):
                if cegma_gene.shares_exon_structure_with( thing=gm_es_gene ) == True:
                    match_found = True
                    break

        if match_found == False:
            cegma_not_matching_gm_es.append(cegma_gene)

    print("{0} CEGMA genes don't have a structural match to a Genemark-ES one".format(len(cegma_not_matching_gm_es)) )

    #############################################################################
    gm_expression_shared_genes = list()
    
    for gm_es_gene in gm_es_genes:
        for tf in transcript_genes:
            if gm_es_gene.shares_CDS_structure_with( tf ):
                gm_expression_shared_genes.append(gm_es_gene)
                break
    
    print("{0} Genemark-ES genes had an exact expression match".format(len(gm_expression_shared_genes)) )
    
    #############################################################################

    gm_cegma_expression_shared_genes = list()

    for shared_gene in genemark_cegma_shared_genes:
        if shared_gene in gm_expression_shared_genes:
            gm_cegma_expression_shared_genes.append( shared_gene )

    print("{0} genes were shared perfectly between Genemark-ES and CEGMA with expression support".format(len(gm_cegma_expression_shared_genes)) )
    ##############################################################################
    
    gm_cegma_expression_aat_shared_genes = list()

    for shared_gene in gm_cegma_expression_shared_genes:
        if shared_gene in genemark_aat_shared_genes:
            gm_cegma_expression_aat_shared_genes.append(shared_gene)
        else:
            type2_best.append(shared_gene)
        
    for gene in gm_cegma_expression_aat_shared_genes:
        type1_best.append(gene)

    print("{0} genes were shared with Genemark-ES, CEGMA, expression, AAT support".format(len(gm_cegma_expression_aat_shared_genes)) )
    ##############################################################################

    for gm_es_gene in gm_es_genes:
        if gm_es_gene not in genemark_cegma_shared_genes:
            if gm_es_gene in gm_expression_shared_genes:
                type2_better.append(gm_es_gene)
    
    ##############################################################################
    cegma_expression_shared_genes = list()
    
    for cegma_gene in cegma_genes:
        for tf in transcript_genes:
            if cegma_gene.shares_CDS_structure_with( tf ):
                cegma_expression_shared_genes.append(cegma_gene)
                break
    
    print("{0} CEGMA genes had an exact expression match".format(len(cegma_expression_shared_genes)) )


    ##############################################################################
    
    cegma_not_gmes_with_aat = list()

    for cegma_gene in cegma_genes:
        if cegma_gene in cegma_not_matching_gm_es:
            
            for aat_gene in aat_muris_genes:
                if cegma_gene.shares_exon_structure_with( thing=aat_gene ) == True:
                    cegma_not_gmes_with_aat.append(cegma_gene)
                    break

            if cegma_gene not in cegma_not_gmes_with_aat:
                for aat_gene in aat_parvum_genes:
                    if cegma_gene.shares_exon_structure_with( thing=aat_gene ) == True:
                        cegma_not_gmes_with_aat.append(gm_es_gene)
                        break

            if cegma_gene in cegma_not_gmes_with_aat:
                if cegma_gene in cegma_expression_shared_genes:
                    type3_still_better.append(cegma_gene)
                

    print("TYPE 1 - BEST: {0}".format(len(type1_best)) )
    print("TYPE 2 - BEST: {0}".format(len(type2_best)) )
    print("TYPE 2 - BETTER: {0}".format(len(type2_better)) )
    print("TYPE 3 - STILL BETTER: {0}".format(len(type3_still_better)) )
    
    html_out.write("<!doctype html>\n")
    html_out.write("<html lang=\"en\">\n")
    html_out.write("<head><meta charset=\"utf-8\"><title>Gene classification</title></head>\n")
    html_out.write("<body>\n")

    html_out.write("<h3>Type 1 - Best ({0})</h3>\n".format(len(type1_best)))
    print_gene_list(html_out, type1_best, flanking_bases)

    html_out.write("<h3>Type 2 - Best ({0})</h3>\n".format(len(type2_best)))
    print_gene_list(html_out, type2_best, flanking_bases)

    html_out.write("<h3>Type 2 - Better ({0})</h3>\n".format(len(type2_better)))
    print_gene_list(html_out, type2_better, flanking_bases)

    html_out.write("<h3>Type 3 - Still better ({0})</h3>\n".format(len(type3_still_better)))
    print_gene_list(html_out, type3_still_better, flanking_bases)

    html_out.write("</body>\n")
    html_out.write("</html>\n")

    #print("They are:\n")
    #for gene in gm_cegma_expression_aat_shared_genes:
    #    print("\t{0}".format(gene.id))
        

def print_gene_list( outfh, genes, flank ):
    base = 'http://jbrowse.igs.umaryland.edu/c_hominis_TU502'
    # /?loc=ctg7180000000305%3A76133..78659
    track_string = '&tracks=DNA%2CReference%20CDS%20(gmap)%2CReference%20transcripts%20(gmap)%2CGenemark-ES%2CCEGMA%2CTrinity%20(gmap)%2CCufflinks%2CAAT%20(C.%20muris)%2CAAT%20(C.%20parvum)%2CNCBI%20fl-cDNAs%20(gmap)%2Ccufflinks_bam%2Ccufflinks_UKH1_bam%2Ctophat_sense%2Ctophat_antisense&highlight='

    outfh.write("<table>\n")

    for gene in genes:
        loc = gene.location()
        
        outfh.write("<tr>\n")
        outfh.write("  <td><a href=\"{1}/?loc={2}%3A{3}..{4}{5}\">{0}</a></td>\n".format(gene.id, base, loc.on.id, loc.fmin + 1 - flank, loc.fmax + flank, track_string))
        outfh.write("  <td>{0}</td>\n".format(loc.on.id))
        outfh.write("  <td>{0}</td>\n".format(loc.fmin + 1))
        outfh.write("  <td>{0}</td>\n".format(loc.fmax))
        outfh.write("<tr>\n")


def get_genes_from_dict( features ):
    genes = list()

    for feat_id in features:
        feat = features[feat_id]
        if feat.__class__.__name__ == 'Gene':
            genes.append(feat)

    return genes
    

if __name__ == '__main__':
    main()
