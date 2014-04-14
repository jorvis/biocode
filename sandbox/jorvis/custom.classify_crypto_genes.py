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

Author:  Joshua Orvis
'''

import os
import biocodegff
import biocodeutils


def main():

    gm_es_file = 'genemark_hmm.gff3'
    cegma_file = 'output.cegma.gff3'
    transcript_file = 'expression_data.gff3'

    debug_gm_es_gene = '2460_g'

    print("INFO: parsing Genemark-ES data")
    (assemblies, gm_es_features) = biocodegff.get_gff3_features( gm_es_file )

    print("INFO: parsing CEGMA data")
    (assemblies, cegma_features) = biocodegff.get_gff3_features( cegma_file, assemblies=assemblies )

    print("INFO: parsing expression data (Trinity, Cufflinks, GMAP cDNAs)")
    (assemblies, transcript_features) = biocodegff.get_gff3_features( transcript_file, assemblies=assemblies)
    
    #biocodeutils.add_assembly_fasta(assemblies, args.masked_fasta)

    genemark_cegma_shared_genes = list()

    for gm_es_feat_id in gm_es_features:
        gm_es_feat = gm_es_features[gm_es_feat_id]

        if gm_es_feat.__class__.__name__ != 'Gene':
            continue

        for cegma_feat_id in cegma_features:
            cegma_feat = cegma_features[cegma_feat_id]
            cegma_feat_loc = cegma_feat.location()

            if cegma_feat.__class__.__name__ != 'Gene':
                continue

            if gm_es_feat.has_same_coordinates_as( thing=cegma_feat ):

                if gm_es_feat.shares_exon_structure_with( thing=cegma_feat ) == True:
                    #print("Found complete agreement between GeneMark-ES:{0} and CEGMA:{1} on {2}, {3}-{4}".format(gm_es_feat_id, cegma_feat_id, \
                    #  cegma_feat_loc.on.id, cegma_feat_loc.fmin, cegma_feat_loc.fmax))
                    genemark_cegma_shared_genes.append(gm_es_feat)
                    break

        if gm_es_feat.id == debug_gm_es_gene:
            if gm_es_feat in genemark_cegma_shared_genes:
                print("GM-ES gene {0} has a CEGMA match".format(gm_es_feat.id))
            else:
                print("GM-ES gene {0} doesn't have a CEGMA match".format(gm_es_feat.id))
                
    print("\n{0} genes were shared perfectly between Genemark-ES and CEGMA".format(len(genemark_cegma_shared_genes)) )

    gm_cegma_expression_shared_genes = list()

    for shared_gene in genemark_cegma_shared_genes:
        for tf_id in transcript_features:
            tf = transcript_features[tf_id]

            if tf.__class__.__name__ != 'Gene':
                continue

            if shared_gene.shares_CDS_structure_with( tf ):
                #print("Good gene: {0}".format(shared_gene.id))
                gm_cegma_expression_shared_genes.append( shared_gene )
                break

        if shared_gene.id == debug_gm_es_gene:
            if shared_gene in gm_cegma_expression_shared_genes:
                print("GM-ES gene {0} has an expression match".format(shared_gene.id))
            else:
                print("GM-ES gene {0} doesn't have an expression match".format(shared_gene.id))


    print("\n{0} genes were shared perfectly between Genemark-ES and CEGMA with expression support".format(len(gm_cegma_expression_shared_genes)) )

    

if __name__ == '__main__':
    main()






