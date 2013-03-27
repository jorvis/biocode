#!/usr/bin/env python3.2

import argparse
import os
import fileinput
import biothings
import biocodegff
import sys

def interface():
    parser = argparse.ArgumentParser( description='Put a description of your script here')
    parser.add_argument('-a1', '--annotation_1',type=str, required=True,
		      help='The first annotation file.')
    parser.add_argument('-a2', '--annotation_2',type=str, required=False,
		      help='The second annotation file.')

    ## output file to be written
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Path to an output directory' )
    parser = parser.parse_args()
    return parser


def process_files(args):
    (assemblies_1, features_1) = biocodegff.get_gff3_features(args.annotation_1)
    (assemblies_2, features_2) = biocodegff.get_gff3_features(args.annotation_2)


    a_exons = set()                                    ## Set contains only uniq exons from known annotation, since multiple same exons can appear in a gff file.  
    p_exons = set()                                    ## For predicted annotation

    a_gene = set()
    p_gene = set()

    a_mrna = set()
    p_mrna = set()

    exon_pred_all = set()
    gene_true = set()
    mrna_true = set()

    chr = []
    
    for asm_id in assemblies_1:                                                                                     ## Iterate through each chromosome from the known ref annotation

        assembly_1 = assemblies_1[asm_id]
        assembly_2 = assemblies_2.get(asm_id,-1)                                                                    ## Find that chromosome in the predicted gff file
        genes_1 = assembly_1.genes()                                                                                ## All genes from known annotation
        
        for gene_1 in genes_1 :                                                                                     ## Add unique gene, mrna , exon features from known annotation to get each known feature total count 
            gene_1_loc = gene_1.location_on(assembly_1)
            cord = asm_id + str(gene_1_loc.fmin) + ":" + str(gene_1_loc.fmax)+ ":"  + str(gene_1_loc.strand)        ## Use chromosome id+start+stop+strand as a string to determine uniqueness.
            a_gene.add(cord)
            
            for mrna_1 in gene_1.mRNAs() :
                mrna_1_loc = mrna_1.location_on(assembly_1)
                cord = asm_id + str(mrna_1_loc.fmin) + ":" + str(mrna_1_loc.fmax) + ":" + str(mrna_1_loc.strand)
                a_mrna.add(cord)

                for exon_1 in mrna_1.exons() :
                    exon_1_loc = exon_1.location_on(assembly_1)
                    cord = asm_id + str(exon_1_loc.fmin) + ":" + str(exon_1_loc.fmax) + ":" + str(exon_1_loc.strand)
                    a_exons.add(cord)
                        
        if (type(assembly_2) is int) :                     ##    If the chromosome is not found in prediected file, move to next chromosome.
            continue
        

        genes_2 = assembly_2.genes()                      ## All genes from predicted annotation.
        chr.append(asm_id)                                ## Append all found chromosome in a list.

        for gene_2 in genes_2 :                           ## Add unique gene, mrna , exon features from predicted annotation to get each predicted feature total count.  
            gene_2_loc = gene_2.location_on(assembly_2)
            cord = asm_id + str(gene_2_loc.fmin) + ":" + str(gene_2_loc.fmax) + ":" +  str(gene_2_loc.strand)
            p_gene.add(cord)
            
            for mrna_2 in gene_2.mRNAs() :
                mrna_2_loc = mrna_2.location_on(assembly_2)
                cord = asm_id + str(mrna_2_loc.fmin) + ":" + str(mrna_2_loc.fmax)+ ":" +  str(mrna_2_loc.strand)
                p_mrna.add(cord)
                
                for exon_2 in mrna_2.exons() :
                    exon_2_loc = exon_2.location_on(assembly_2)
                    cord = asm_id + str(exon_2_loc.fmin) + ":" + str(exon_2_loc.fmax)+ ":" + str(exon_2_loc.strand)
                    p_exons.add(cord)
                    
                    

        for gene_2 in sorted(genes_2) :                                         ## From the predicted feature determine the true once. Iterate through each predicted gene sorted by cordinate
            gene_2_loc = gene_2.location_on(assembly_2)
            cord_g = asm_id + str(gene_2_loc.fmin) + ":" +  str(gene_2_loc.fmax) + ":" + str(gene_2_loc.strand)

            if (cord_g in gene_true) :                                          ## To prevent duplication, check if the feature already exists in the set of truly predicted gene.
                continue

            true_pred_mrna_per_gene = 0
            true_exon_cord = []
            found_exon = 0
            

            for gene_1 in sorted(genes_1) :                                    ## Iterate through each known gene
                gene_1_loc = gene_1.location_on(assembly_1)
                if (gene_1_loc.strand != gene_2_loc.strand) :                  ## Check for strand
                    continue

                if (gene_2_loc.fmax < gene_1_loc.fmin) :                       ## Save time by ignoring the known genes whose start is after predicted gene's end.
                    break
                
                if(gene_2.overlaps_with(gene_1)) :                             ## If known and predicted overlap proceed
                    
                    exon_pred = set();                                         
                    exon_anno = set();

                    for mrna_2 in sorted(gene_2.mRNAs()) :                     ## Add all exons for the predicted and known genes to a set 
                        for exon_2 in sorted(mrna_2.exons()) :
                            exon_2_loc = exon_2.location_on(assembly_2)
                            cord = asm_id + str(exon_2_loc.fmin) + ":" + str(exon_2_loc.fmax) + ":" + str(exon_2_loc.strand) 
                            exon_pred.add(cord)
                            
                                                
                    for mrna_1 in sorted(gene_1.mRNAs()) :
                        for exon_1 in sorted(mrna_1.exons()) :
                            exon_1_loc = exon_1.location_on(assembly_1)
                            cord = asm_id + str(exon_1_loc.fmin) + ":" + str(exon_1_loc.fmax) + ":" + str(exon_1_loc.strand)
                            exon_anno.add(cord)
                            

                    for exon_2 in exon_pred :                                      
                        if (exon_2 in exon_pred_all) :
                            continue
                        for exon_1 in exon_anno :
                                                                                  ## If exon boundaries match, the predicted exon is TRUE...
                             if (exon_1 == exon_2) : 
                                 exon_pred_all.add(exon_2)
                                 true_exon_cord.append(exon_2)                    # Append each true exon per gene in a list  to check for true mrna
                                 found_exon = 1
                                 break
                    
            
            for mrna_2 in sorted(gene_2.mRNAs()) :                                ## Iterate through each predicted mrna , if all of its exon is true , then the predicted mRNA is true
                mrna_2_loc = mrna_2.location_on(assembly_2)
                cord_m = asm_id + str(mrna_2_loc.fmin) + ":" + str(mrna_2_loc.fmax)  + ":" + str(mrna_2_loc.strand)

                if (cord_m in mrna_true) :
                    continue

                if (found_exon == 0) :
                    break

                count = 0
                pred_exon = set()

                for exon_2 in sorted(mrna_2.exons()) :
                    exon_2_loc = exon_2.location_on(assembly_2)
                    cord = asm_id + str(exon_2_loc.fmin) + ":" + str(exon_2_loc.fmax) + ":" +  str(exon_2_loc.strand)
                    if cord in pred_exon :
                        continue
                    pred_exon.add(cord)
                    for true_exon in true_exon_cord :
                        if (cord == true_exon) :
                            count += 1
                            break
                        
                if (len(pred_exon) == count) :
                    mrna_true.add(cord_m)
                    true_pred_mrna_per_gene += 1

            if (true_pred_mrna_per_gene >= 1) :                                  ## If the predicted gene has atleast one true predicted mrna, then the gene is true.
                gene_true.add(cord_g)


    for asm_id in assemblies_2:                                                  ## Iterate through each chromosome from the predicted annotation
        if asm_id not in chr :
            assembly_2 = assemblies_2.get(asm_id,-1)                             ## Find that chromosome in the predicted gff file which is not found in known annotation
            genes_2 = assembly_2.genes()                                         ## Add  genes, mrna, exon features from predicted annotation to total predicted feature set.
        
            for gene_2 in genes_2 :
                gene_2_loc = gene_2.location_on(assembly_2)
                cord = asm_id + str(gene_2_loc.fmin) + ":" + str(gene_2_loc.fmax)  + ":"+ str(gene_2_loc.strand)
                p_gene.add(cord)
            
                for mrna_2 in gene_2.mRNAs() :
                    mrna_2_loc = mrna_2.location_on(assembly_2)
                    cord = asm_id + str(mrna_2_loc.fmin) + ":" + str(mrna_2_loc.fmax) + ":" + str(mrna_2_loc.strand)
                    p_mrna.add(cord)

                    for exon_2 in mrna_2.exons() :
                        exon_2_loc = exon_2.location_on(assembly_2)
                        cord = asm_id + str(exon_2_loc.fmin) + ":" + str(exon_2_loc.fmax) + ":" + str(exon_2_loc.strand)
                        p_exons.add(cord)
   


    #Calculate SN/SP for exons 
    annotated_exon = len(a_exons)
    predicted_exon = len(p_exons)
    true_pred_exon = len(exon_pred_all)
    
    exon_sn = (true_pred_exon/annotated_exon) * 100                                 
    exon_sp = (true_pred_exon/predicted_exon) * 100

    #Calculate SN/SP for transcript
    
    annotated_mrna = len(a_mrna)
    predicted_mrna = len(p_mrna)
    true_pred_mrna = len(mrna_true)
    
    mrna_sn = (true_pred_mrna/annotated_mrna) * 100
    mrna_sp = (true_pred_mrna/predicted_mrna) * 100
       
    #Calculate SN/SP for genes 

    annotated_gene = len(a_gene)
    predicted_gene = len(p_gene)
    true_pred_gene = len(gene_true)
    
    gene_sn = (true_pred_gene/annotated_gene) * 100                                 
    gene_sp = (true_pred_gene/predicted_gene) * 100

    
    out_file = args.output_dir + '/summary.txt'
    if not (os.path.exists(out_file)) :
        sys.exit("Directory does not exist.")
    fout = open(out_file,'w')

    fout.write("Feature\tKnown\tPredicted\tTrue_Predicted\tSN\tSP\n")
    fout.write("Gene\t"+str(annotated_gene)+"\t"+str(predicted_gene)+"\t"+str(true_pred_gene)+"\t"+str(gene_sn)+"\t"+str(gene_sp)+"\n")
    fout.write("mRNA\t"+str(annotated_mrna)+"\t"+str(predicted_mrna)+"\t"+str(true_pred_mrna)+"\t"+str(mrna_sn)+"\t"+str(mrna_sp)+"\n")
    fout.write("Exon\t"+str(annotated_exon)+"\t"+str(predicted_exon)+"\t"+str(true_pred_exon)+"\t"+str(exon_sn)+"\t"+str(exon_sp)+"\n")
 
    
   
if __name__ == '__main__':
    args = interface()
    process_files(args)






