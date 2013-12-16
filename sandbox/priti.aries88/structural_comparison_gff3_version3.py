#!/usr/bin/env python3.2

"""
========
OVERVIEW
========

The script compares two gff3 file at base, exon and gene level and outputs the positive 
predictive value(PPV) and sensitivity(SN). The first input will be treated as a known gff3 file
and second as the predicted gff3 file. The PPV/SN in the output will be for predicted gff3 with
respect to the known gff3.


=====
INPUT
=====

Options
-------
-a1 --ANNOTATION_1       required
    Path to first annotation file(known).

-a2 --ANNOTATION_2       required
    Path to second annotation file(predicted).

-o  --output_dir          required
    Path to output directory.


======
OUTPUT
======

Script generates a summary.txt file in the output directory.
Summary.txt is a tab delimited file with following headers :

Feature    Known    Predicted    True_predicted    SN    PPV
Gene
Exon
Base

No. of predicted gene overlapping  0 known gene (new gene):  0
No. of predicted gene overlapping > 1 known gene:  115
No. of predicted gene overlaping 1 known gene :  4057
No. of known gene overlapping > 1 predicted gene :  115
No. of known gene overlapping 1 predicted gene :  4057
No. of known gene overlapping 0 predicted gene (gene missing) :  0


The summary will also be printed in the console.


====
NOTE
====

Both the input annotation file should be in the correct gff3 format.
http://www.sequenceontology.org/gff3.shtml

Examples :
chr1   genbank   gene    1    3000    .    -    .    ID=gene1
chr1   genbank   mRNA    10   2900    .    -    .    ID=mrna1;Parent=gene1
chr1   genbank   exon    10   1000    .    -    .    ID=exon1;Parent=mrna1


======
Author
======

Priti Kumari
Bioinformatics Software Engineer 
Institute for Genome Sciences
University of Maryland
Baltimore, Maryland 21201

"""


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


    a_exons = []                                    ## Set contains only uniq exons from known annotation, since multiple same exons can appear in a gff file.  
    p_exons = []                                    ## For predicted annotation

    a_gene = []
    p_gene = []

    a_mrna = []
    p_mrna = []

    exon_pred_all = set()
    gene_true = set()
    mrna_true = set()



    a_base = 0
    p_base = 0
    true_base = 0
    
    chr = []
    
    for asm_id in assemblies_1:                                                                                     ## Iterate through each chromosome from the known ref annotation        
        assembly_1 = assemblies_1[asm_id]
        assembly_2 = assemblies_2.get(asm_id,-1)                                                                    ## Find that chromosome in the predicted gff file
        genes_1 = assembly_1.genes()                                                                                ## All genes from known annotation
        anno_exons = set()

        for gene_1 in sorted(genes_1) :                                                                                     ## Add unique gene, mrna , exon features from known annotation to get each known feature total count 
            gene_1_loc = gene_1.location_on(assembly_1)
            cord = asm_id  + ":" + str(gene_1_loc.fmin) + ":" + str(gene_1_loc.fmax)+ ":"  + str(gene_1_loc.strand)        ## Use chromosome id+start+stop+strand as a string to determine uniqueness.
            if (cord not in a_gene) :
                a_gene.append(cord)
            
            for mrna_1 in sorted(gene_1.mRNAs()) :
                mrna_1_loc = mrna_1.location_on(assembly_1)
                cord = asm_id  + ":" + str(mrna_1_loc.fmin) + ":" + str(mrna_1_loc.fmax) + ":" + str(mrna_1_loc.strand)
                if (cord not in a_mrna) :
                    a_mrna.append(cord)

                for exon_1 in sorted(mrna_1.exons()) :
                    exon_1_loc = exon_1.location_on(assembly_1)
                    cord = asm_id + ":"  + str(exon_1_loc.fmin) + ":" + str(exon_1_loc.fmax) + ":" + str(exon_1_loc.strand)
                    if (cord not in a_exons) :
                        a_exons.append(cord)
                    anno_exons.add(cord)

        if (type(assembly_2) is int) :                     ##    If the chromosome is not found in prediected file, move to next chromosome.
            continue
        

        genes_2 = assembly_2.genes()                      ## All genes from predicted annotation.
        chr.append(asm_id)                                ## Append all found chromosome in a list.
        pred_exons = set()

        for gene_2 in sorted(genes_2) :                           ## Add unique gene, mrna , exon features from predicted annotation to get each predicted feature total count.  
            gene_2_loc = gene_2.location_on(assembly_2)
            cord = asm_id + ":" + str(gene_2_loc.fmin) + ":" + str(gene_2_loc.fmax) + ":" +  str(gene_2_loc.strand)
            if (cord not in p_gene) :
                p_gene.append(cord)
            
            for mrna_2 in sorted(gene_2.mRNAs()) :
                mrna_2_loc = mrna_2.location_on(assembly_2)
                cord = asm_id  + ":" + str(mrna_2_loc.fmin) + ":" + str(mrna_2_loc.fmax)+ ":" +  str(mrna_2_loc.strand)
                if (cord not in p_mrna) :
                    p_mrna.append(cord)
                
                for exon_2 in sorted(mrna_2.exons()) :
                    exon_2_loc = exon_2.location_on(assembly_2)
                    cord = asm_id  + ":" + str(exon_2_loc.fmin) + ":" + str(exon_2_loc.fmax)+ ":" + str(exon_2_loc.strand)
                    pred_exons.add(cord)
                    if (cord not in p_exons) :
                        p_exons.append(cord)
                        
        exon_pred_all.update(pred_exons.intersection(anno_exons)) # true exons
        
        
        for gene_2 in sorted(genes_2) :                                         ## From the predicted feature determine the true once. Iterate through each predicted gene sorted by cordinate
            gene_2_loc = gene_2.location_on(assembly_2)
            cord_g = asm_id  + ":"+ str(gene_2_loc.fmin) + ":" +  str(gene_2_loc.fmax) + ":" + str(gene_2_loc.strand)
            
            if (cord_g in gene_true) :                                          ## To prevent duplication, check if the feature already exists in the set of truly predicted gene.
                continue
            ex_mrna1 = set()
            ex_mrna2 = set()
			
        
            for gene_1 in sorted(genes_1) :
                gene_1_loc = gene_1.location_on(assembly_1)
                if (gene_1_loc.strand != gene_2_loc.strand) :
                    continue
                if (gene_2.overlaps_with(gene_1)) :
                    for mrna_2 in sorted(gene_2.mRNAs()) :
                        for exon_2 in sorted(mrna_2.exons()) :
                            exon_2_loc = exon_2.location_on(assembly_2)
                            cord2 = asm_id + ":" + str(exon_2_loc.fmin) + ":" + str(exon_2_loc.fmax) + ":" +  str(exon_2_loc.strand)
                            ex_mrna2.add(cord2)
                            
                    for mrna_1 in sorted(gene_1.mRNAs()) :
                        for exon_1 in sorted(mrna_1.exons()) :
                            exon_1_loc = exon_1.location_on(assembly_1)
                            cord1 = asm_id + ":" + str(exon_1_loc.fmin) + ":" + str(exon_1_loc.fmax) + ":" +  str(exon_1_loc.strand)
                            ex_mrna1.add(cord1)
                    
                    ex_union = ex_mrna1.union(ex_mrna2)
                    if (len(ex_union) ==  len(ex_mrna1) and len(ex_union) == len(ex_mrna2)) :
                    	gene_true.add(cord_g)
                    	break
          
    for asm_id in assemblies_2:                                                  ## Iterate through each chromosome from the predicted annotation
        if asm_id not in chr :
            assembly_2 = assemblies_2.get(asm_id,-1)                             ## Find that chromosome in the predicted gff file which is not found in known annotation
            genes_2 = assembly_2.genes()                                         ## Add  genes, mrna, exon features from predicted annotation to total predicted feature set.
            
            for gene_2 in sorted(genes_2) :
                gene_2_loc = gene_2.location_on(assembly_2)
                cord = asm_id + ":" + str(gene_2_loc.fmin) + ":" + str(gene_2_loc.fmax)  + ":"+ str(gene_2_loc.strand)
                if (cord not in p_gene) :
                    p_gene.append(cord)
            
                for mrna_2 in sorted(gene_2.mRNAs()) :
                    mrna_2_loc = mrna_2.location_on(assembly_2)
                    cord = asm_id  + ":" + str(mrna_2_loc.fmin) + ":" + str(mrna_2_loc.fmax) + ":" + str(mrna_2_loc.strand)
                    if (cord not in p_mrna) :
                        p_mrna.append(cord)
                    
                    for exon_2 in sorted(mrna_2.exons()) :
                        exon_2_loc = exon_2.location_on(assembly_2)
                        cord = asm_id  + ":" + str(exon_2_loc.fmin) + ":" + str(exon_2_loc.fmax) + ":" + str(exon_2_loc.strand)
                        if (cord not in p_exons) :
                            p_exons.append(cord)

    exon2_bed = args.output_dir + '/exon_2.bed'
    e_bed = open(exon2_bed, 'w')
    for exon in p_exons :
        chrom = (exon.split(':'))[0]
        start = int((exon.split(':'))[1])
        stop = int((exon.split(':'))[2])
        strand = (exon.split(':'))[3]
        if (strand == str(1)) :
            strand = "+"
        else :
            strand = "-"
        e_bed.write(chrom+"\t"+str(start)+"\t"+str(stop)+"\texon\t"+str(0)+"\t"+strand+"\n")

    e_bed.close()
    
    out2 = args.output_dir + '/exon_2_merged.bed'
    cmd = "bedtools merge -nms -scores sum -i " + exon2_bed + " -s >"+out2
    #print(cmd)
    os.system(cmd)
    
    exon1_bed = args.output_dir + '/exon_1.bed'
    e_bed = open(exon1_bed, 'w')
    for exon in a_exons :
        chrom = (exon.split(':'))[0]
        start = int((exon.split(':'))[1])
        stop = int((exon.split(':'))[2])
        strand = (exon.split(':'))[3]
        if (strand == str(1)) :
            strand = "+"
        else :
            strand = "-"
        e_bed.write(chrom+"\t"+str(start)+"\t"+str(stop)+"\texon\t"+str(0)+"\t"+strand+"\n")
    e_bed.close()

    out1 = args.output_dir + '/exon_1_merged.bed'
    cmd = "bedtools merge -nms -scores sum -i " + exon1_bed + " -s >"+out1
    #print(cmd)
    os.system(cmd)
    
    out_intersect = args.output_dir + '/exon_1_2_intersect.bed'
    cmd = "bedtools intersect -s -wo -a " + out1 + " -b " + out2 + " >" + out_intersect
    #print(cmd)
    os.system(cmd)
    
    a_base_file = open(out1,'r')
    for line in a_base_file :
        arr = line.split("\t")
        a_base = a_base + (int(arr[2]) - int(arr[1]))
    a_base_file.close()
    
    p_base_file = open(out2,'r')
    for line in p_base_file :
        arr = line.split("\t")
        p_base = p_base + (int(arr[2]) - int(arr[1]))
    p_base_file.close()

    true_base_file = open(out_intersect,'r')
    for line in true_base_file :
        arr = line.split("\t")
        true_base = true_base + int(arr[12])
    true_base_file.close()

    #Calculate SN/SP for bases 

    base_sn = (true_base/a_base) * 100                                 
    base_sp = (true_base/p_base) * 100

    #Calculate SN/SP for exons 
    annotated_exon = len(a_exons)
    predicted_exon = len(p_exons)
    true_pred_exon = len(exon_pred_all)
    
    exon_sn = (true_pred_exon/annotated_exon) * 100                                 
    exon_sp = (true_pred_exon/predicted_exon) * 100

    #Calculate SN/SP for transcript
    
    #annotated_mrna = len(a_mrna)
    #predicted_mrna = len(p_mrna)
    #true_pred_mrna = len(mrna_true)
    
    #mrna_sn = (true_pred_mrna/annotated_mrna) * 100
    #mrna_sp = (true_pred_mrna/predicted_mrna) * 100
       
    #Calculate SN/SP for genes 

    annotated_gene = len(a_gene)
    predicted_gene = len(p_gene)
    true_pred_gene = len(gene_true)

    
    temp_file7 = args.output_dir + '/true_gene.temp7.txt'
    ft7 = open(temp_file7,'w')
    for g in gene_true :
        ft7.write(g+"\n")
        
        


    
    gene_sn = (true_pred_gene/annotated_gene) * 100                                 
    gene_sp = (true_pred_gene/predicted_gene) * 100
    print("Feature\tKnown\tPredicted\tTrue_Predicted\tSN\tPPV\n")
    print("Gene\t"+str(annotated_gene)+"\t"+str(predicted_gene)+"\t"+str(true_pred_gene)+"\t"+str(gene_sn)+"\t"+str(gene_sp))
    #print("mRNA\t"+str(annotated_mrna)+"\t"+str(predicted_mrna)+"\t"+str(true_pred_mrna)+"\t"+str(mrna_sn)+"\t"+str(mrna_sp))
    print("Exon\t"+str(annotated_exon)+"\t"+str(predicted_exon)+"\t"+str(true_pred_exon)+"\t"+str(exon_sn)+"\t"+str(exon_sp))
    print("Base\t"+str(a_base)+"\t"+str(p_base)+"\t"+str(true_base)+"\t"+str(base_sn)+"\t"+str(base_sp))
    
    out_file = args.output_dir + '/summary.txt'
    if not (os.path.exists(args.output_dir)) :
        sys.exit("Directory does not exist.")
    fout = open(out_file,'w')

    fout.write("Feature\tKnown\tPredicted\tTrue_Predicted\tSN\tPPV\n")
    fout.write("Gene\t"+str(annotated_gene)+"\t"+str(predicted_gene)+"\t"+str(true_pred_gene)+"\t"+str(gene_sn)+"\t"+str(gene_sp)+"\n")
   # fout.write("mRNA\t"+str(annotated_mrna)+"\t"+str(predicted_mrna)+"\t"+str(true_pred_mrna)+"\t"+str(mrna_sn)+"\t"+str(mrna_sp)+"\n")
    fout.write("Exon\t"+str(annotated_exon)+"\t"+str(predicted_exon)+"\t"+str(true_pred_exon)+"\t"+str(exon_sn)+"\t"+str(exon_sp)+"\n")
    fout.write("Base\t"+str(a_base)+"\t"+str(p_base)+"\t"+str(true_base)+"\t"+str(base_sn)+"\t"+str(base_sp)+"\n\n")
    
    new_gene = 0
    gene_merge = 0
    gene_found = 0
    gene_split = 0
    gene_missing = 0
    altered_pred = 0
    altered_known = 0
    gene = 0



    temp_file1 = args.output_dir + '/pred_new.txt'
    temp_file2 = args.output_dir + '/pred_merged.txt'
    temp_file3 = args.output_dir + '/pred_1.txt'
    temp_file4 = args.output_dir + '/known_split.txt'
    temp_file5 = args.output_dir + '/known_1.txt'
    temp_file6 = args.output_dir + '/known_missed.txt'
    temp_file8 = args.output_dir + '/pred_altered.txt'
    temp_file9 = args.output_dir + '/known_altered.txt'
    
    ft1 = open(temp_file1,'w')
    ft2 = open(temp_file2,'w')
    ft3 = open(temp_file3,'w')
    ft4 = open(temp_file4,'w')
    ft5 = open(temp_file5,'w')
    ft6 = open(temp_file6,'w')
    ft8 = open(temp_file8,'w')
    ft9 = open(temp_file9,'w')
    
    for gene2 in p_gene :
        gene_overlap = []
        chrom2 = (gene2.split(':'))[0]
        start2 = int((gene2.split(':'))[1])
        stop2 = int((gene2.split(':'))[2])
        strand2 = (gene2.split(':'))[3]
        for gene1 in a_gene:
            chrom1 = (gene1.split(':'))[0]
            start1 = int((gene1.split(':'))[1])
            stop1 = int((gene1.split(':'))[2])
            strand1 = (gene1.split(':'))[3]
            if (chrom1 != chrom2) :
                continue
            if (strand1 != strand2) :
                continue
            if (start1 > stop2) :
                break
            if(start1 <= stop2 and start2 <= stop1) :
                arr = [start1,stop1,start2,stop2]
                arr.sort()
                len_overlap = arr[2] - arr[1]
                per_overlap = (len_overlap/(stop1 - start1)) * 100
                gene_overlap.append(per_overlap)

        if (len(gene_overlap) == 0) :
            new_gene += 1
            ft1.write(gene2+"\n")
            
        if (len(gene_overlap) > 1) :
            true_overlap = 0
            for overlap in gene_overlap :
                if(overlap >= 50) :
                    true_overlap += 1;
            if (true_overlap >= 2) :
                gene_merge += 1
                ft2.write(gene2+"\n")
            else :
                altered_pred += 1;
                ft8.write(gene2+"\n")
                
        if (len(gene_overlap) == 1) :
            gene_found += 1
            ft3.write(gene2+"\n")
            
        
    for gene1 in a_gene :
        gene_overlap = []
        chrom1 = (gene1.split(':'))[0]
        start1 = int((gene1.split(':'))[1])
        stop1 = int((gene1.split(':'))[2])
        strand1 = (gene1.split(':'))[3]
        for gene2 in p_gene:
            chrom2 = (gene2.split(':'))[0]
            start2 = int((gene2.split(':'))[1])
            stop2 = int((gene2.split(':'))[2])
            strand2 = (gene2.split(':'))[3]
            if (chrom1 != chrom2) :
                continue
            if (strand1 != strand2) :
                continue
            if (start2 > stop1) :
                break
            if(start1 <= stop2 and start2 <= stop1) :
                arr = [start1,stop1,start2,stop2]
                arr.sort()
                len_overlap = arr[2] - arr[1]
                per_overlap = (len_overlap/(stop2 - start2)) * 100
                gene_overlap.append(per_overlap)

                
        if (len(gene_overlap) > 1) :
            true_overlap = 0
            for overlap in gene_overlap :
                if(overlap >= 50) :
                    true_overlap += 1;
            if (true_overlap >= 2) :
                gene_split += 1
                ft4.write(gene1+"\n")
            else :
                altered_known += 1
                ft9.write(gene1+"\n")
        
        if (len(gene_overlap) == 1) :
            gene += 1
            ft5.write(gene1+"\n")
        if (len(gene_overlap) == 0) :
            gene_missing += 1
            ft6.write(gene1+"\n")
            
    print ("1. No. of predicted gene overlapping  0 known gene (new gene): ",new_gene)
    print ("2. No. of predicted gene overlapping > 1 known gene by at least 50%: ",gene_merge)
    print ("3. No. of altered predicted gene: ",altered_pred)
    print ("4. No. of predicted gene overlaping 1 known gene : ",gene_found)
    print ("5. No. of known gene overlapping > 1 predicted gene by at least 50% : ",gene_split)
    print ("6. No. of altered known gene: ",altered_known)
    print ("7. No. of known gene overlapping 1 predicted gene : ",gene)
    print ("8. No. of known gene overlapping 0 predicted gene (gene missing) : ",gene_missing)


    fout.write ("1. No. of predicted gene overlapping  0 known gene (new gene): "+str(new_gene)+"\n")
    fout.write ("2. No. of predicted gene overlapping > 1 known gene by at least 50%: "+str(gene_merge)+"\n")
    fout.write ("3. No. of altered predicted gene: "+str(altered_pred)+"\n")
    fout.write ("4. No. of predicted gene overlaping 1 known gene : "+str(gene_found)+"\n")
    fout.write ("5. No. of known gene overlapping > 1 predicted gene by at least 50% : "+str(gene_split)+"\n")
    fout.write ("6. No. of altered known gene: "+str(altered_known)+"\n")
    fout.write ("7. No. of known gene overlapping 1 predicted gene : "+str(gene)+"\n")
    fout.write ("8. No. of known gene overlapping 0 predicted gene (gene missing) : "+str(gene_missing)+"\n")

    
    fout.close()
    ft1.close()
    ft2.close()
    ft3.close()
    ft4.close()
    ft5.close()
    ft6.close()
    ft7.close()
    ft8.close()
    ft9.close()


    #Clean up
    cmd = "rm " + args.output_dir + "/*.bed"
    os.system(cmd)



    
if __name__ == '__main__':
    args = interface()
    process_files(args)






