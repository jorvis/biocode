#!/usr/bin/env python3.2

import argparse
import os
import biothings
import biocodegff
import sys
import re

## for debugging
from pprint import pprint

## tracked in SILLAB-unk, this script was written to answer a specific set of questions, namely:

# 1.1 How many genes in the current annotation have evidence from the assembled transcripts?
# What I would say counts for evidence is the whole gene included within the start/end of the transcript coordinates (regardless of whether or not the transcript covers multiple genes); I say for now we ignore partial overlaps.

# 1.2 How many genes structures will be modified as a result of the transcriptome sequencing? 
# This is obviously a very broad and involved question (or the answer is, rather, because of all this issue of multigene transcripts). So, we need a simplified version of this question that might give us the lower boundary of that number, while being reasonable doable. So, I was thinking two things:

# 1.2.1 How many transcripts cover no currently annotated genes? I have seen a few in IGV while randomly browsing, so I assume this is a non-trivial number…

# 1.2.2 Would it be possible to determine, of the genes with transcript coverage, how many introns have incorrect boundaries? Essentially, the %incorrect introns, and how many genes that maps to.

# Here, I assume 1.2.1 is easier than 1.2.2, but I don't know if you already have the necessary tools to actually do 1.2.2. If you do, fine, if not, we could see what number 1.2.1 comes up to, and see if that shows enough of an impact, before determining if we need to go on to 1.2.2.


def main():
    parser = argparse.ArgumentParser( description='Reports statistics of reference gene coverage and extension by aligned RNA-seq transcript data.')

    ## output file to be written
    parser.add_argument('-r', '--reference_file', type=str, required=True, help='GFF3 file of a reference annotation' )
    parser.add_argument('-q', '--alignment_file', type=str, required=True, help='GFF3 file with RNA-seq assembly transcript features aligned to the same reference genome.  Usually with something like GMAP.' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    (ref_assemblies, ref_features) = parse_gff3( args.reference_file )
    (qry_assemblies, qry_features) = parse_gff3( args.alignment_file )

    ## for stats
    total_ref_gene_count = 0
    total_qry_gene_count = 0
    qry_genes_not_overlapping_a_ref_gene = list()
    nonoverlapping_qry_genes_not_overlapping_a_ref_gene = list()
    

    print("STATS: Found {0} reference assemblies with {1} features".format(len(ref_assemblies), len(ref_features) ) )
    print("STATS: Found {0} query assemblies with {1} features".format(len(qry_assemblies), len(qry_features) ) )

    for asm_id in ref_assemblies:
        ref_assembly = ref_assemblies[asm_id]
        qry_assembly = None

        ## only continue processing if we have both qry and ref assemblies with the same name
        if asm_id in qry_assemblies:
            qry_assembly = qry_assemblies[asm_id]
        else:
            continue

        print("DEBUG: comparing assembly {0}".format(asm_id) )
        ref_gene_count = len(ref_assembly.genes())
        qry_gene_count = len(qry_assembly.genes())
        total_ref_gene_count += ref_gene_count
        total_qry_gene_count += qry_gene_count
        print("STATS: assembly:{2} - ref genes: ({0}), qry genes: ({1})".format( ref_gene_count, qry_gene_count, asm_id ) )
        
        for ref_gene in ref_assembly.genes():
            print("REF_GENE {0}".format(ref_gene.id) )
            for qry_gene in qry_assembly.genes():
                overlap = ref_gene.overlaps_with( qry_gene )

                if overlap is True:
                    print("DEBUG: {0} and {1} appear to overlap".format(ref_gene.id, qry_gene.id) )

        ## now the opposite comparison (for 1.2.1 above)
        for qry_gene in qry_assembly.genes():
            found_overlap = False
            
            for ref_gene in ref_assembly.genes():
                if ref_gene.overlaps_with( qry_gene ):
                    found_overlap = True

            if found_overlap == False:
                qry_genes_not_overlapping_a_ref_gene.append(qry_gene)
        
    ## Now, for those qry genes which didn't map to a reference compare within
    #   the same set to see how many of these are not overlapping themselves
    nonoverlapping_qry_genes_not_overlapping_a_ref_gene = reduce_overlaps( qry_genes_not_overlapping_a_ref_gene )

    print("STATS: total reference genes found: {0}".format(total_ref_gene_count) )
    print("STATS: total query genes found: {0}".format(total_qry_gene_count) )
    print("STATS: total query genes not overlapping a reference gene: {0}".format(len(qry_genes_not_overlapping_a_ref_gene)) )
    print("STATS:   of these, count which are non-overlapping: {0}".format(len(nonoverlapping_qry_genes_not_overlapping_a_ref_gene)) )


def reduce_overlaps( things ):
    handled_ids = dict()
    nonoverlapping_set = list()
    
    for qry_gene in things:
        if qry_gene.id in handled_ids:
            continue
        
        ## mark this one as handled
        handled_ids[qry_gene.id] = 1
        nonoverlapping_set.append(qry_gene)

        ## then mark any that align to it (except self)
        for sbj_gene in things:
            if qry_gene.id == sbj_gene.id:
                continue
            elif qry_gene.overlaps_with(sbj_gene):
                handled_ids[sbj_gene.id] = 1
            
    return nonoverlapping_set

def parse_gff3(gff3_file):

    assemblies = dict()
    features   = dict()
    
    for line in open(gff3_file):
        cols = line.split("\t")

        if len(cols) != 9:
            continue

        mol_id = cols[0]
        
        ## initialize this assembly if we haven't seen it yet
        if mol_id not in assemblies:
            assemblies[mol_id] = biothings.Assembly( id=mol_id )

        current_assembly = assemblies[mol_id]
        rfmin = int(cols[3]) - 1
        rfmax = int(cols[4])
        rstrand = None
        feat_id = biocodegff.column_9_value(cols[8], 'ID')
        parent_id = biocodegff.column_9_value(cols[8], 'Parent')
        parent_feat = None
        
        if parent_id is not None:
            if parent_id in features:
                parent_feat = features[parent_id]
            else:
                raise Exception("Error in GFF3: Parent {0} referenced by a child feature before it was defined".format(parent_id) )

        if cols[6] == '-':
            rstrand = -1
        elif cols[6] == '+':
            rstrand = 1
        else:
            rstrand = 0
            
        if cols[2] == 'gene':
            gene = biothings.Gene(id=feat_id)
            gene.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            features[feat_id] = gene
            current_assembly.add_gene(gene)

        
        elif cols[2] == 'mRNA':
            mRNA = biothings.mRNA(id=feat_id, parent=parent_feat)
            mRNA.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            parent_feat.add_mRNA(mRNA)
            features[feat_id] = mRNA

        elif cols[2] == 'rRNA':
            rRNA = biothings.rRNA(id=feat_id, parent=parent_feat)
            rRNA.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            parent_feat.add_rRNA(rRNA)
            features[feat_id] = rRNA
            
        elif cols[2] == 'tRNA':
            tRNA = biothings.tRNA(id=feat_id, parent=parent_feat)
            tRNA.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            parent_feat.add_tRNA(tRNA)
            features[feat_id] = tRNA

        elif cols[2] == 'exon':
            exon = biothings.Exon(id=feat_id, parent=parent_feat)
            exon.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            parent_feat.add_exon(exon)
            features[feat_id] = exon

        elif cols[2] == 'CDS':
            CDS = biothings.CDS(id=feat_id, parent=parent_feat)
            CDS.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            parent_feat.add_CDS(CDS)
            features[feat_id] = CDS

    return (assemblies, features)



if __name__ == '__main__':
    main()




