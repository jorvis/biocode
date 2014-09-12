#!/usr/bin/env python3

import argparse
import os
import re
import biothings
import biocodegff

'''
This script converts the GFF-ish output of Metagenemark into legal GFF3 with full,
canonical gene models.

http://www.sequenceontology.org/gff3.shtml

EXPECTED INPUT example (the commented lines are ignored):

855 length:1510	GeneMark.hmm	gene	1	852	.	-	0	gene_id 1
##Protein 1
##MAAKYDDVQELLRKKADINARLSLLAYDGTPEIKNRGNGKYLYTRKRVSGKLTSTYVGVY
##SDDLYNLLLRNASESRQLRKELRAVVRQLANAGYSNSELSADLVNNIAFARANLKANIYD
##QAVLEGVATTFPQTEEIIDNGKVFGVSASDVQKILNLKHAWEFILDEDVVLSRSDYYMLS
##HIAKIVNEGFFLDGGRIRGVPVAIGGSTYIPPIPNEIDVKEKIRNIVEDEGDSENVEGRV
##RRKEPIDIAIELCVYCMKSQIFLDGNKRASVIFANHYLISHGIG
##end-Protein

35 length:9081	GeneMark.hmm	gene	1	3378	.	-	0	gene_id 2
35 length:9081	GeneMark.hmm	gene	3974	4357	.	-	0	gene_id 3
35 length:9081	GeneMark.hmm	gene	4398	4652	.	-	0	gene_id 4
35 length:9081	GeneMark.hmm	gene	4792	9081	.	-	0	gene_id 5
##Protein 2
##MNDKMQMHRLIEQKRIEGADKKPRYGMRKLTIGTVSCLLGFASLLAFTTPNFSQAAESGN
##TGGGTNSLITGTVPKENKSSSEEGEKLRAPQDVSGQLENLKIKLSGDNHENASLIHPLRP
##ADDNDESSDQQMKVNFSFDVDGSEIKEGDYFDLNLSNNLNLYGATSKKADIQTKLYVGND
##LVAKGVYDVENHRMRYTFTKKASEYGHFSQSVMEPVFIDAKGVPTNNKNVAVEASIGNHK
##SQKEVEVTYDLQPQAGQDNLNSNGTANLFDIDESTGTYKETIYVNNKQREQNNTRILIEN


EXAMPLE OUTPUT:

855	GeneMark.hmm	gene	1	852	.	-	.	ID=HUZ239124.gene.1
855	GeneMark.hmm	mRNA	1	852	.	-	.	ID=HUZ239124.mRNA.1;Parent=HUZ239124.gene.1
855	GeneMark.hmm	CDS	1	852	.	-	0	ID=HUZ239124.CDS.1;Parent=HUZ239124.mRNA.1
855	GeneMark.hmm	exon	1	852	.	-	.	ID=HUZ239124.exon.1;Parent=HUZ239124.mRNA.1

'''

def main():
    parser = argparse.ArgumentParser( description='Metagenemark GFF -> GFF3 conversion script')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to a GFF file created by Metagenemark' )
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-p', '--prefix', type=str, required=True, help='Prefix to use in ID generation')
    args = parser.parse_args()

    assemblies = dict()
    current_assembly = None
    
    gene = None
    mRNAs = dict()
    current_sequence = None
    current_gene_comment_lines = list()

    fout = open(args.output, mode='wt', encoding='utf-8')
    fout.write("##gff-version 3\n")

    for line in open(args.input):
        if line.startswith("#"):
            current_gene_comment_lines.append(line)

        else:
            cols = line.split("\t")

            if len(cols) != 9:
                continue

            mol_id = cols[0]
            mol_id_m = re.match('^(\S+) ', mol_id)

            if mol_id_m:
                print("MATCH!")
                mol_id = mol_id_m.group(1)
            
            feat_type = cols[2]

            ## we expect only gene types here
            if feat_type != 'gene':
                raise Exception("ERROR: expected only 'gene' feature types as input.")

            m_gene = re.match('gene_id (\d+)', cols[8])

            if m_gene:
                gene_num = m_gene.group(1)
            else:
                raise Exception("ERROR: expected 9th column to have gene ids like: gene_id 5")

            ## initialize this assembly if we haven't seen it yet
            if mol_id not in assemblies:
                assemblies[mol_id] = biothings.Assembly( id=mol_id )

            current_assembly = assemblies[mol_id]

            gene = biothings.Gene( id="{0}.gene.{1}".format(args.prefix, gene_num) )
            gene.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6] )

            mRNA = biothings.mRNA( id="{0}.mRNA.{1}".format(args.prefix, gene_num), parent=gene.id )
            mRNA.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6] )
            gene.add_mRNA(mRNA)

            CDS = biothings.CDS( id="{0}.CDS.{1}".format(args.prefix, gene_num), parent=mRNA.id )
            CDS.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6], phase=int(cols[7]) )
            mRNA.add_CDS(CDS)

            exon = biothings.Exon( id="{0}.exon.{1}".format(args.prefix, gene_num), parent=mRNA.id )
            exon.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6] )
            mRNA.add_exon(exon)

            gene.print_as(fh=fout, source='GeneMark.hmm', format='gff3')
            fout.write( "".join(current_gene_comment_lines) )
            current_gene_comment_lines = list()
                    

if __name__ == '__main__':
    main()







