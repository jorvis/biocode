#!/usr/bin/env python3

import argparse
import re

from biocode import things

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


EXPECTED INPUT (v3.25+)

In my tests of v 3.25 the use of 'gene' in the column above has been replaced with CDS instead.
The rest appears to be the same.

SRS019986_Baylor_scaffold_53    GeneMark.hmm    CDS     3       512     -635.997977     +       0       gene_id=18, length=510, gene_score=-635.997977, rbs_score=-0.020000, rbs_spacer=-1, stop_enforced=N, start_codon=0, logodd=37.242660
SRS019986_Baylor_scaffold_53    GeneMark.hmm    CDS     530     1501    -1202.233893    +       0       gene_id=19, length=972, gene_score=-1202.233893, rbs_score=-0.020000, rbs_spacer=-1, stop_enforced=N, start_codon=0, logodd=62.393607
SRS019986_Baylor_scaffold_53    GeneMark.hmm    CDS     1603    2109    -608.550058     +       0       gene_id=20, length=507, gene_score=-608.550058, rbs_score=-0.020000, rbs_spacer=-1, stop_enforced=N, start_codon=0, logodd=52.060246


EXAMPLE OUTPUT:

855	GeneMark.hmm	gene	1	852	.	-	.	ID=HUZ239124.gene.1
855	GeneMark.hmm	mRNA	1	852	.	-	.	ID=HUZ239124.mRNA.1;Parent=HUZ239124.gene.1
855	GeneMark.hmm	CDS	1	852	.	-	0	ID=HUZ239124.CDS.1;Parent=HUZ239124.mRNA.1
855	GeneMark.hmm	exon	1	852	.	-	.	ID=HUZ239124.exon.1;Parent=HUZ239124.mRNA.1
855	GeneMark.hmm	polypeptide	1	852	.	-	.	ID=HUZ239124.polypeptide.1;Parent=HUZ239124.mRNA.1

'''

def main():
    parser = argparse.ArgumentParser( description='Metagenemark GFF -> GFF3 conversion script')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to a GFF file created by Metagenemark' )
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-p', '--prefix', type=str, required=True, help='Prefix to use in ID generation')
    parser.add_argument('-pf', '--protein_fasta', type=str, required=False, help='Optional protein FASTA to be written')
    args = parser.parse_args()

    assemblies = dict()
    current_assembly = None

    # key like 2 = SRS014890.polypeptide.2
    polypeptide_lookup = dict()
    writing_protein = False
    
    gene = None
    mRNAs = dict()
    current_sequence = None
    current_gene_comment_lines = list()

    fout = open(args.output, mode='wt', encoding='utf-8')
    fout.write("##gff-version 3\n")

    if args.protein_fasta is not None:
        protein_out = open(args.protein_fasta, mode='wt', encoding='utf-8')

    for line in open(args.input):
        if line.startswith("#"):
            if line.startswith("##FASTA"):
                current_gene_comment_lines.append("#{0}".format(line))
                
            elif line.startswith("##end-Protein"):
                writing_protein = False
                current_gene_comment_lines.append(line)
                
            # since we're already doing our own header, don't duplicate the old one
            elif line.startswith("##gff-version"):
                continue
            else:
                if line.startswith("##Protein "):
                    m = re.match("##Protein (\d+)", line)
                    if m:
                        writing_protein = True
                        protein_out.write(">{0}\n".format(polypeptide_lookup[m.group(1)]))
                    else:
                        raise Exception("ERROR: Expected line to match: ##Protein N")
                elif writing_protein == True:
                    protein_out.write(line[2:])
                    
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
            if feat_type not in ['gene', 'CDS']:
                raise Exception("ERROR: expected only 'gene' or 'CDS' feature types as input (depending on metagenemark version).")

            m_gene = re.match('gene_id[ =](\d+)', cols[8])

            if m_gene:
                gene_num = m_gene.group(1)
            else:
                raise Exception("ERROR: expected 9th column to have gene ids like: gene_id 5")

            ## initialize this assembly if we haven't seen it yet
            if mol_id not in assemblies:
                assemblies[mol_id] = things.Assembly(id=mol_id)

            current_assembly = assemblies[mol_id]

            gene = things.Gene(id="{0}.gene.{1}".format(args.prefix, gene_num))
            gene.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6] )

            mRNA = things.mRNA(id="{0}.mRNA.{1}".format(args.prefix, gene_num), parent=gene.id)
            mRNA.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6] )
            gene.add_mRNA(mRNA)

            CDS = things.CDS(id="{0}.CDS.{1}".format(args.prefix, gene_num), parent=mRNA.id)
            CDS.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6], phase=int(cols[7]) )
            mRNA.add_CDS(CDS)

            exon = things.Exon(id="{0}.exon.{1}".format(args.prefix, gene_num), parent=mRNA.id)
            exon.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6] )
            mRNA.add_exon(exon)

            polypeptide_id = "{0}.polypeptide.{1}".format(args.prefix, gene_num)
            polypeptide = things.Polypeptide(id=polypeptide_id, parent=mRNA.id)
            polypeptide.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6] )
            mRNA.add_polypeptide(polypeptide)
            polypeptide_lookup[gene_num] = polypeptide_id

            gene.print_as(fh=fout, source='GeneMark.hmm', format='gff3')
            fout.write( "".join(current_gene_comment_lines) )
            current_gene_comment_lines = list()
                    

if __name__ == '__main__':
    main()







