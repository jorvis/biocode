#!/usr/bin/env python3.2

import argparse
import os
import re
import biothings
import biocodegff

'''
Admittedly, when you run augustus with --gff3=on you will get legal GFF3 with lots of
comments.  What I often need is GFF3 that also follows the gene model graph conventions
defined in the GFF3, which the default output does not.  There are no exon or mRNA
features, for example.

This script reads the native GFF output and makes it conform to the suggested canonical
gene model, defined here:

http://www.sequenceontology.org/gff3.shtml

WARNING: The augustus output format appears to change pretty frequently.  I've provided
an example of the expected input below, and the version I tested with is 2.7

EXPECTED INPUT example:

###
# start gene g2
FO082871        AUGUSTUS        gene    16991   18816   0.04    +       .       ID=g2
FO082871        AUGUSTUS        transcript      16991   18816   0.04    +       .       ID=g2.t1;Parent=g2
FO082871        AUGUSTUS        transcription_start_site        16991   16991   .       +       .       Parent=g2.t1
FO082871        AUGUSTUS        five_prime_utr  16991   17193   0.04    +       .       Parent=g2.t1
FO082871        AUGUSTUS        start_codon     17194   17196   .       +       0       Parent=g2.t1
FO082871        AUGUSTUS        intron  17317   17557   0.15    +       .       Parent=g2.t1
FO082871        AUGUSTUS        CDS     17194   17316   0.15    +       0       ID=g2.t1.cds;Parent=g2.t1
FO082871        AUGUSTUS        CDS     17558   17647   0.15    +       0       ID=g2.t1.cds;Parent=g2.t1
FO082871        AUGUSTUS        stop_codon      17645   17647   .       +       0       Parent=g2.t1
FO082871        AUGUSTUS        three_prime_utr 17648   18816   0.07    +       .       Parent=g2.t1
FO082871        AUGUSTUS        transcription_end_site  18816   18816   .       +       .       Parent=g2.t1
# protein sequence = [MQDNYTDGQPGIHRAIAKLEAIINEVDQQLYDYLLERQIHLSPPSKYWEISDIDELVAKAYVCFYAGLCC]
# Evidence for and against this transcript:
# % of transcript supported by hints (any source): 80
# CDS exons: 2/2
#      W:   2 
# CDS introns: 0/1
# 5'UTR exons and introns: 1/1
#      W:   1 
# 3'UTR exons and introns: 1/1
#      W:   1 
# hint groups fully obeyed: 148
#      W: 148 
# incompatible hint groups: 11
#      W:  11 
# end gene g2

'''

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to a GFF file created by Augustus' )
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    assemblies = dict()
    current_assembly = None
    
    gene = None
    mRNAs = dict()
    in_sequence = False
    current_sequence = None
    current_gene_comment_lines = list()

    fout = open(args.output, mode='wt', encoding='utf-8')
    fout.write("##gff-version 3")

    for line in open(args.input):
        if line.startswith("#"):
            current_gene_comment_lines.append(line)
            
            if line.startswith("# end gene "):
                ## purge the comments, then write the gene
                fout.write( "".join(current_gene_comment_lines) )
                gene.print_as(fh=fout, source='AUGUSTUS', format='gff3')

                gene = None
                mRNAs = dict()
                in_sequence = False
                current_sequence = None
                current_gene_comment_lines = list()

            elif line.startswith("# protein sequence = ["):
                pass
            elif in_sequence is True:
                # build 'current_sequence'
                pass

        else:
            cols = line.split("\t")

            if len(cols) != 9:
                continue

            mol_id = cols[0]
            feat_type = cols[2]
            feat_id = biocodegff.column_9_value(cols[8], 'ID')

            ## initialize this assembly if we haven't seen it yet
            if mol_id not in assemblies:
                assemblies[mol_id] = biothings.Assembly( id=mol_id )

            current_assembly = assemblies[mol_id]

            if feat_type == "gene":
                gene = biothings.Gene( id=feat_id )
                gene.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6] )
            elif feat_type == "transcript":
                mRNA = biothings.mRNA( id=feat_id, parent=gene )
                mRNA.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6] )
                gene.add_mRNA(mRNA)
                mRNAs[mRNA.id] = mRNA
            elif feat_type == "CDS":
                parent_id = biocodegff.column_9_value( cols[8], 'Parent' )

                ## sanity check that we've seen this parent
                if parent_id not in mRNAs:
                    raise Exception("ERROR: Found CDS column with parent ({0}) mRNA not yet in the file".format(parent_id))

                CDS = biothings.CDS( id=feat_id, parent=mRNAs[parent_id] )
                CDS.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6], phase=int(cols[7]) )
                mRNA.add_CDS(CDS)

                exon = biothings.Exon( id=feat_id, parent=mRNAs[parent_id] )
                exon.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6] )
                mRNA.add_exon(exon)
        

if __name__ == '__main__':
    main()







