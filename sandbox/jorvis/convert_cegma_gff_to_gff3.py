#!/usr/bin/env python3

"""
The default CEGMA output, while technically GFF, has only a subsection of the conventional
features and lacks relationships, so any displays reading the GFF show very segmented matches.  This
script corrects this.

Note: This was written/tested using the file created by cegma named like 'output.cegma.gff', not the
variant with 'local' in the file name.


Example input:
    jcf7180000788954        cegma   First   113875  113877  -6.54   +       0       KOG0002.9
    jcf7180000788954        cegma   Exon    113875  113877  -6.54   +       .       KOG0002.9
    jcf7180000788954        cegma   Internal        113937  114040  36.65   +       0       KOG0002.9
    jcf7180000788954        cegma   Exon    113937  114040  36.65   +       .       KOG0002.9
    jcf7180000788954        cegma   Terminal        114093  114141  19.71   +       1       KOG0002.9
    jcf7180000788954        cegma   Exon    114093  114141  19.71   +       .       KOG0002.9

    or

    jcf7180000788614        cegma   Single  100757  101674  133.99  +       0       KOG0003.5
    jcf7180000788614        cegma   Exon    100757  101674  133.99  +       .       KOG0003.5

Example output:

    jcf7180000788954        cegma   gene    113875  114141  .       +       .       ID=cegma.gene.1
    jcf7180000788954        cegma   mRNA    113875  114141  .       +       .       ID=cegma.mRNA.1;Parent=cegma.gene.1
    jcf7180000788954        cegma   CDS     113875  113877  .       +       0       ID=cegma.CDS.1;Parent=cegma.mRNA.1
    jcf7180000788954        cegma   CDS     113937  114040  .       +       0       ID=cegma.CDS.2;Parent=cegma.mRNA.1
    jcf7180000788954        cegma   CDS     114093  114141  .       +       1       ID=cegma.CDS.3;Parent=cegma.mRNA.1
    jcf7180000788954        cegma   exon    113875  113877  .       +       .       ID=cegma.exon.1;Parent=cegma.mRNA.1
    jcf7180000788954        cegma   exon    113937  114040  .       +       .       ID=cegma.exon.2;Parent=cegma.mRNA.1
    jcf7180000788954        cegma   exon    114093  114141  .       +       .       ID=cegma.exon.3;Parent=cegma.mRNA.1

    or 

    jcf7180000788614        cegma   gene    100757  101674  .       +       .       ID=cegma.gene.2
    jcf7180000788614        cegma   mRNA    100757  101674  .       +       .       ID=cegma.mRNA.2;Parent=cegma.gene.2
    jcf7180000788614        cegma   CDS     100757  101674  .       +       0       ID=cegma.CDS.4;Parent=cegma.mRNA.2
    jcf7180000788614        cegma   exon    100757  101674  .       +       .       ID=cegma.exon.4;Parent=cegma.mRNA.2


Author: Joshua Orvis
Contact: jorvis AT gmail
"""

import argparse
import os
import biothings


def main():
    parser = argparse.ArgumentParser( description='Converts CEGMA GFF output to GFF3')

    # output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to parse' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    
    args = parser.parse_args()

    fout = open(args.output_file, 'w')
    fout.write("##gff-version 3\n")

    assemblies = dict()
    current_assembly = None
    current_gene = None
    current_mRNA = None

    current_gene_fmin = None
    current_gene_fmax = None
    current_gene_strand = None
    
    last_phase = None

    next_id_nums = {'gene':1, 'mRNA':1, 'CDS':1, 'exon':1}

    for line in open(args.input_file, 'r'):
        if line.startswith('#'):
            continue
        
        cols = line.split("\t")

        if len(cols) != 9:
            continue

        mol_id = cols[0]
        feat_type = cols[2]
        feat_fmin = int(cols[3]) - 1
        feat_fmax = int(cols[4])
    
        if feat_type == 'Single' or feat_type == 'First':
            # If there's an existing gene already, print it out
            if current_gene is not None:
                current_gene.locate_on( target=current_assembly, fmin=current_gene_fmin, fmax=current_gene_fmax, strand=current_gene_strand )
                current_mRNA.locate_on( target=current_assembly, fmin=current_gene_fmin, fmax=current_gene_fmax, strand=current_gene_strand )

                #current_gene.print_as(format='text')
                current_gene.print_as(fh=fout, source='cegma', format='gff3')

            # initialize this assembly if we haven't seen it yet
            if mol_id not in assemblies:
                assemblies[mol_id] = biothings.Assembly( id=mol_id )
                last_assembly = current_assembly
                current_assembly = assemblies[mol_id]
                
            feat_id = "cegma.gene.{0}".format(next_id_nums['gene'])
            next_id_nums['gene'] += 1
            gene = biothings.Gene( id=feat_id )
            current_gene = gene
            current_gene_strand = cols[6]
            current_gene_fmin = feat_fmin
            current_gene_fmax = feat_fmax

            mRNA_id = "cegma.mRNA.{0}".format(next_id_nums['mRNA'])
            next_id_nums['mRNA'] += 1
            
            mRNA = biothings.mRNA( id=mRNA_id, parent=gene )
            gene.add_mRNA(mRNA)
            current_mRNA = mRNA

            last_phase = int(cols[7])
            
        # CEGMA is pretty annoying here in that it has two rows for all exonic features.  In the first
        #  its type classifies them and the second is actually the 'Exon' with the same coordintes.
        #  The biggest pain of it is that the phase is correct on the classifying feature, but omitted
        #  from the Exon one.  So we track that and use the last one.
        elif feat_type == 'Internal' or feat_type == 'Terminal':
            last_phase = int(cols[7])
            
        elif feat_type == 'Exon':
            CDS_id = "cegma.CDS.{0}".format(next_id_nums['CDS'])
            next_id_nums['CDS'] += 1
            CDS = biothings.CDS( id=CDS_id, parent=current_mRNA )
            CDS.locate_on( target=current_assembly, fmin=feat_fmin, fmax=feat_fmax, strand=cols[6], phase=last_phase )
            current_mRNA.add_CDS(CDS)

            #print("DEBUG: Added CDS {0}".format(CDS))

            exon_id = "cegma.exon.{0}".format(next_id_nums['exon'])
            next_id_nums['exon'] += 1
            exon = biothings.Exon( id=exon_id, parent=current_mRNA )
            exon.locate_on( target=current_assembly, fmin=feat_fmin, fmax=feat_fmax, strand=cols[6] )
            mRNA.add_exon(exon)

            if feat_fmin < current_gene_fmin:
                current_gene_fmin = feat_fmin

            if feat_fmax > current_gene_fmax:
                current_gene_fmax = feat_fmax
            

    # don't forget the last gene
    if current_gene is not None:
        current_gene.locate_on( target=current_assembly, fmin=current_gene_fmin, fmax=current_gene_fmax, strand=current_gene_strand )
        current_mRNA.locate_on( target=current_assembly, fmin=current_gene_fmin, fmax=current_gene_fmax, strand=current_gene_strand )
        current_gene.print_as(fh=fout, source='cegma', format='gff3')
            
            


if __name__ == '__main__':
    main()





















