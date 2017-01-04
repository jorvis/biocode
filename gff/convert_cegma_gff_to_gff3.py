#!/usr/bin/env python3

"""
The default CEGMA output, while technically GFF, has only a subsection of the conventional
features and lacks relationships, so any displays reading the GFF show very segmented matches.  This
script corrects this.

Note: This was written/tested using the file created by cegma named like 'output.cegma.gff', not the
variant with 'local' in the file name.


Example input (CEGMA 2.4):
    jcf7180000788954        cegma   First   113875  113877  -6.54   +       0       KOG0002.9
    jcf7180000788954        cegma   Exon    113875  113877  -6.54   +       .       KOG0002.9
    jcf7180000788954        cegma   Internal        113937  114040  36.65   +       0       KOG0002.9
    jcf7180000788954        cegma   Exon    113937  114040  36.65   +       .       KOG0002.9
    jcf7180000788954        cegma   Terminal        114093  114141  19.71   +       1       KOG0002.9
    jcf7180000788954        cegma   Exon    114093  114141  19.71   +       .       KOG0002.9

    or

    jcf7180000788614        cegma   Single  100757  101674  133.99  +       0       KOG0003.5
    jcf7180000788614        cegma   Exon    100757  101674  133.99  +       .       KOG0003.5

Example input (CEGMA 2.5):
    Scaffold9	cegma	First	1282294	1282324	2.54	-	0	KOG0018.9
    Scaffold9	cegma	Internal	1282135	1282228	42.86	-	2	KOG0018.9
    Scaffold9	cegma	Internal	1280282	1282085	512.31	-	1	KOG0018.9
    Scaffold9	cegma	Internal	1279142	1280233	308.01	-	0	KOG0018.9
    Scaffold9	cegma	Terminal	1278352	1279071	304.09	-	0	KOG0018.9

    or

    Scaffold9	cegma	First	870770	871027	142.35	-	0	KOG0019.4
    Scaffold9	cegma	Terminal	868869	870713	1029.84	-	0	KOG0019.4

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

from biocode import things


def main():
    parser = argparse.ArgumentParser( description='Converts CEGMA GFF output to spec-legal GFF3')

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
    
    next_id_nums = {'gene':1, 'mRNA':1, 'CDS':1, 'exon':1}
    exon_column_types = ['First', 'Internal', 'Terminal', 'Single']

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
                assemblies[mol_id] = things.Assembly(id=mol_id)

            current_assembly = assemblies[mol_id]
                
            feat_id = "cegma.gene.{0}".format(next_id_nums['gene'])
            next_id_nums['gene'] += 1
            gene = things.Gene(id=feat_id)
            current_gene = gene
            current_gene_strand = cols[6]
            current_gene_fmin = feat_fmin
            current_gene_fmax = feat_fmax

            mRNA_id = "cegma.mRNA.{0}".format(next_id_nums['mRNA'])
            next_id_nums['mRNA'] += 1
            
            mRNA = things.mRNA(id=mRNA_id, parent=gene)
            gene.add_mRNA(mRNA)
            current_mRNA = mRNA

        # CEGMA versions < 2.5 had two rows for each exon.  We don't need to process both of them, so
        #  we skip the Exon one because its phase information is incorrect.
        if feat_type in exon_column_types:
            CDS_id = "cegma.CDS.{0}".format(next_id_nums['CDS'])
            next_id_nums['CDS'] += 1
            CDS = things.CDS(id=CDS_id, parent=current_mRNA)
            CDS.locate_on( target=current_assembly, fmin=feat_fmin, fmax=feat_fmax, strand=cols[6], phase=cols[7] )
            current_mRNA.add_CDS(CDS)

            exon_id = "cegma.exon.{0}".format(next_id_nums['exon'])
            next_id_nums['exon'] += 1
            exon = things.Exon(id=exon_id, parent=current_mRNA)
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





















