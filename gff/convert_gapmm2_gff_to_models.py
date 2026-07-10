#!/usr/bin/env python3

'''
The gapmm2 tool creates GFF3 output which is valid, but not in the canonical gene model 
format.  This script converts the output of gapmm2 into a more canonical gene model format.

EXAMPLE INPUT:

HiC_scaffold_9:52722461-53134832        gapmm2  cDNA_match      197149  197558  99.51   -       .       ID=gapmm2_1;Target=CNS.TRINITY_DN69598_c1_g1_i1 11 408
HiC_scaffold_10:39319993-41941512       gapmm2  cDNA_match      2520858 2521179 91.56   +       .       ID=gapmm2_2;Target=CNS.TRINITY_DN69600_c0_g4_i1 1 320
HiC_scaffold_10:39319993-41941512       gapmm2  cDNA_match      2520085 2520979 95.51   +       .       ID=gapmm2_3;Target=CNS.TRINITY_DN69600_c0_g1_i1 1 847
HiC_scaffold_14:31731738-32950677       gapmm2  cDNA_match      1109768 1110490 99.14   +       .       ID=gapmm2_5;Target=CNS.TRINITY_DN28092_c0_g1_i1 1 723
HiC_scaffold_14:31731738-32950677       gapmm2  cDNA_match      1112819 1112906 99.14   +       .       ID=gapmm2_5;Target=CNS.TRINITY_DN28092_c0_g1_i1 724 811
HiC_scaffold_14:31731738-32950677       gapmm2  cDNA_match      1113531 1113669 99.14   +       .       ID=gapmm2_5;Target=CNS.TRINITY_DN28092_c0_g1_i1 812 950
HiC_scaffold_14:31731738-32950677       gapmm2  cDNA_match      1119671 1120132 99.14   +       .       ID=gapmm2_5;Target=CNS.TRINITY_DN28092_c0_g1_i1 951 1393

This is perfectly valid GFF3, but some needs call for this to be transformed into
something which looks more like predicted gene models.

This script reads the native GFF output and makes it conform to the suggested canonical
gene model, defined here:

http://www.sequenceontology.org/gff3.shtml

For the input example above, this would be:

EXAMPLE OUTPUT:

HiC_scaffold_9:52722461-53134832        gapmm2  gene    197149  197558  .       -       .       ID=gapmm2_1;locus_tag=TRINITY_DN69598_c1_g1_i1
HiC_scaffold_9:52722461-53134832        gapmm2  mRNA    197149  197558  .       -       .       ID=gapmm2_1.mRNA;Parent=gapmm2_1;locus_tag=TRINITY_DN69598_c1_g1_i1
HiC_scaffold_9:52722461-53134832        gapmm2  CDS     197149  197558  .       -       .       ID=gapmm2_1.CDS;Parent=gapmm2_1.mRNA
HiC_scaffold_9:52722461-53134832        gapmm2  exon    197149  197558  .       -       .       ID=gapmm2_1.mRNA.exon1;Parent=gapmm2_1.mRNA
HiC_scaffold_10:39319993-41941512       gapmm2  gene    2520858 2521179 .       +       .       ID=gapmm2_2;locus_tag=TRINITY_DN69600_c0_g4_i1
HiC_scaffold_10:39319993-41941512       gapmm2  mRNA    2520858 2521179 .       +       .       ID=gapmm2_2.mRNA;Parent=gapmm2_2;locus_tag=TRINITY_DN69600_c0_g4_i1
HiC_scaffold_10:39319993-41941512       gapmm2  CDS     2520858 2521179 .       +       .       ID=gapmm2_2.CDS;Parent=gapmm2_2.mRNA
HiC_scaffold_10:39319993-41941512       gapmm2  exon    2520858 2521179 .       +       .       ID=gapmm2_2.mRNA.exon1;Parent=gapmm2_2.mRNA
HiC_scaffold_10:39319993-41941512       gapmm2  gene    2520085 2520979 .       +       .       ID=gapmm2_3;locus_tag=TRINITY_DN69600_c0_g1_i1
HiC_scaffold_10:39319993-41941512       gapmm2  mRNA    2520085 2520979 .       +       .       ID=gapmm2_3.mRNA;Parent=gapmm2_3;locus_tag=TRINITY_DN69600_c0_g1_i1
HiC_scaffold_10:39319993-41941512       gapmm2  CDS     2520085 2520979 .       +       .       ID=gapmm2_3.CDS;Parent=gapmm2_3.mRNA
HiC_scaffold_10:39319993-41941512       gapmm2  exon    2520085 2520979 .       +       .       ID=gapmm2_3.mRNA.exon1;Parent=gapmm2_3.mRNA
HiC_scaffold_14:31731738-32950677       gapmm2  gene    1109768 1120132 .       +       .       ID=gapmm2_5;locus_tag=TRINITY_DN28092_c0_g1_i1
HiC_scaffold_14:31731738-32950677       gapmm2  mRNA    1109768 1120132 .       +       .       ID=gapmm2_5.mRNA;Parent=gapmm2_5;locus_tag=TRINITY_DN28092_c0_g1_i1
HiC_scaffold_14:31731738-32950677       gapmm2  CDS     1109768 1110490 .       +       .       ID=gapmm2_5.CDS;Parent=gapmm2_5.mRNA
HiC_scaffold_14:31731738-32950677       gapmm2  CDS     1112819 1112906 .       +       .       ID=gapmm2_5.CDS;Parent=gapmm2_5.mRNA
HiC_scaffold_14:31731738-32950677       gapmm2  CDS     1113531 1113669 .       +       .       ID=gapmm2_5.CDS;Parent=gapmm2_5.mRNA
HiC_scaffold_14:31731738-32950677       gapmm2  CDS     1119671 1120132 .       +       .       ID=gapmm2_5.CDS;Parent=gapmm2_5.mRNA
HiC_scaffold_14:31731738-32950677       gapmm2  exon    1109768 1110490 .       +       .       ID=gapmm2_5.mRNA.exon1;Parent=gapmm2_5.mRNA
HiC_scaffold_14:31731738-32950677       gapmm2  exon    1112819 1112906 .       +       .       ID=gapmm2_5.mRNA.exon2;Parent=gapmm2_5.mRNA
HiC_scaffold_14:31731738-32950677       gapmm2  exon    1113531 1113669 .       +       .       ID=gapmm2_5.mRNA.exon3;Parent=gapmm2_5.mRNA
HiC_scaffold_14:31731738-32950677       gapmm2  exon    1119671 1120132 .       +       .       ID=gapmm2_5.mRNA.exon4;Parent=gapmm2_5.mRNA

NOTE:
The 2nd column will be made 'gapmm2' unless you pass a different value with the -s option.

The transcript ID from the input Target column will be used as the locus_tag for the gene and mRNA features.
'''

import argparse

from biocode import gff, things


def main():
    parser = argparse.ArgumentParser( description='Convert gapmm2 GFF file to canonical gene models')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to a GFF file created by gapmm2' )
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-s', '--source', type=str, required=False, default='gapmm2', help='Value to use for the 2nd (source) column' )
    args = parser.parse_args()

    assemblies = dict()
    current_assembly = None
    
    gene = None
    mRNA = None
    gene_fmin = None
    gene_fmax = None
    gene_strand = None

    ## Used for tracking the exon count for each gene (for ID purposes)
    exon_count_by_mRNA = dict()
    
    fout = open(args.output, mode='wt', encoding='utf-8')
    fout.write("##gff-version 3\n")

    for line in open(args.input):
        cols = line.split("\t")

        if len(cols) != 9:
            continue

        mol_id = cols[0]
        feat_type = cols[2]
        feat_id = gff.column_9_value(cols[8], 'ID')
        feat_target = gff.column_9_value(cols[8], 'Target')

        # If the target has a space in it, we only want the first part (the actual target ID)
        if feat_target is not None and ' ' in feat_target:
            feat_target = feat_target.split(' ')[0]

        # we expect all columns to be cDNA_match
        if feat_type != 'cDNA_match':
            raise Exception("ERROR: expected all columns to be of type 'cDNA_match' but found a {0}".format(feat_type))

        ## initialize this assembly if we haven't seen it yet
        if mol_id not in assemblies:
            assemblies[mol_id] = things.Assembly(id=mol_id)

        if gene is None or feat_id != gene.id:
            if gene is not None:
                # finish the previous one first
                mRNA.locate_on( target=current_assembly, fmin=gene_fmin, fmax=gene_fmax, strand=gene_strand )
                gene.locate_on( target=current_assembly, fmin=gene_fmin, fmax=gene_fmax, strand=gene_strand )
                gene.add_mRNA(mRNA)
                current_assembly.add_gene( gene )
                gene.print_as(fh=fout, source=args.source, format='gff3')

            # now start a new one
            gene = things.Gene(id=feat_id, locus_tag=feat_target)
            mRNA = things.mRNA(id="{0}.mRNA".format(feat_id), parent=gene, locus_tag=feat_target)
            exon_count_by_mRNA[mRNA.id] = 0
            
            gene_fmin = int(cols[3]) - 1
            gene_fmax = int(cols[4])
            gene_strand = cols[6]

        current_assembly = assemblies[mol_id]
            
        # each row is a new CDS/exon for the current mRNA
        CDS = things.CDS(id="{0}.CDS".format(feat_id), parent=mRNA.id)
        CDS.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6], phase='.' )
        mRNA.add_CDS(CDS)
        
        exon_count_by_mRNA[mRNA.id] += 1
        exon_id = "{0}.exon{1}".format(mRNA.id, exon_count_by_mRNA[mRNA.id])
        exon = things.Exon(id=exon_id, parent=mRNA.id)
        exon.locate_on( target=current_assembly, fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6] )
        mRNA.add_exon(exon)

        if int(cols[3]) - 1 < gene_fmin:
            gene_fmin = int(cols[3]) - 1

        if int(cols[4]) > gene_fmax:
            gene_fmax = int(cols[4])

    # don't orphan the last one
    if gene is not None:
        # finish the previous one first
        mRNA.locate_on( target=current_assembly, fmin=gene_fmin, fmax=gene_fmax, strand=gene_strand )
        gene.locate_on( target=current_assembly, fmin=gene_fmin, fmax=gene_fmax, strand=gene_strand )
        gene.add_mRNA(mRNA)
        current_assembly.add_gene( gene )
        gene.print_as(fh=fout, source=args.source, format='gff3')

if __name__ == '__main__':
    main()







