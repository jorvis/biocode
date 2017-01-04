#!/usr/bin/env python3

'''
By default PASA will create an output file named like ?.pasa_assemblies.gff3, which
looks like this:

EXAMPLE INPUT:
ChromosomeI  assembler-pasa_foo      cDNA_match      14050   14471   .       +       .       ID=align_113453;Target=asmbl_9 1 422 +
ChromosomeI  assembler-pasa_foo      cDNA_match      14480   14749   .       +       .       ID=align_113454;Target=asmbl_10 1 270 +
ChromosomeI  assembler-pasa_foo      cDNA_match      14771   15824   .       +       .       ID=align_113454;Target=asmbl_10 271 1324 +

This is perfectly valid GFF3, but some needs call for this to be transformed into
something which looks more like predicted gene models.

This script reads the native GFF output and makes it conform to the suggested canonical
gene model, defined here:

http://www.sequenceontology.org/gff3.shtml

For the input example above, this would be:

EXAMPLE OUTPUT:
ChromosomeI  PASA    gene    14050   14471   .       +       .       ID=align_113453
ChromosomeI  PASA    mRNA    14050   14471   .       +       .       ID=align_113453.mRNA;Parent=align_113453
ChromosomeI  PASA    CDS     14050   14471   .       +       .       ID=align_113453.CDS;Parent=align_113453.mRNA
ChromosomeI  PASA    exon    14050   14471   .       +       .       ID=align_113453.mRNA.exon1;Parent=align_113453.mRNA
ChromosomeI  PASA    gene    14480   15824   .       +       .       ID=align_113454
ChromosomeI  PASA    mRNA    14480   15824   .       +       .       ID=align_113454.mRNA;Parent=align_113454
ChromosomeI  PASA    CDS     14480   14749   .       +       .       ID=align_113454.CDS;Parent=align_113454.mRNA
ChromosomeI  PASA    CDS     14771   15824   .       +       .       ID=align_113454.CDS;Parent=align_113454.mRNA
ChromosomeI  PASA    exon    14480   14749   .       +       .       ID=align_113454.mRNA.exon1;Parent=align_113454.mRNA
ChromosomeI  PASA    exon    14771   15824   .       +       .       ID=align_113454.mRNA.exon2;Parent=align_113454.mRNA


NOTE:
The 2nd column will be made 'PASA' unless you pass a different value with the -s option.

'''

import argparse

from biocode import gff, things


def main():
    parser = argparse.ArgumentParser( description='Convert PASA GFF file to canonical gene models')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to a GFF file created by PASA' )
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-s', '--source', type=str, required=False, default='PASA', help='Value to use for the 2nd (source) column' )
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
            gene = things.Gene(id=feat_id)
            mRNA = things.mRNA(id="{0}.mRNA".format(feat_id), parent=gene)
            exon_count_by_mRNA[mRNA.id] = 0
            
            gene_fmin = int(cols[3]) - 1
            gene_fmax = int(cols[4])
            gene_strand = cols[6]

        current_assembly = assemblies[mol_id]
            
        # each row is a new CDS/exon for the current mRNA
        CDS = things.CDS(id="{0}.CDS".format(feat_id), parent=mRNA.id)
        # FIX THIS PHASE
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







