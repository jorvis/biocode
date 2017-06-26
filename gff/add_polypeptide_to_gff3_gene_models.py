#!/usr/bin/env python3

'''
The biocode libraries primarly work with functional annotation attached to the polypeptide feature
of the canonical gene model.  Many tools, though, produce gff with gene/mRNA/CDS/exon features but
no polypeptide entry.  

This script checks the children of each mRNA and creates a polypeptide entry (using the CDS
coordinate span as a guide) for each which doesn't already have one.  If a there is a polypeptide
already, that mRNA is skipped.

Example input:

    unitig_50|quiver        .       gene    7189    7643    .       -       .       ID=9531_g
    unitig_50|quiver        .       mRNA    7189    7643    .       -       .       ID=9531_g.mRNA.9531_t;Parent=9531_g
    unitig_50|quiver        .       CDS     7189    7534    .       -       1       ID=9531_g.cds.1;Parent=9531_g.mRNA.9531_t
    unitig_50|quiver        .       CDS     7591    7643    .       -       0       ID=9531_g.cds.1;Parent=9531_g.mRNA.9531_t
    unitig_50|quiver        .       exon    7189    7534    .       -       .       ID=9531_g.exon.1;Parent=9531_g.mRNA.9531_t
    unitig_50|quiver        .       exon    7591    7643    .       -       .       ID=9531_g.exon.2;Parent=9531_g.mRNA.9531_t

And after this script:

    unitig_50|quiver        .       gene    7189    7643    .       -       .       ID=9531_g
    unitig_50|quiver        .       mRNA    7189    7643    .       -       .       ID=9531_g.mRNA.9531_t;Parent=9531_g
    unitig_50|quiver        .       CDS     7189    7534    .       -       1       ID=9531_g.cds.1;Parent=9531_g.mRNA.9531_t
    unitig_50|quiver        .       CDS     7591    7643    .       -       0       ID=9531_g.cds.1;Parent=9531_g.mRNA.9531_t
    unitig_50|quiver        .       exon    7189    7534    .       -       .       ID=9531_g.exon.1;Parent=9531_g.mRNA.9531_t
    unitig_50|quiver        .       exon    7591    7643    .       -       .       ID=9531_g.exon.2;Parent=9531_g.mRNA.9531_t
    unitig_50|quiver        .       polypeptide    7189    7643    .       -       .       ID=9531_g.exon.2.polypeptide;Parent=9531_g.mRNA.9531_t

Follow the GFF3 specification!

Author:  Joshua Orvis
'''

import argparse

from biocode import utils, gff, things

def main():
    parser = argparse.ArgumentParser( description='Adds polypeptide rows in GFF3 for gene models which are missing them')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-g', '--genome_fasta', type=str, required=False, help='Optional.  You must specify this unless the FASTA sequences for the molecules are embedded in the GFF')
    parser.add_argument('-s', '--source', type=str, required=False, default='.', help='Optional.  Sets the value for column 2 in all rows.  Default = .' )
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_file)

    fout = open(args.output_file, mode='wt', encoding='utf-8')
    fout.write("##gff-version 3\n")

    # deal with the FASTA file if the user passed one
    if args.genome_fasta is not None:
        process_assembly_fasta(assemblies, args.genome_fasta)

    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            for mRNA in gene.mRNAs():
                check_and_add_polypeptide(mRNA)

            gene.print_as(fh=fout, source=args.source, format='gff3')

    fasta_header_written = False

    for assembly_id in assemblies:
        if assemblies[assembly_id].length is not None and assemblies[assembly_id].length > 0:
            if fasta_header_written is False:
                fout.write("##FASTA\n")
                fasta_header_written = True

            fout.write(">{0}\n".format(assemblies[assembly_id].id) )
            fout.write("{0}\n".format(utils.wrapped_fasta(assemblies[assembly_id].residues)))

def check_and_add_polypeptide(mRNA):
    polypeptides = mRNA.polypeptides()
    if len(polypeptides) == 0:
        polypeptide_id = "{0}.polypeptide".format(mRNA.id)
        polypeptide = things.Polypeptide(id=polypeptide_id, parent=mRNA.id)
        mRNA.add_polypeptide(polypeptide)

def process_assembly_fasta(mols, fasta_file):
    fasta_seqs = utils.fasta_dict_from_file(fasta_file)

    for mol_id in mols:
        # check if the FASTA file provides sequence for this
        if mol_id in fasta_seqs:
            mol = mols[mol_id]
            mol.residues = fasta_seqs[mol_id]['s']
            mol.length   = len(mol.residues)

if __name__ == '__main__':
    main()






