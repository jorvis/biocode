#!/usr/bin/env python3

"""

Pipelines like GALES produce GFF3 models for protein-coding genes, then secondary
files for non-coding features like tRNAs and rRNAs.

This script combines these and does a bit of correction where the tools don't 
precisely produce spec-GFF3.

Current inputs:

models: Any spec-compliant GFF3 describing the protein-coding genes (or more).  These
        features are just retained and passed forward.

barrnap: Optional. rRNA prediction 'gff' file from Barrnap

tRNAscan-SE: Support coming soon

"""

import argparse
import os
import uuid

from biocode import gff, things, utils


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-m', '--model_gff', type=str, required=True, help='Input (pass-through) GFF file' )
    parser.add_argument('-o', '--output_gff', type=str, required=False, help='Output file to be written.  Default=STDOUT' )
    parser.add_argument('-b', '--barrnap_gff', type=str, required=False, help='GFF file from Barrnap prediction' )
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.model_gff)

    if args.barrnap_gff:
        add_barrnap_features(assemblies, features, args.barrnap_gff)

    utils.serialize_gff3(path=args.output_gff, assemblies=assemblies, features=features)


def add_barrnap_features(assemblies, features, barrnap_gff):
    for line in open(barrnap_gff):
        if line.startswith('#'):
            continue

        cols = line.split("\t")

        if len(cols) == 9:
            if cols[0] in assemblies:
                current_assembly = assemblies[cols[0]]
            else:
                raise Exception("ERROR: encountered assembly ID ({0}) in Barrnap file that was not found in the model GFF3 file".format(cols[0]))

            if cols[2] == 'rRNA':
                atts = gff.column_9_dict(cols[8])
                feat_base = "rRNA_{0}".format(uuid.uuid4())
                gene_id = "{0}_gene".format(feat_base)
                rRNA_id = "{0}_rRNA".format(feat_base)

                rfmin = int(cols[3]) - 1
                rfmax = int(cols[4])

                if cols[6] == '-':
                    rstrand = -1
                elif cols[6] == '+':
                    rstrand = 1
                else:
                    rstrand = 0

                gene = things.Gene(id=gene_id)
                gene.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
                features[gene_id] = gene
                current_assembly.add_gene(gene)

                rRNA = things.rRNA(id=rRNA_id, parent=gene)
                rRNA.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
                gene.add_rRNA(rRNA)
                rRNA.annotation = gff.parse_annotation_from_column_9(cols[8])
                features[rRNA_id] = rRNA
    


if __name__ == '__main__':
    main()







