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

ARAGORN: Optional. tRNA, mtRNA, and tmRNA prediction.  Not actually GFF3 but convenient
         to include here.  Assumes you ARAGORN has been run with the -w option

"""

import argparse
import os
import re
import uuid

from biocode import annotation, gff, things, utils

def main():
    parser = argparse.ArgumentParser( description='Creates a single GFF from the output of a few different model prediction tools (coding and non-coding)')

    ## output file to be written
    parser.add_argument('-m', '--model_gff', type=str, required=True, help='Input (pass-through) GFF file' )
    parser.add_argument('-o', '--output_gff', type=str, required=False, help='Output file to be written.  Default=STDOUT' )
    parser.add_argument('-b', '--barrnap_gff', type=str, required=False, help='GFF file from Barrnap prediction' )
    parser.add_argument('-g', '--genomic_fasta', type=str, required=True, help='Source genomic FASTA file' )
    parser.add_argument('-a', '--aragorn_out', type=str, required=False, help='Raw output file (with -w) from ARAGORN prediction' )
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.model_gff)
    utils.add_assembly_fasta(assemblies, args.genomic_fasta)

    if args.barrnap_gff:
        add_barrnap_features(assemblies, features, args.barrnap_gff)

    if args.aragorn_out:
        add_aragorn_features(assemblies, features, args.aragorn_out)

    with open(args.output_gff, 'wt') as f:
        gff.print_gff3_from_assemblies(ofh=f, assemblies=assemblies)


def add_aragorn_features(assemblies, features, aragorn_file):
    current_assembly_id = None

    for line in open(aragorn_file):
        line = line.rstrip()

        if line.startswith('>'):
            m = re.match('>(\S+)', line)
            current_assembly_id = m.group(1)

            if current_assembly_id not in assemblies:
                assembly = things.Assembly(id=current_assembly_id, residues='')
                assemblies[current_assembly_id] = assembly

        else:
            cols = line.split()

            if len(cols) == 5:
                if cols[1].startswith('tRNA'):
                    feat_type = 'tRNA'
                elif cols[1].startswith('tmRNA'):
                    feat_type = 'tmRNA'
                elif cols[1].startswith('mtRNA'):
                    feat_type = 'mtRNA'
                else:
                    raise Exception("Unexpected type in ARAGORN, column value: {0}".format(cols[1]))
                    
                feat_base = "{0}_{1}".format(feat_type, uuid.uuid4())
                gene_id = "{0}_gene".format(feat_base)
                RNA_id = "{0}_{1}".format(feat_base, feat_type)

                m = re.match('(c*)\[(\d+),(\d+)\]', cols[2])
                if m:
                    rfmin = int(m.group(2)) - 1
                    rfmax = int(m.group(3))

                    # For predictions spanning the origin of circular molecules, fmax needs to be adjusted
                    #  https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
                    # Sub-heading: Circular Genomes
                    #
                    # Example input lines:
                    #  62  tRNA-Met                 c[2764659,18]      29      (cat)
                    #  1   tRNA-Ile                  [2697442,47]      35      (gat)
                    if rfmax < rfmin:
                        assemblies[current_assembly_id].is_circular = True
                        rfmax = rfmin + rfmax + 1

                    if m.group(1):
                        rstrand = -1
                    else:
                        rstrand = 1
                else:
                    raise Exception("ERROR: unexpected coordinate format: {0}".format(cols[2]))
                
                current_assembly = assemblies[current_assembly_id]
                gene = things.Gene(id=gene_id)
                gene.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
                features[gene_id] = gene
                current_assembly.add_gene(gene)

                if feat_type == 'tRNA':
                    RNA = things.tRNA(id=RNA_id, parent=gene, anticodon=cols[4][1:4].upper())
                    gene.add_tRNA(RNA)
                else:
                    RNA = things.tmRNA(id=RNA_id, parent=gene)
                    gene.add_tmRNA(RNA)

                RNA.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
                RNA.annotation = annotation.FunctionalAnnotation(product_name=cols[1])
                features[RNA_id] = RNA
    

def add_barrnap_features(assemblies, features, barrnap_gff):
    for line in open(barrnap_gff):
        if line.startswith('#'):
            continue

        cols = line.split("\t")

        if len(cols) == 9:
            if cols[0] in assemblies:
                current_assembly = assemblies[cols[0]]
            else:
                current_assembly = things.Assembly(id=cols[0], residues='')
                assemblies[cols[0]] = current_assembly

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







