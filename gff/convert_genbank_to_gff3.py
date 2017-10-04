#!/usr/bin/env python3

"""
This is a script to convert GenBank flat files to GFF3 format with a specific focus on
initially maintaining as much structural annotation as possible, then expanding into
functional annotation support.

This is not guaranteed to convert all features, but warnings will be printed wherever possible
for features which aren't included.

Currently supported:
  Structural features:  gene, CDS, mRNA, tRNA, rRNA
  Annotations: primary identifiers, gene product name

This is written to handle multi-entry GBK files

Caveats:
- Because the GBK flatfile format doesn't explicitly model parent/child features, this script
  links them using the expected format convention of shared /locus_tag entries for each feature
  of the gene graph (gene, mRNA, CDS)
- It has only been tested with prokaryotic (non-spliced) genes

Author: Joshua Orvis (jorvis AT gmail)
"""

import argparse
import sys
from collections import defaultdict

from Bio import SeqIO
from biocode import annotation, things, utils


def main():
    parser = argparse.ArgumentParser( description='Convert GenBank flat files to GFF3 format')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input GBK file' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output GFF file to be created' )
    parser.add_argument('--with_fasta', dest='fasta', action='store_true', help='Include the FASTA section with genomic sequence at end of file.  (default)' )
    parser.add_argument('--no_fasta', dest='fasta', action='store_false' )
    parser.set_defaults(fasta=True)
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    ofh = sys.stdout
    if args.output_file is not None:
        ofh = open(args.output_file, 'wt')

    ofh.write("##gff-version 3\n")

    assemblies = dict()
    current_assembly = None
    current_gene = None
    current_RNA = None

    rna_count_by_gene = defaultdict(int)
    exon_count_by_RNA = defaultdict(int)

    seqs_pending_writes = False

    features_skipped_count = 0

    # each gb_record is a SeqRecord object
    for gb_record in SeqIO.parse(open(args.input_file, "r"), "genbank"):
        mol_id = gb_record.name

        if mol_id not in assemblies:
            assemblies[mol_id] = things.Assembly(id=mol_id)

        if len(str(gb_record.seq)) > 0:
            seqs_pending_writes = True
            assemblies[mol_id].residues = str(gb_record.seq)
            assemblies[mol_id].length = len(str(gb_record.seq))

        current_assembly = assemblies[mol_id]
            
        # each feat is a SeqFeature object
        for feat in gb_record.features:
            #print(feat)
            fmin = int(feat.location.start)
            fmax = int(feat.location.end)

            if feat.location.strand == 1:
                strand = '+'
            elif feat.location.strand == -1:
                strand = '-'
            else:
                raise Exception("ERROR: unstranded feature encountered: {0}".format(feat))

            #print("{0} located at {1}-{2} strand:{3}".format( locus_tag, fmin, fmax, strand ) )
            if feat.type == 'source':
                continue
            
            if feat.type == 'gene':
                # print the previous gene (if there is one)
                if current_gene is not None:
                    gene.print_as(fh=ofh, source='GenBank', format='gff3')
                
                locus_tag = feat.qualifiers['locus_tag'][0]
                gene = things.Gene(id=locus_tag, locus_tag=locus_tag)
                gene.locate_on( target=current_assembly, fmin=fmin, fmax=fmax, strand=strand )
                current_gene = gene
                current_RNA = None

            elif feat.type == 'mRNA':
                locus_tag = feat.qualifiers['locus_tag'][0]
                rna_count_by_gene[locus_tag] += 1
                feat_id = "{0}.mRNA.{1}".format( locus_tag, rna_count_by_gene[locus_tag] )
                
                mRNA = things.mRNA(id=feat_id, parent=current_gene, locus_tag=locus_tag)
                mRNA.locate_on( target=current_assembly, fmin=fmin, fmax=fmax, strand=strand )
                gene.add_mRNA(mRNA)
                current_RNA = mRNA

                if feat_id in exon_count_by_RNA:
                    raise Exception( "ERROR: two different RNAs found with same ID: {0}".format(feat_id) )
                else:
                    exon_count_by_RNA[feat_id] = 0

            elif feat.type == 'tRNA':
                locus_tag = feat.qualifiers['locus_tag'][0]
                rna_count_by_gene[locus_tag] += 1
                feat_id = "{0}.tRNA.{1}".format( locus_tag, rna_count_by_gene[locus_tag] )
                
                tRNA = things.tRNA(id=feat_id, parent=current_gene)
                tRNA.locate_on( target=current_assembly, fmin=fmin, fmax=fmax, strand=strand )
                gene.add_tRNA(tRNA)
                current_RNA = tRNA

                if feat_id in exon_count_by_RNA:
                    raise Exception( "ERROR: two different RNAs found with same ID: {0}".format(feat_id) )
                else:
                    exon_count_by_RNA[feat_id] = 0

            elif feat.type == 'rRNA':
                locus_tag = feat.qualifiers['locus_tag'][0]
                rna_count_by_gene[locus_tag] += 1
                feat_id = "{0}.rRNA.{1}".format( locus_tag, rna_count_by_gene[locus_tag] )
                
                rRNA = things.rRNA(id=feat_id, parent=current_gene)
                rRNA.locate_on( target=current_assembly, fmin=fmin, fmax=fmax, strand=strand )
                gene.add_rRNA(rRNA)
                current_RNA = rRNA

                if feat_id in exon_count_by_RNA:
                    raise Exception( "ERROR: two different RNAs found with same ID: {0}".format(feat_id) )
                else:
                    exon_count_by_RNA[feat_id] = 0

            elif feat.type == 'CDS':
                locus_tag = feat.qualifiers['locus_tag'][0]
                # If processing a prokaryotic GBK, we'll encounter CDS before mRNA, so we have to
                #  manually make one
                if current_RNA is None:
                    feat_id = "{0}.mRNA.{1}".format( locus_tag, rna_count_by_gene[locus_tag] )
                    mRNA = things.mRNA(id=feat_id, parent=current_gene)
                    mRNA.locate_on( target=current_assembly, fmin=fmin, fmax=fmax, strand=strand )
                    gene.add_mRNA(mRNA)
                    current_RNA = mRNA

                    if 'product' in feat.qualifiers:
                        product = feat.qualifiers['product'][0]
                    else:
                        product = None

                    if 'gene' in feat.qualifiers:
                        gene_symbol = feat.qualifiers['gene'][0]
                    else:
                        gene_symbol = None
                        
                    annot = annotation.FunctionalAnnotation(product_name=product, gene_symbol=gene_symbol)

                    if 'db_xref' in feat.qualifiers:
                        for dbxref in feat.qualifiers['db_xref']:
                            annot.add_dbxref(dbxref)
                    
                    polypeptide_id = "{0}.polypeptide.{1}".format( locus_tag, rna_count_by_gene[locus_tag] )
                    polypeptide = things.Polypeptide(id=polypeptide_id, parent=mRNA, annotation=annot)
                    mRNA.add_polypeptide(polypeptide)
                
                exon_count_by_RNA[current_RNA.id] += 1
                cds_id = "{0}.CDS.{1}".format( current_RNA.id, exon_count_by_RNA[current_RNA.id] )
                current_CDS_phase = 0
                
                for loc in feat.location.parts:
                    subfmin = int(loc.start)
                    subfmax = int(loc.end)
                    
                    CDS = things.CDS(id=cds_id, parent=current_RNA)
                    CDS.locate_on( target=current_assembly, fmin=subfmin, fmax=subfmax, strand=strand, phase=current_CDS_phase )
                    current_RNA.add_CDS(CDS)

                    # calculate the starting phase for the next CDS feature (in case there is one)
                    # 0 + 6 = 0     TTGCAT
                    # 0 + 7 = 2     TTGCATG
                    # 1 + 6 = 1     TTGCAT
                    # 2 + 7 = 1     TTGCATG
                    # general: 3 - ((length - previous phase) % 3)
                    current_CDS_phase = 3 - (((subfmax - subfmin) - current_CDS_phase) % 3)
                    if current_CDS_phase == 3:
                        current_CDS_phase = 0

                    exon_id = "{0}.exon.{1}".format( current_RNA.id, exon_count_by_RNA[current_RNA.id] )
                    exon = things.Exon(id=exon_id, parent=current_RNA)
                    exon.locate_on( target=current_assembly, fmin=subfmin, fmax=subfmax, strand=strand )
                    current_RNA.add_exon(exon)
                    exon_count_by_RNA[current_RNA.id] += 1
                
            else:
                print("WARNING: The following feature was skipped:\n{0}".format(feat))
                features_skipped_count += 1

    # don't forget to do the last gene, if there were any
    if current_gene is not None:
        gene.print_as(fh=ofh, source='GenBank', format='gff3')

    if args.fasta is True:
        if seqs_pending_writes is True:
            ofh.write("##FASTA\n")
            for assembly_id in assemblies:
                ofh.write(">{0}\n".format(assembly_id))
                ofh.write("{0}\n".format(utils.wrapped_fasta(assemblies[assembly_id].residues)))

    if features_skipped_count > 0:
        print("Warning: {0} unsupported feature types were skipped".format(features_skipped_count))

if __name__ == '__main__':
    main()







