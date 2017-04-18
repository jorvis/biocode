#!/usr/bin/env python3

'''
This script reads a GFF3 file and FASTA file (or FASTA embedded in the GFF) and
checks the ends of CDS features for start/stop codons.

This is splice-aware, and works for prokaryotic or eukaryotic models.

Example input gene model:
AAGK01000001      .       gene    128156  128682  .       -       .       ID=210D90E7BF533A9A08EF22D9F1F1CED9_1_2;locus_tag=TpMuguga_01g00060
AAGK01000001      .       mRNA    128156  128682  .       -       .       ID=0D5133EC9A534184443BC6C9D05E565B;Parent=210D90E7BF533A9A08EF22D9F1F1CED9_1_2;locus_tag=TpMuguga_01g00060
AAGK01000001      .       CDS     128349  128616  .       -       0       ID=5BBE8CD7940726C5808483FC5EC65D28;Parent=0D5133EC9A534184443BC6C9D05E565B
AAGK01000001      .       CDS     128156  128278  .       -       0       ID=5BBE8CD7940726C5808483FC5EC65D28;Parent=0D5133EC9A534184443BC6C9D05E565B
AAGK01000001      .       exon    128156  128278  .       -       .       ID=3E36768E591D34F98C62893C736848FB;Parent=0D5133EC9A534184443BC6C9D05E565B
AAGK01000001      .       exon    128349  128682  .       -       .       ID=3A9F94FDF76013E7216D47419026FDF1;Parent=0D5133EC9A534184443BC6C9D05E565B
AAGK01000001      .       polypeptide     128156  128682  .       -       .       ID=0D5133EC9A534184443BC6C9D05E565B_p;Parent=0D5133EC9A534184443BC6C9D05E565B;product_name=hypothetical protein

Marking the features as partial in GFF follows these guidelines:
   http://www.sequenceontology.org/so_wiki/index.php/PartialFeatures

In short, partiality is maintained between related features.  If a CDS is partial then the
corresponding exon, mRNA and gene features in the graph are also marked as partial.  The
different key/values added in the 9th column for partiality are:

   Partial=5prime
   Partial=3prime
   Partial=5prime,3prime

A short report is printed on STDOUT when completed.  Example:

   Genes marked as 5' partial: 15
   Genes marked as 3' partial: 67

Codons allowed:
   Start: ATG, GTG, TTG
   Stop : TAG, TAA, TGA

A use case for needing this is in downstream exporting of the GFF for submission
to GenBank.  Their partiality requirements:

   https://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation/#Partialcodingregionsinincompletegenomes

Follow the GFF3 specification!

Author:  Joshua Orvis
'''

import argparse
from biocode import utils, gff

def main():
    parser = argparse.ArgumentParser( description='Extends GFF gene models to the first in-frame stop')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-g', '--genome_fasta', type=str, required=False, help='Optional.  You must specify this unless the FASTA sequences for the molecules are embedded in the GFF')
    parser.add_argument('-o', '--output_gff', type=str, required=False, help='Optional.  Writes an output GFF3 file with CDS (and containing features) extended to nearest stop')
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_file)

    # deal with the FASTA file if the user passed one
    if args.genome_fasta is not None:
        utils.add_assembly_fasta(assemblies, args.genome_fasta)

    start_codons = ['ATG', 'GTG', 'TTG']
    stop_codons = ['TAG', 'TAA', 'TGA']

    newly_marked_5prime_partial = 0
    newly_marked_3prime_partial = 0

    for assembly_id in sorted(assemblies):
        for gene in sorted(assemblies[assembly_id].genes()):
            gene_loc = gene.location_on(assemblies[assembly_id])

            for mRNA in gene.mRNAs():
                mRNA_loc = mRNA.location_on(assemblies[assembly_id])
                coding_seq = mRNA.get_CDS_residues()
                translation = utils.translate(coding_seq)

                if not translation.endswith('*'):
                    newly_marked_3prime_partial += 1
                    CDSs = sorted(mRNA.CDSs())

                    if mRNA_loc.strand == 1:
                        mRNA_loc.fmax_partial = True
                        CDSs[-1].location_on(assemblies[assembly_id]).fmax_partial = True
                        gene_loc.fmax_partial = True

                        # The exon is tricky, as there's no direct link between the CDS fragment
                        #  and the corresponding exon.  The assumption here is that there won't
                        #  be terminal non-coding exons if the CDS is partial.
                        mRNA.exons()[-1].location_on(assemblies[assembly_id]).fmax_partial = True

                    else:
                        mRNA_loc.fmin_partial = True
                        gene_loc.fmin_partial = True
                        CDSs[0].location_on(assemblies[assembly_id]).fmin_partial = True
                        mRNA.exons()[0].location_on(assemblies[assembly_id]).fmin_partial = True

                start_codon = coding_seq[0:3].upper().replace('U', 'T')
                if start_codon not in start_codons:
                    newly_marked_5prime_partial += 1
                    CDSs = sorted(mRNA.CDSs())

                    if mRNA_loc.strand == 1:
                        mRNA_loc.fmin_partial = True
                        CDSs[0].location_on(assemblies[assembly_id]).fmin_partial = True
                        gene_loc.fmin_partial = True

                        # The exon is tricky, as there's no direct link between the CDS fragment
                        #  and the corresponding exon.  The assumption here is that there won't
                        #  be terminal non-coding exons if the CDS is partial.
                        mRNA.exons()[0].location_on(assemblies[assembly_id]).fmin_partial = True

                    else:
                        mRNA_loc.fmax_partial = True
                        gene_loc.fmax_partial = True
                        CDSs[-1].location_on(assemblies[assembly_id]).fmax_partial = True
                        mRNA.exons()[-1].location_on(assemblies[assembly_id]).fmax_partial = True

    print ("Genes marked as 5' partial: {0}".format(newly_marked_5prime_partial))
    print ("Genes marked as 3' partial: {0}".format(newly_marked_3prime_partial))

    ofh = open(args.output_gff, 'wt')
    gff.print_gff3_from_assemblies(assemblies=assemblies, ofh=ofh)


if __name__ == '__main__':
    main()






