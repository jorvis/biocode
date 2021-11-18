#!/usr/bin/env python3

"""
This is a script to convert GTF flat files from Stringtie GTF to GFF3 format.

Tested with StringTie version 2.1.1

INPUT example:

contig_724	StringTie	transcript	40257	46234	1000	+	.	gene_id "STRG.4"; transcript_id "STRG.4.1"; cov "100.958267"; FPKM "5.136298"; TPM "25.765736";
contig_724	StringTie	exon	40257	40357	1000	+	.	gene_id "STRG.4"; transcript_id "STRG.4.1"; exon_number "1"; cov "104.381264";
contig_724	StringTie	exon	43165	46234	1000	+	.	gene_id "STRG.4"; transcript_id "STRG.4.1"; exon_number "2"; cov "100.845650";
contig_724	StringTie	transcript	40271	46384	1000	+	.	gene_id "STRG.4"; transcript_id "STRG.4.2"; cov "60.526722"; FPKM "3.079324"; TPM "15.447130";
contig_724	StringTie	exon	40271	40336	1000	+	.	gene_id "STRG.4"; transcript_id "STRG.4.2"; exon_number "1"; cov "59.643803";
contig_724	StringTie	exon	43165	46384	1000	+	.	gene_id "STRG.4"; transcript_id "STRG.4.2"; exon_number "2"; cov "60.544823";

OUTPUT example (--export_mode=model):

contig_724      StringTie       gene    40257   46234   .       +       .       ID=STRG.4
contig_724      StringTie       mRNA    40257   46234   .       +       .       ID=STRG.4.1;Parent=STRG.4
contig_724      StringTie       CDS     40257   40357   .       +       0       ID=STRG.4.1.CDS.1;Parent=STRG.4.1
contig_724      StringTie       CDS     43165   46234   .       +       1       ID=STRG.4.1.CDS.2;Parent=STRG.4.1
contig_724      StringTie       exon    40257   40357   .       +       .       ID=STRG.4.1.exon.1;Parent=STRG.4.1
contig_724      StringTie       exon    43165   46234   .       +       .       ID=STRG.4.1.exon.2;Parent=STRG.4.1
contig_724      StringTie       mRNA    40271   46384   .       +       .       ID=STRG.4.2;Parent=STRG.4
contig_724      StringTie       CDS     40271   40336   .       +       0       ID=STRG.4.2.CDS.1;Parent=STRG.4.2
contig_724      StringTie       CDS     43165   46384   .       +       0       ID=STRG.4.2.CDS.2;Parent=STRG.4.2
contig_724      StringTie       exon    40271   40336   .       +       .       ID=STRG.4.2.exon.1;Parent=STRG.4.2
contig_724      StringTie       exon    43165   46384   .       +       .       ID=STRG.4.2.exon.2;Parent=STRG.4.2

Caveats:
- Notice that the input defined above describes multiple isoforms for the same gene.  Some downstream
  GFF3 parsers may not support this (though they should.)
- Currently not exporting any of the other attributes from column 9 of the input files


Author: Joshua Orvis (jorvis AT gmail)
"""

import argparse
import sys
from collections import defaultdict

from biocode import gff, things


def main():
    parser = argparse.ArgumentParser( description='A GTF -> GFF3 conversion script for StringTie output')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input GTF file' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output GFF file to be created' )
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
    
    current_match = None

    rna_count_by_gene = defaultdict(int)
    exon_count_by_RNA = defaultdict(int)

    for line in open(args.input_file, "r"):
        cols = line.split("\t")

        if len(cols) != 9:
            print("SKIPPING: {0}".format(line))
            continue
        
        mol_id = cols[0]

        if mol_id not in assemblies:
            assemblies[mol_id] = things.Assembly(id=mol_id)

        current_assembly = assemblies[mol_id]
        ftype  = cols[2]
        fmin = int(cols[3]) - 1
        fmax = int(cols[4])
        strand = cols[6]
        col9 = cols[8]

        # this makes it look like GFF column 9 so I can use biocodeutils.column_9_value(str, key)
        col9 = col9.replace(' "', '="')
        gene_id       = gff.column_9_value(col9, 'gene_id').replace('"', '')
        transcript_id = gff.column_9_value(col9, 'transcript_id').replace('"', '')
        cov = gff.column_9_value(col9, 'cov').replace('"', '')
        
        if ftype == 'transcript':
            if current_gene is not None and current_gene.id != gene_id:
                gene.print_as(fh=ofh, source='StringTie', format='gff3')

            if current_gene is None or current_gene.id != gene_id:
                gene = things.Gene(id=gene_id)
                gene.locate_on( target=current_assembly, fmin=fmin, fmax=fmax, strand=strand )
                current_gene = gene

            mRNA = things.mRNA(id=transcript_id, parent=current_gene)
            mRNA.locate_on( target=current_assembly, fmin=fmin, fmax=fmax, strand=strand )
            gene.add_mRNA(mRNA)
            current_RNA = mRNA
            exon_count_by_RNA[transcript_id] = 0
            current_CDS_phase = 0

        elif ftype == 'exon':
            exon_number = gff.column_9_value(col9, 'exon_number').replace('"', '')
            exon_count_by_RNA[transcript_id] += 1

            cds_id = "{0}.CDS.{1}".format( current_RNA.id, exon_count_by_RNA[current_RNA.id] )
            CDS = things.CDS(id=cds_id, parent=current_RNA)
            CDS.locate_on( target=current_assembly, fmin=fmin, fmax=fmax, strand=strand, phase=current_CDS_phase )
            current_RNA.add_CDS(CDS)

             # calculate the starting phase for the next CDS feature (in case there is one)
            current_CDS_phase = 3 - (((fmax - fmin) - current_CDS_phase) % 3)
            if current_CDS_phase == 3:
                current_CDS_phase = 0

            exon_id = "{0}.exon.{1}".format( current_RNA.id, exon_count_by_RNA[current_RNA.id] )
            exon = things.Exon(id=exon_id, parent=current_RNA)
            exon.locate_on( target=current_assembly, fmin=fmin, fmax=fmax, strand=strand )
            current_RNA.add_exon(exon)
                
    # don't forget to do the last gene, if there were any
    if current_gene is not None:
        gene.print_as(fh=ofh, source='StringTie', format='gff3')


if __name__ == '__main__':
    main()







