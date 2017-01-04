#!/usr/bin/env python3

"""
This is a script to convert GTF flat files from Cufflinks to GFF3 format.

INPUT example:

tp.assembly.567468735.1	Cufflinks	transcript	20970	21966	1000	-	.	gene_id "CUFF.4"; transcript_id "CUFF.4.1"; FPKM "41.3914553346"; frac "1.000000"; conf_lo "38.477359"; conf_hi "44.305552"; cov "218.850840";
tp.assembly.567468735.1	Cufflinks	exon	20970	21966	1000	-	.	gene_id "CUFF.4"; transcript_id "CUFF.4.1"; exon_number "1"; FPKM "41.3914553346"; frac "1.000000"; conf_lo "38.477359"; conf_hi "44.305552"; cov "218.850840";
tp.assembly.567468735.1	Cufflinks	transcript	22213	23755	1000	-	.	gene_id "CUFF.5"; transcript_id "CUFF.5.1"; FPKM "143.4241721939"; frac "0.561266"; conf_lo "139.073452"; conf_hi "147.774892"; cov "755.889535";
tp.assembly.567468735.1	Cufflinks	exon	22213	22259	1000	-	.	gene_id "CUFF.5"; transcript_id "CUFF.5.1"; exon_number "1"; FPKM "143.4241721939"; frac "0.561266"; conf_lo "139.073452"; conf_hi "147.774892"; cov "755.889535";
tp.assembly.567468735.1	Cufflinks	exon	22294	23755	1000	-	.	gene_id "CUFF.5"; transcript_id "CUFF.5.1"; exon_number "2"; FPKM "143.4241721939"; frac "0.561266"; conf_lo "139.073452"; conf_hi "147.774892"; cov "755.889535";
tp.assembly.567468735.1	Cufflinks	transcript	22213	28599	124	-	.	gene_id "CUFF.5"; transcript_id "CUFF.5.2"; FPKM "17.7850894771"; frac "0.322877"; conf_lo "16.994111"; conf_hi "18.576068"; cov "93.732896";
tp.assembly.567468735.1	Cufflinks	exon	22213	23755	124	-	.	gene_id "CUFF.5"; transcript_id "CUFF.5.2"; exon_number "1"; FPKM "17.7850894771"; frac "0.322877"; conf_lo "16.994111"; conf_hi "18.576068"; cov "93.732896";
tp.assembly.567468735.1	Cufflinks	exon	23887	25044	124	-	.	gene_id "CUFF.5"; transcript_id "CUFF.5.2"; exon_number "2"; FPKM "17.7850894771"; frac "0.322877"; conf_lo "16.994111"; conf_hi "18.576068"; cov "93.732896";
tp.assembly.567468735.1	Cufflinks	exon	25287	28599	124	-	.	gene_id "CUFF.5"; transcript_id "CUFF.5.2"; exon_number "3"; FPKM "17.7850894771"; frac "0.322877"; conf_lo "16.994111"; conf_hi "18.576068"; cov "93.732896";

OUTPUT example (--export_mode=model):

tp.assembly.567468735.1 Cufflinks       gene    20970   21966   .       -       .       ID=CUFF.4
tp.assembly.567468735.1 Cufflinks       mRNA    20970   21966   .       -       .       ID=CUFF.4.1;Parent=CUFF.4
tp.assembly.567468735.1 Cufflinks       CDS     20970   21966   .       -       0       ID=CUFF.4.1.CDS.1;Parent=CUFF.4.1
tp.assembly.567468735.1 Cufflinks       exon    20970   21966   .       -       .       ID=CUFF.4.1.exon.1;Parent=CUFF.4.1
tp.assembly.567468735.1 Cufflinks       gene    22213   23755   .       -       .       ID=CUFF.5
tp.assembly.567468735.1 Cufflinks       mRNA    22213   23755   .       -       .       ID=CUFF.5.1;Parent=CUFF.5
tp.assembly.567468735.1 Cufflinks       CDS     22213   22259   .       -       0       ID=CUFF.5.1.CDS.1;Parent=CUFF.5.1
tp.assembly.567468735.1 Cufflinks       CDS     22294   23755   .       -       1       ID=CUFF.5.1.CDS.2;Parent=CUFF.5.1
tp.assembly.567468735.1 Cufflinks       exon    22213   22259   .       -       .       ID=CUFF.5.1.exon.1;Parent=CUFF.5.1
tp.assembly.567468735.1 Cufflinks       exon    22294   23755   .       -       .       ID=CUFF.5.1.exon.2;Parent=CUFF.5.1
tp.assembly.567468735.1 Cufflinks       mRNA    22213   28599   .       -       .       ID=CUFF.5.2;Parent=CUFF.5
tp.assembly.567468735.1 Cufflinks       CDS     22213   23755   .       -       0       ID=CUFF.5.2.CDS.1;Parent=CUFF.5.2
tp.assembly.567468735.1 Cufflinks       CDS     23887   25044   .       -       2       ID=CUFF.5.2.CDS.2;Parent=CUFF.5.2
tp.assembly.567468735.1 Cufflinks       CDS     25287   28599   .       -       2       ID=CUFF.5.2.CDS.3;Parent=CUFF.5.2
tp.assembly.567468735.1 Cufflinks       exon    22213   23755   .       -       .       ID=CUFF.5.2.exon.1;Parent=CUFF.5.2
tp.assembly.567468735.1 Cufflinks       exon    23887   25044   .       -       .       ID=CUFF.5.2.exon.2;Parent=CUFF.5.2
tp.assembly.567468735.1 Cufflinks       exon    25287   28599   .       -       .       ID=CUFF.5.2.exon.3;Parent=CUFF.5.2

Caveats:
- Notice that the input defined above describes multiple isoforms for the same gene.  Some downstream
  GFF3 parsers may not support this (though they should.)
- This script currently skips all attributes in the last column of the input GTF except for
  the gene_id and transcript_id.

OUTPUT example (--export_mode=cDNA_match):

ChromosomeIII_BmicrotiR1        Cufflinks       cDNA_match      4       3995    .       -       .       ID=CUFF.1.1
ChromosomeIII_BmicrotiR1        Cufflinks       cDNA_match      5052    6630    .       -       .       ID=CUFF.2.1
ChromosomeIII_BmicrotiR1        Cufflinks       cDNA_match      7857    10480   .       -       .       ID=CUFF.3.1
ChromosomeIII_BmicrotiR1        Cufflinks       cDNA_match      10101   10426   .       +       .       ID=CUFF.4.1
ChromosomeIII_BmicrotiR1        Cufflinks       cDNA_match      10447   10752   .       +       .       ID=CUFF.4.1
ChromosomeIII_BmicrotiR1        Cufflinks       cDNA_match      10775   10894   .       +       .       ID=CUFF.4.1
ChromosomeIII_BmicrotiR1        Cufflinks       cDNA_match      10916   11185   .       +       .       ID=CUFF.4.1
ChromosomeIII_BmicrotiR1        Cufflinks       cDNA_match      11207   11275   .       +       .       ID=CUFF.4.1

This type of output is required for tools like EVM.  It has only cDNA_match features which are
segmented.  Each segment shares the same ID but has no explicit Parent groupings.

Author: Joshua Orvis (jorvis AT gmail)
"""

import argparse
import sys
from collections import defaultdict

from biocode import gff, things


def main():
    parser = argparse.ArgumentParser( description='A GTF -> GFF3 conversion script for Cufflinks output')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input GTF file' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output GFF file to be created' )
    parser.add_argument('-e', '--export_mode', type=str, required=False, default='model', help='Export mode for results (model or cDNA_match)' )
    args = parser.parse_args()

    if args.export_mode not in ['model', 'cDNA_match']:
        raise Exception("ERROR: the only valid values for --export_mode are 'model' or 'cDNA_match'")
    
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

    # each gb_record is a SeqRecord object
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
        
        if ftype == 'transcript':
            if args.export_mode == 'model':
                if current_gene is not None and current_gene.id != gene_id:
                    gene.print_as(fh=ofh, source='Cufflinks', format='gff3')

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

            elif args.export_mode == 'cDNA_match':
                if current_match is not None and current_match.id != transcript_id:
                    match.print_as( fh=ofh, source='Cufflinks', format='gff3' )
                
                match = things.Match(id=transcript_id, subclass='cDNA_match', length=fmax - fmin)
                match.locate_on( target=current_assembly, fmin=fmin, fmax=fmax, strand=strand )
                current_match = match
            
        elif ftype == 'exon':
            exon_number = gff.column_9_value(col9, 'exon_number').replace('"', '')
            
            if args.export_mode == 'model':
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
                
            elif args.export_mode == 'cDNA_match':
                mp_id = "{0}.match_part.{1}".format(transcript_id, exon_number)
                mp = things.MatchPart(id=mp_id, parent=current_match, length=fmax - fmin)
                mp.locate_on( target=current_assembly, fmin=fmin, fmax=fmax, strand=strand )
                current_match.add_part(mp)

    # don't forget to do the last gene, if there were any
    if args.export_mode == 'model':
        if current_gene is not None:
            gene.print_as(fh=ofh, source='GenBank', format='gff3')
            
    elif args.export_mode == 'cDNA_match':
        if current_match is not None:
            match.print_as( fh=ofh, source='Cufflinks', format='gff3' )


if __name__ == '__main__':
    main()







