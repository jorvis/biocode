#!/usr/bin/env python3.2

"""
INPUT:

The input is a tab-delimited output file that is the product of parsing eukaryotic
annotation evidence.  It has the following columns:

   0. Transcript ID
   1. ORF ID
   2. ORF length
   3. BLAST description
   4. BLAST E-value
   5. HMM description
   6. HMM total score
   7. HMM total E-value
   8. GO terms (derived from HMM), comma-seprated
   9. EC numbers (dervied from HMM), comma-separated

OUTPUT:

WARNING: This script creates a lot of files in the specified output directory.
To be more specific, it will create (N*2)+1 where N is the number of molecular
fragments in the input nucleotide set.  This is true whether that number is
10 supercontigs of a large genome or 50,000 transcripts annotated without
a genome.

Pathologic format described:
http://bioinformatics.ai.sri.com/ptools/pathologic-format.pdf

Example (genetic elements list):
http://bioinformatics.ai.sri.com/ptools/sample-genetic-elements.dat

Example (per molecule):
http://bioinformatics.ai.sri.com/ptools/tpal.pf

"""
import argparse
import os
import biocodeutils
import biothings

def main():
    parser = argparse.ArgumentParser( description='Transforms a tab-delimited annotation file to PathoLogic format')

    ## output file to be written
    parser.add_argument('-a', '--annotation_tab', type=str, required=True, help='Path to an input file to be parsed' )
    parser.add_argument('-g', '--genomic_fasta', type=str, required=True, help='Underlying nucleotide FASTA file for the annotated proteins' )
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    ntseqs = biocodeutils.fasta_dict_from_file( args.genomic_fasta )

    write_elements_file(ntseqs, args.output_dir )
    write_seq_files(ntseqs, args.output_dir)

    for line in open(args.annotation_tab):
        if line.startswith("#"):
            continue


def write_seq_files( seqs, outdir ):
    for seq in seqs:
        ## any changes to this path convention need to be also updated in write_elements_file()
        seqout_path = "{0}/{1}.fna".format(outdir, seq)
        seqout = open(seqout_path, mode='wt', encoding='utf-8')
        seqout.write( ">{0} {1}\n".format(seq, seqs[seq]['h']) )

        for i in range(0, len(seqs[seq]['s']), 60):
            seqout.write(seqs[seq]['s'][i : i + 60] + "\n")
        
        seqout.close()

        
def write_elements_file( seqs, outdir ):
    elements_file_path = "{0}/genetic-elements.dat".format(outdir)

    fout = open(elements_file_path, mode='wt', encoding='utf-8')

    for seq in seqs:
        fout.write( "ID\t{0}\n".format(seq) )

        if not seqs[seq]['h'].isspace():
            fout.write( "NAME\t{0} {1}\n".format(seq, seqs[seq]['h']) )

        fout.write( "CIRCULAR?\tN\n" )
        fout.write( "ANNOT-FILE\t{0}/{1}.pf\n".format(outdir, seq) ) ## must be full path
        fout.write( "SEQ-FILE\t{0}/{1}.fna\n".format(outdir, seq) ) ## must be full path
        fout.write("//\n")
    
    fout.close()

if __name__ == '__main__':
    main()







