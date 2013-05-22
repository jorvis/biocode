#!/usr/bin/env python3.2

"""
This evolved into a very specialized script that I'm only saving to tease apart
like Lego bits to make something more permanent later.

INPUT:

Three total inputs are required.  

The first is a tab-delimited output file (-a) that is the product of parsing eukaryotic
annotation evidence.  It has the following columns, the first of which must strictly
match this naming convention.

   0. Transcript ID
   1. ORF ID
   2. ORF length
   3. BLAST description
   4. BLAST E-value
   5. HMM description
   6. HMM total score
   7. HMM total E-value
   8. GO terms (derived from HMM), comma-separated
   9. EC numbers (dervied from HMM), comma-separated

Example:

comp104220_c0_seq1      m.14260 1160    ATP-dependent zinc metalloprotease FtsH OS=Shigella flexneri GN=ftsH PE=3 SV=1  0       ATP-dependent metallopeptidase HflB     613.7   2.30E-184       "GO:0004222,GO:0006508,GO:0016887,GO:0051301"   3.4.24.-

The second, defined by -g, is the genomic FASTA file.  Usually genomic, this could even
be assembled RNA-seq transcripts, with headers like:

   >comp20_c0_seq1 len=234 path=[1:0-233]

Finally, the protein (-p) file is required that is currently the output from transdecoder:

>m.13585 g.13585  ORF g.13585 m.13585 type:3prime_partial len:76 (+) comp100033_c0_seq1:118-348(+)
MGVGILTVGCGGQLAVNNPPPSLPLGILREITCIFSHLSMNSGQTRLRSSFLHTGVNTGY
FLDTGVDTGHNLDTGV


Development testing command:
$ /opt/python-3.3.1/bin/python3.3 ~/svn/biocode/sandbox/jorvis/convert_euk_autoannotate_to_pathologic.py -a best_candidates.eclipsed_orfs_removed.evidence_annotation_EC_transcripts.tab -g combined_Trinity.fasta -o pathologic_out -p best_candidates.eclipsed_orfs_removed.pep


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
import re
import biocodeutils
import biothings
import bioannotation

def main():
    parser = argparse.ArgumentParser( description='Transforms a tab-delimited annotation file to PathoLogic format')

    ## output file to be written
    parser.add_argument('-a', '--annotation_tab', type=str, required=True, help='Path to an input file to be parsed' )
    parser.add_argument('-g', '--genomic_fasta', type=str, required=True, help='Underlying nucleotide FASTA file for the annotated proteins' )
    parser.add_argument('-p', '--protein_fasta', type=str, required=True, help='Protein input sequences to the pipeline, with specific headers required' )
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    molecules = biocodeutils.fasta_dict_from_file( args.genomic_fasta )

    #write_elements_file(molecules, args.output_dir )
    #write_seq_files(molecules, args.output_dir)

    genes = dict()
    
    for line in open(args.annotation_tab):
        if line.startswith("#"):
            continue

        parse_annotation_line( line, genes, molecules )
        
    write_annotation_files( genes, molecules, args.output_dir )



def write_annotation_files( genes, molecules, outdir ):
    for gene_id in genes:
        gene = genes[gene_id]

        for mRNA in gene.mRNAs():
            ## any changes to this path convention need to be also updated in write_elements_file()
            annot_path = "{0}/{1}.pf".format(outdir, mRNA.id)
            annotout = open(annot_path, mode='wt', encoding='utf-8')

            for CDS in mRNA.CDSs():
                
                annotout.write("ID\t{0}\n".format(CDS.id))
                annotout.write("PRODUCT-TYPE\tP\n")

                if CDS.annotation.product_name is None:
                    annotout.write("FUNCTION\tHypothetical protein\n")
                else:
                    annotout.write("FUNCTION\t{0}\n".format(CDS.annotation.product_name))
                    
                annotout.write("//\n")

            annotout.close()


def get_gene_id_from_transcript( transcript_id ):
    m = re.match( "(comp\d+_c\d+)_seq\d+", transcript_id )
    if m:
        return m.group(1)
    else:
        raise Exception("ERROR: expected transcript ID to be like 'comp20_c0_seq1' but got: {0}".format(transcript_id) )
        

def parse_annotation_line(line, genes, molecules):
    cols = line.split("\t")

    transcript_id = cols[0]
    CDS_id = cols[1]
    gene_id = get_gene_id_from_transcript( transcript_id )
    
    if cols[5] is None:
        gene_product_name = cols[3]
    else:
        gene_product_name = cols[5]

    if transcript_id not in molecules:
        raise Exception("ERROR: found molecule {0} in referenced in annotation tab file but not in genomic_fasta file".format(transcript_id))

    if gene_id in genes:
        gene = genes[gene_id]
    else:
        gene = biothings.Gene( id=gene_id )
        genes[gene_id] = gene

    mRNA = biothings.mRNA( id=transcript_id )
    gene.add_mRNA( mRNA )

    annotation = bioannotation.FunctionalAnnotation( product_name=gene_product_name )
    
    CDS = biothings.CDS( id=CDS_id, annotation=annotation )
    mRNA.add_CDS( CDS )

    

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







