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
   9. EC numbers (derived from HMM), comma-separated

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

From EUK-337: When completed, please write results here:
/local/projects-t2/e_affinis/rna_seq/analysis/combined_assembly/pathway_tools

OUTPUT:

WARNING: This script creates a lot of files in the specified output directory.
To be more specific, it will create (N*2)+1 where N is the number of molecular
fragments in the input nucleotide set.  This is true whether that number is
10 supercontigs of a large genome or 50,000 transcripts annotated without
a genome.

I try to limit the number of files showing up in the base directory by creating
subdirectories for the files in a repeatable way.  How?  A bit of a hack, but I
take the md5sum of each molecules name and then only use the first two characters
of it.

Pathologic format described:
http://bioinformatics.ai.sri.com/ptools/pathologic-format.pdf

Example (genetic elements list):
http://bioinformatics.ai.sri.com/ptools/sample-genetic-elements.dat

Example (per molecule):
http://bioinformatics.ai.sri.com/ptools/tpal.pf


"""
import argparse
import hashlib
import os
import re

from biocode import utils, things


def main():
    parser = argparse.ArgumentParser( description='Transforms a tab-delimited annotation file to PathoLogic format')

    ## output file to be written
    parser.add_argument('-a', '--annotation_tab', type=str, required=True, help='Path to an input file to be parsed' )
    parser.add_argument('-g', '--genomic_fasta', type=str, required=True, help='Underlying nucleotide FASTA file for the annotated proteins' )
    parser.add_argument('-p', '--protein_fasta', type=str, required=True, help='Protein input sequences to the pipeline, with specific headers required' )
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    molecules = utils.fasta_dict_from_file(args.genomic_fasta)
    protein_coords = get_protein_coordinates_from_FASTA(args.protein_fasta)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    create_subdirectories(molecules, args.output_dir)
    write_elements_file(molecules, args.output_dir)
    write_seq_files(molecules, args.output_dir)

    genes = dict()
    
    for line in open(args.annotation_tab):
        if line.startswith("#"):
            continue

        parse_annotation_line( line, genes, molecules )
        
    write_annotation_files( genes, molecules, protein_coords, args.output_dir )




def get_base_directory_name( mol_name ):
    m = hashlib.md5()
    m.update(mol_name.encode())
    digest = m.hexdigest()
    return digest[:2]

def create_subdirectories( mols, output_dir ):
    for mol in mols:
        base_name = get_base_directory_name( mol )
        subdir_path = "{0}/{1}".format(output_dir, base_name)

        if not os.path.exists(subdir_path):
            os.makedirs(subdir_path)


def get_protein_coordinates_from_FASTA(protein_fasta):
    '''
    This function is probably only relevant to a limited number of tasks where the protein input
    FASTA file to the pipeline was produced by transdecoder, which incorporates the predicted
    ORF coordinates into the FASTA header, like this:

    >m.13585 g.13585  ORF g.13585 m.13585 type:3prime_partial len:76 (+) comp100033_c0_seq1:118-348(+)

    Of all that, the only part we care about is that the first model number 'm.13585' matches that
    of the second column in the annotation file, along with the matching genomic molecule name at
    the end of the header 'comp100033_c0_seq1:118-348(+)

    Returns a dict keyed on the model name (like 'm.13585') with the followed keyed values:
       'mol'    = molecule ID (like 'comp100033_c0_seq1')
       'fmin'   = 0-interbase start coordinate (117, from example above)
       'fmax'   = 0-interbase stop coordinate (348, from example above)
       'strand' = 1, 0 or -1 direction
    '''
    protein_locs = dict()

    fasta_dict = utils.fasta_dict_from_file(protein_fasta)

    pattern = re.compile('(comp\d+_c\d+_seq\d+)\:(\d+)\-(\d+)\(\+\)')

    for model_id in fasta_dict:
        if model_id in protein_locs:
            raise Exception("ERROR: found duplicate model ID in file: {0}".format(protein_fasta) )

        m = pattern.search(fasta_dict[model_id]['h'])

        if m:
            protein_locs[model_id] = { 'mol': m.group(1), 'fmin': int(m.group(2)) - 1, 'fmax': int(m.group(3)), 'strand': 1}
        else:
            raise Exception("ERROR: unexpected header format.  Expected to parse something like comp100033_c0_seq1:118-348(+).  Got: {0}".format(fasta_dict[model_id]['h']))

    return protein_locs
    



def write_annotation_files( genes, molecules, locs, outdir ):
    for gene_id in genes:
        gene = genes[gene_id]

        for mRNA in gene.mRNAs():
            ## any changes to this path convention need to be also updated in write_elements_file()
            base_dir = get_base_directory_name(mRNA.id)
            annot_path = "{0}/{1}/{2}.pf".format(outdir, base_dir, mRNA.id)
            annotout = open(annot_path, mode='wt', encoding='utf-8')

            for CDS in mRNA.CDSs():
                
                annotout.write("ID\t{0}\n".format(CDS.id))

                if CDS.id in locs:
                    annotout.write("STARTBASE\t{0}\n".format(str(locs[CDS.id]['fmin'] + 1)) )
                    annotout.write("ENDBASE\t{0}\n".format(str(locs[CDS.id]['fmax'])) )

                if len(CDS.annotation.ec_numbers):
                    for ec_annot in CDS.annotation.ec_numbers:
                        annotout.write("EC\t{0}\n".format(ec_annot.number))

                if len(CDS.annotation.go_annotations):
                    for go_annot in CDS.annotation.go_annotations:
                        annotout.write("DBLINK\tGO:{0}\n".format(go_annot.go_id))

                if CDS.annotation.product_name is None:
                    annotout.write("FUNCTION\tHypothetical protein\n")
                else:
                    annotout.write("FUNCTION\t{0}\n".format(CDS.annotation.product_name))


                annotout.write("PRODUCT-TYPE\tP\n")
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

    if len(cols) != 10:
        print("WARNING: Ignoring the following line because I expected 10 columns:\n{0}".format(line))
        return False
    
    cols[9] = cols[9].rstrip()

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
        gene = things.Gene(id=gene_id)
        genes[gene_id] = gene

    mRNA = things.mRNA(id=transcript_id)
    gene.add_mRNA( mRNA )

    annotation = annotation.FunctionalAnnotation(product_name=gene_product_name)

    ec_num_pattern = re.compile('\d+.')

    if cols[9] is not None:
        ec_nums = cols[9].split(',')

        for ec_num in ec_nums:
            m = ec_num_pattern.search(ec_num)

            if m:
                ec = annotation.ECAnnotation(number=ec_num)
                annotation.add_ec_number( ec )

    go_pattern = re.compile('(\d+)')
    if cols[8] is not None:
        go_terms = cols[8].split(',')

        for go_term in go_terms:
            m = go_pattern.search(go_term)

            if m:
                go = annotation.GOAnnotation(go_id=go_term)
                annotation.add_go_annotation( go )
                
    CDS = things.CDS(id=CDS_id, annotation=annotation)
    mRNA.add_CDS( CDS )

    

def write_seq_files( seqs, outdir ):
    for seq in seqs:
        ## any changes to this path convention need to be also updated in write_elements_file()
        base_dir = get_base_directory_name(seq)
        seqout_path = "{0}/{1}/{2}.fna".format(outdir, base_dir, seq)
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
        base_dir = get_base_directory_name(seq)

        if not seqs[seq]['h'].isspace():
            fout.write( "NAME\t{0} {1}\n".format(seq, seqs[seq]['h']) )

        fout.write( "CIRCULAR?\tN\n" )
        fout.write( "ANNOT-FILE\t{0}/{1}/{2}.pf\n".format(outdir, base_dir, seq) ) ## must be full path
        fout.write( "SEQ-FILE\t{0}/{1}/{2}.fna\n".format(outdir, base_dir, seq) ) ## must be full path
        fout.write("//\n")
    
    fout.close()

if __name__ == '__main__':
    main()







