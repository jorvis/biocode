#!/usr/bin/env python3

"""
This script can be used to transform a GFF3 file into a TBL file suitable for submitting
to NCBI.  Its encoding is meant to include prokaryotes and eukaryotes, though it was
written/tested first for eukaryotes.

The options:

--input_file:  Must follow the GFF3 specification

--output_base:  Will be used to write NCBI table and fasta formatted files.  TBL defined here:
                http://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation

--lab_name:  A short, unique identifier, no spaces, which will form part of the exported
             protein and transcript IDs.  See the documentation for examples:
             Documentation: http://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation#protein_id

--ncbi_acc_prefix:  Assigned by NCBI, this is the accession prefix for this submission.

--fasta:  NCBI has strict requirements for feature identifiers, and this script will reformat them
          unless you have already prepared them correctly.  If corrections need to be made, a new
          corresponding FASTA file will also be written, which requires this input source fasta file.
                
--go_obo:  Optional.  If your GFF3 annotation includes GO attributes, this is required because
           the GFF3 has only the ID but the TBL file also requires the descriptor and class, so
           this enables a lookup of these using the annotated GO ids.

"""

import argparse
import biocodegff
import biocodetbl
import biothings
import biocodeutils
import os


def main():
    parser = argparse.ArgumentParser( description='Create a TBL file for submission to NCBI from GFF3')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_base', type=str, required=True, help='Base name of output files to be created' )
    parser.add_argument('-ln', '--lab_name', type=str, required=True, help='Required by NCBI to identify the submitting group' )
    parser.add_argument('-nap', '--ncbi_acc_prefix', type=str, required=True, help='Required and assigned by NCBI' )
    parser.add_argument('-gf', '--genomic_fasta', type=str, required=False, help='FASTA file of genomic sequence, if not embedded in GFF' )
    parser.add_argument('-go', '--go_obo', type=str, required=False, help='GO terms will not be exported unless you pass the path to a GO OBO file')
    
    args = parser.parse_args()

    (assemblies, features) = biocodegff.get_gff3_features( args.input_file )
    
    if args.genomic_fasta is not None:
        biocodeutils.add_assembly_fasta(assemblies, args.genomic_fasta)
        
    new_assemblies = dict() 

    ## We need to first check the ID format
    reformat_IDs = True

    ## maps old IDs (like tp.assembly.567468735.1) to new ones (like AAGK01000001)
    asm_id_map = dict()
    asm_num = 1

    for asm_id in assemblies:
        # pre-formatted IDs are like this: gnl|WGS:XXXX|SeqID|gb|XXXX01xxxxxx
        if asm_id.startswith('gnl|WGS:'):
            reformat_IDs = False
            break
        else:
            new_id = "gnl|WGS:{0}|SeqID|gb|{0}01{1:06d}".format(args.ncbi_acc_prefix, asm_num)
            asm_id_map[asm_id] = new_id
            asm_num += 1
            new_assemblies[new_id] = assemblies[asm_id]
            new_assemblies[new_id].id = new_id

    if reformat_IDs == True:
        assemblies = new_assemblies

    # >gi|68352484|gb|AAGK01000001.1|
    # AAGK01000001	NC_007344.1	tp.assembly.567468735.1

    ofh = open("{0}.tbl".format(args.output_base), 'wt')
    biocodetbl.print_tbl_from_assemblies(assemblies=assemblies, ofh=ofh, go_obo=args.go_obo, lab_name=args.lab_name)

    mset = biothings.AssemblySet()
    mset.load_from_dict(assemblies)
    mset.write_fasta(path="{0}.fna".format(args.output_base))

if __name__ == '__main__':
    main()







