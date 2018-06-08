#!/usr/bin/env python3

"""
Converts GFF3 representing gene models to Genbank flat-file format.

GFF3 specification:
http://www.sequenceontology.org/gff3.shtml

Genbank flat file specification:
https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html

--molecule_type:
    http://www.ncbi.nlm.nih.gov/Sequin/sequin.hlp.html#Molecule

--genbank_division:
    http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html#GenBankDivisionB

MAJOR assumptions:
- Functional annotation is attached to a polypeptide feature

"""

import argparse
import os
import pickle
import sys

from biocode import utils, genbank, gff
from jinja2 import Environment, FileSystemLoader

from pkg_resources import Requirement, resource_filename

# If biocode is installed via a GitHub checkout, the path to the template will be referential
template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/genbank_flat_file_header.template')

# if referential discovery didn't work, this was probably a PyPi install.  Use pkg_resources to find it instead.
if os.path.isfile(template_path):
    template_dir = os.path.dirname(template_path)
else:
    template_path = resource_filename(Requirement.parse("biocode"), "genbank_flat_file_header.template")
    template_dir = "{0}/{1}".format(os.path.dirname(template_path), '/biocode/data')

# for now, these are the same
data_dir = template_dir
    
TEMPLATE_ENVIRONMENT = Environment(
    autoescape=False,
    loader=FileSystemLoader(template_dir),
    trim_blocks=False)

def main():
    parser = argparse.ArgumentParser( description='Converts GFF3 into a GenBank flat file')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input GFF3 file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to a Genbank flat file to be created. Supersedes --output_dir if both are specified.' )
    parser.add_argument('-od', '--output_dir', type=str, required=False, help='Path to an output directory. If this option is specified then each input assembly will be written to a separate GenBank output file, named with the assembly_id.' )
    parser.add_argument('-g', '--genome_fasta', type=str, required=False, help='Optional.  You must specify this unless the FASTA sequences for the molecules are embedded in the GFF')
    parser.add_argument('-go', '--go_index', type=str, required=False, default='', help='Pickled GO index (created by biocode/general/make_go_index.py). By default, reads from within the data directory within the biocode distribution')
    parser.add_argument('-mt', '--molecule_type', type=str, required=False, default='DNA', help='Molecule type' )
    parser.add_argument('-gbd', '--genbank_division', type=str, required=False, default='.', help='GenBank Division (3-letter abbreviation)' )
    parser.add_argument('-md', '--modification_date', type=str, required=False, default='DD-MMM-YYYY', help='The modification date for header in format like 21-JUN-1999' )
    parser.add_argument('-org', '--organism', type=str, required=False, default='.', help='Full organism name (including strain)' )
    parser.add_argument('-str', '--strain', type=str, required=False, help="Only the strain designation, which is written to the FEATURES.source element" )
    parser.add_argument('-d', '--definition', type=str, required=False, default='.', help='Brief description of sequence; includes information such as source organism, gene name/protein name, or some description of the sequence\'s function.' )
    parser.add_argument('-s', '--source', type=str, required=False, default='.', help='Free-format information including an abbreviated form of the organism name, sometimes followed by a molecule type.' )
    parser.add_argument('-t', '--taxon_id', type=int, required=False, help='NCBI taxon ID, if known' )
    parser.add_argument('-l', '--lineage', type=str, required=False, default='Unknown', help='Semicolon-delimited lineage of the organism e.g., "Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Piroplasmida; Theileriidae; Theileria"' )
    parser.add_argument('-seq', '--include_sequence', action='store_true', help='Include sequence (if present) in the output GenBank flat file(s).' )
    parser.add_argument('-p', '--locus_id_prefix', required=False, default='', help='Prefix to add to the GenBank LOCUS id in the output GenBank flat file(s).' )
    args = parser.parse_args()

    # check that output directory exists
    if args.output_dir is not None:
        if not os.path.isdir(args.output_dir):
            sys.stderr.write("FATAL: the specified output directory (" + args.output_dir + ") does not exist\n");
            exit(1)

    # line-wrap lineage to stay below 79 character GenBank flat file width
    lineage = genbank.line_wrap_lineage_string(args.lineage)

    if args.go_index == '':
        go_index_path = "{0}/go.pickle".format(data_dir)
    else:
        go_index_path = args.go_index

    if os.path.isfile(go_index_path):
        go_index = pickle.load(open(go_index_path, 'rb'))
    else:
        raise Exception("ERROR: Expected to find a pickled GO index at the following path: {0}".format(go_index_path))

    (assemblies, features) = gff.get_gff3_features(args.input_file)
    ofh = sys.stdout
    if args.output_file is not None:
        if args.output_dir is None:
            ofh = open(args.output_file, 'wt')
        else:
            sys.stderr.write("WARN: both -o/--output_file and -od/--output_dir were passed so the former will be ignored\n")

    # deal with the FASTA file if the user passed one
    if args.genome_fasta is not None:
        process_assembly_fasta(assemblies, args.genome_fasta)

    for assembly_id in assemblies:
        locus_id = args.locus_id_prefix + assembly_id
        if args.output_dir is not None:
            ofn = args.output_dir + "/" + locus_id + ".gbk"
            ofh = open(ofn, 'wt')
        assembly = assemblies[assembly_id]

        context = { 'locus':locus_id, 'molecule_size':assembly.length, 'molecule_type':args.molecule_type,
                    'division':args.genbank_division, 'modification_date':args.modification_date,
                    'accession':'.', 'version':'.', 
                    'source':args.source, 'definition':args.definition, 'organism':args.organism,
                    'lineage':lineage
        }
        
        header = TEMPLATE_ENVIRONMENT.get_template('genbank_flat_file_header.template').render(context)
        ofh.write(header)
        ofh.write("\nFEATURES             Location/Qualifiers\n")
        ofh.write("     source          1..{0}\n".format(assembly.length))
        ofh.write("                     /organism=\"{0}\"\n".format(args.organism))
        ofh.write("                     /mol_type=\"genomic DNA\"\n")

        if args.strain is not None:
            ofh.write("                     /strain=\"{0}\"\n".format(args.strain))

        if args.taxon_id is not None:
            ofh.write("                     /db_xref=\"taxon:{0}\"\n".format(args.taxon_id))
        
        for gene in assemblies[assembly_id].genes():
            genbank.print_biogene(gene=gene, fh=ofh, on=assembly, go_index=go_index)

        if args.include_sequence:
            ofh.write("ORIGIN\n")
            genbank.print_sequence(seq=assembly.residues, fh=ofh)

        ofh.write("//\n")
        # there may be multiple output files
        if args.output_dir is not None:
            ofh.close()

    # there is only one output file
    if args.output_dir is None:
        ofh.close()

def process_assembly_fasta(mols, fasta_file):
    fasta_seqs = utils.fasta_dict_from_file(fasta_file)

    for mol_id in mols:
        # check if the FASTA file provides sequence for this
        if mol_id in fasta_seqs:
            mol = mols[mol_id]
            mol.residues = fasta_seqs[mol_id]['s']
            mol.length   = len(mol.residues)
        

if __name__ == '__main__':
    main()







