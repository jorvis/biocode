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

"""

import argparse
import biocodegenbank
import biocodegff
import os
import sys

from jinja2 import Environment, FileSystemLoader

PATH = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_ENVIRONMENT = Environment(
    autoescape=False,
    loader=FileSystemLoader(os.path.join(PATH, '../data')),
    trim_blocks=False)

def main():
    parser = argparse.ArgumentParser( description='Converts GFF3 into a GenBank flat file')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input GFF3 file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to a Genbank flat file to be created' )
    parser.add_argument('-mt', '--molecule_type', type=str, required=False, default='DNA', help='Molecule type' )
    parser.add_argument('-gbd', '--genbank_division', type=str, required=False, default='.', help='GenBank Division (3-letter abbreviation)' )
    parser.add_argument('-md', '--modification_date', type=str, required=False, default='DD-MMM-YYYY', help='The modification date for header in format like 21-JUN-1999' )
    parser.add_argument('-org', '--organism', type=str, required=False, default='.', help='Full organism name (including strain)' )
    parser.add_argument('-str', '--strain', type=str, required=False, help="Only the strain designation, which is written to the FEATURES.source element" )
    parser.add_argument('-d', '--definition', type=str, required=False, default='.', help='Brief description of sequence; includes information such as source organism, gene name/protein name, or some description of the sequence\'s function.' )
    parser.add_argument('-s', '--source', type=str, required=False, default='.', help='Free-format information including an abbreviated form of the organism name, sometimes followed by a molecule type.' )
    parser.add_argument('-t', '--taxon_id', type=int, required=False, help='NCBI taxon ID, if known' )
    parser.add_argument('-l', '--lineage', type=str, required=False, default='Unknown', help='Semicolon-delimited lineage of the organism e.g., "Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Piroplasmida; Theileriidae; Theileria"' )
    args = parser.parse_args()

    # line-wrap lineage to stay below 79 character GenBank flat file width
    lineage = biocodegenbank.line_wrap_lineage_string(args.lineage)

    (assemblies, features) = biocodegff.get_gff3_features( args.input_file )
    ofh = sys.stdout
    if args.output_file is not None:
        ofh = open(args.output_file, 'wt')

    for assembly_id in assemblies:
        assembly = assemblies[assembly_id]

        context = { 'locus':assembly_id, 'molecule_size':assembly.length, 'molecule_type':args.molecule_type,
                    'division':args.genbank_division, 'modification_date':args.modification_date,
                    'accession':'.', 'version':'.', 'gi':'.',
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
            biocodegenbank.print_biogene( gene=gene, fh=ofh, on=assembly )

        ofh.write("//\n")


if __name__ == '__main__':
    main()







