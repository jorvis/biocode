#!/usr/bin/env python3

"""
Reads one or more evidence files and reports coverage information for user-chosen
feature types in a GFF3 file.  This allows you to answer, as an example, the following
question, "What is the coverage of all mRNAs in my source GFF3 file given a set
of reference annotations (gff3), aligned reads (bam) and probe locations (bed)?"

INPUT

Primary input type:

   gff3 - Source features upon which the output will be reported.
          This should be defined like FILE:TYPE, such as:
          /path/to/yourfile.gff3:mRNA

Supported evidence types (type relies on matching the file extension):

   gff, gff3 - Such as a reference annotation to compare
   bed - Such as aligned probe locations
   bam - Such as RNA-Seq read alignments

OUTPUT

Put output example here
   
Contact: jorvis@gmail.com
"""

import argparse
import biocodegff
import biocodeutils
import codecs
import os
import sys

def main():
    parser = argparse.ArgumentParser( description='Provides coverage information for features in a GFF3 file')

    ## output file to be written
    parser.add_argument('evidence_files', metavar='N', type=str, nargs='+', help='Path to one or more evidence files, separated by spaces' )
    parser.add_argument('-r', '--reference', type=str, required=True, help='Input path to the reference GFF3 file. So we know what feature type to report on, format should be like FILE:TYPE' )
    parser.add_argument('-f', '--fasta', type=str, required=True, help='Input path to the reference FASTA file.' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Optional path to an output file to be created, else prints on STDOUT' )
    args = parser.parse_args()

    ## parse the fasta
    fasta = biocodeutils.fasta_dict_from_file(args.fasta)

    ## initialize the coverage array for each file (and do a sanity check so we don't waste the user's time)
    cov = dict()  # structure here is like cov['ev_file'][pos:coverage]
    allowed_extensions = ['gff3', 'bed', 'bam']
    for ev_file in args.evidence_files:
        
        for ext in allowed_extensions:

    ## open the output file
    fout = None
    if args.output_file is None:
        fout = codecs.getwriter('utf8')(sys.stdout.buffer)
    else:
        fout = open(args.output_file, "w")

    ## The input file should be defined as $path:$feattype
    if ':' not in args.reference:
        raise Exception("ERROR: input_file must be like /path/to/some.gff3:mRNA")
        
    ref_file_parts = args.reference.split(':')
    print("DEBUG: part count: {0}".format(len(ref_file_parts)))
        
    if ref_file_parts[0].endswith('.gff3'):
        (ref_assemblies, ref_features) = biocodegff.get_gff3_features( ref_file_parts[0] )
    else:
        raise Exception("ERROR: Expected input file (-i) to have a gff3 extension, got {0}".format(ref_file_parts[0]))

    for ev_file in args.evidence_files:
        if ev_file.endswith('bed'):
            pass
        else:
            print("INFO: ignoring evidence file {0} because it doesn't have a recognized file extension")
        




if __name__ == '__main__':
    main()







