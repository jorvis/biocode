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
   sam - Such as RNA-Seq read alignments
   pileup - As created by "samtools pileup"

OUTPUT

Put output example here
   
Contact: jorvis@gmail.com
"""

import argparse
import codecs
import sys

from biocode import utils, gff


def main():
    parser = argparse.ArgumentParser( description='Provides coverage information for features in a GFF3 file')

    ## output file to be written
    parser.add_argument('evidence_files', metavar='N', type=str, nargs='+', help='Path to one or more evidence files, separated by spaces' )
    parser.add_argument('-r', '--reference', type=str, required=True, help='Input path to the reference GFF3 file. So we know what feature type to report on, format should be like FILE:TYPE' )
    parser.add_argument('-f', '--fasta', type=str, required=True, help='Input path to the reference FASTA file.' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Optional path to an output file to be created, else prints on STDOUT' )
    args = parser.parse_args()

    ## parse the fasta
    fasta = utils.fasta_dict_from_file(args.fasta)

    ## open the output file
    fout = None
    if args.output_file is None:
        fout = codecs.getwriter('utf8')(sys.stdout.buffer)
    else:
        fout = open(args.output_file, "w")

    ####################################################
    ## Sanity checks

    allowed_extensions = ['bed', 'gff3', 'pileup', 'sam']
    for ev_file in args.evidence_files:
        valid_ext_found = False
        
        for ext in allowed_extensions:
            if ev_file.endswith(ext):
                valid_ext_found = True

        if valid_ext_found == False:
            raise Exception("ERROR: Evidence file passed with unsupported file extension: {0}.  Supported extensions are {1}".format(ev_file, allowed_extensions))

    ## The input file should be defined as $path:$feattype
    if ':' not in args.reference:
        raise Exception("ERROR: input_file must be like /path/to/some.gff3:mRNA")
        
    ref_file_parts = args.reference.split(':')
    print("DEBUG: part count: {0}".format(len(ref_file_parts)))
        
    if ref_file_parts[0].endswith('.gff3'):
        (ref_assemblies, ref_features) = gff.get_gff3_features(ref_file_parts[0])
    else:
        raise Exception("ERROR: Expected input file (-i) to have a gff3 extension, got {0}".format(ref_file_parts[0]))

    ####################################################
    ## Initialize the coverage arrays

    fasta_cov = dict()
    for seq_id in fasta:
        # create a list of 0s the length of the molecule
        fasta_cov[seq_id] = [0] * len(fasta[seq_id]['s'])

    ####################################################
    ## Now parse the evidence files
        
    for ev_file in args.evidence_files:
        if ev_file.endswith('pileup'):
            parse_pileup(fasta_cov, ev_file)
        elif ev_file.endswith('sam'):
            parse_sam(fasta_cov, ev_file)
        else:
            print("INFO: ignoring evidence file {0} because code to handle its file type isn't currently implemented".format(ev_file))
        

    for id in fasta_cov:
        covered_bases = 0

        for i in fasta_cov[id]:
            if fasta_cov[id][i] > 0:
                covered_bases += 1

        fout.write("{0}\t{1}\t{2}\n".format(id, len(fasta[id]['s']), covered_bases))


def parse_pileup(cov, ev_file):
    for line in open(ev_file):
        cols = line.split("\t")
        
        if len(cols) < 5: continue
        pos = int(cols[1]) - 1

        cov[cols[0]][pos] += int(cols[3])


def parse_sam(cov, ev_file):
    line_number = 1
    
    for line in open(ev_file):
        cols = line.split("\t")
        if len(cols) < 11: continue

        cols[3] = int(cols[3])

        #if cols[3] != 0 and cols[3] != '*':
        # http://samtools.github.io/hts-specs/SAMv1.pdf
        # https://wiki.python.org/moin/BitwiseOperators
        #if int(cols[1]) & 4:
        if cols[2] != '*':
            for i in range(cols[3] - 1, cols[3] + len(cols[9]) - 1):
                try:
                    cov[cols[2]][i] += 1
                except IndexError:
                    print("WARNING: attempted to access position {0} of molecule {1}".format(i, cols[2]))
                    pass
                except KeyError:
                    print("ERROR: Molecule ID {0} wasn't found on line {1}".format(cols[2], line_number))

        line_number += 1
                    




if __name__ == '__main__':
    main()







