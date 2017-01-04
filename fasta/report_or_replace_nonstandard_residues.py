#!/usr/bin/env python3

import argparse
import sys

from biocode import utils

'''
Description:

Some software is not written to handle the ambiguity residues defined by IUPAC.  The
residue "J" in a protein sequence, for example, can represent leucine or isoleucine but
can cause some software to fail.  You can use this script to either produce a report
on the non-standard residues in a multifasta file, or to optionally replace those 
residues with user-specified ones.


Input:

This script can essentially be used in two modes:

    1. Report: Define the -i (input), -t (type) and -l (list) options to read through your
       file and generate an output file showing each sequence ID which contains ambiguous
       bases and a list of those bases found (see Output section below.)

       EXAMPLE:  Create foo.report listing the non-standard residues in a nucleotide multi-fasta
       
       ./report_or_replace_nonstandard_residues.py -i file.fna -t n -l foo.report

    2. Replace: You can pass additional options to replace (-r) any given residue with (-w)
       another residue of your choosing.

       EXAMPLE: Replace all J residues with L in the protein file 'file.faa' and write results
                to a new 'file.fixed.faa'
                
       ./report_or_replace_nonstandard_residues.py -i file.faa -t p -r J -w L -o file.fixed.faa

Output:

If you define the --list option, an output file will be created where each line represents
a sequence containing ambiguous residues.  Tab separated, each additional column contains
the ambiguous residues found within along with the count of each.  Example:

    cgd6_270        B:1     J:3
    cgd3_1180       B:2
    cgd4_4160       B:1     U:1
    PVX_173270      J:1

If replacing characters using the -r and -w options, a new file defined by -o is created
with the user-specified residues replaced.  All IDs and headers are retained and sequences
residues are written 60-characters per line.

If the --print_locations flag is passed, the locations of every non-standard residue is
printed to STDERR:

EXAMPLE: $ ./report_or_replace_nonstandard_residues.py -i test.fna -t n -pl
Molecule gnl|WGS:AAGK|NW_876245.1|gb|AAGK01000005 contains residue K at position 36211
Molecule gnl|WGS:AAGK|NW_876247.1|gb|AAGK01000003 contains residue R at position 4066
Molecule gnl|WGS:AAGK|NW_876247.1|gb|AAGK01000003 contains residue W at position 4096

'''

def main():
    parser = argparse.ArgumentParser( description='Reports on non-standard characters in multifasta files and can optionally replace residues')

    parser.add_argument('-i', '--input', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-t', '--type', type=str, required=True, choices=('n', 'p'), help='Either n for nucleotide or p for protein')
    parser.add_argument('-o', '--output', type=str, required=False, help='Path to an output FASTA file to be created if doing replacement' )
    parser.add_argument('-pl', '--print_locations', dest='print_locations', action='store_true', help='If passed, will report coordinate of each non-standard residue on STDERR' )
    parser.add_argument('-r', '--replace', type=str, required=False, help='Replace this character with the one defined by --with_' )
    parser.add_argument('-w', '--with_', type=str, required=False, help='This character or set replaces all instances of the one found in --replace' )
    parser.add_argument('-l', '--list', type=str, required=False, help='Optional file of IDs where non-standard residues were detected or replaced' )
    parser.add_argument('-g', '--ignore', type=str, required=False, default='N*X', help='List of characters to not report as non-standard.  Default = the universal ambiguity bases (N, X) or the end-of-translation stop for proteins (*)' )
    parser.set_defaults(print_locations=False)
    args = parser.parse_args()

    if args.output is None:
        out_fh = sys.stdout
    else:
        out_fh = open( args.output, 'wt' )

    ## if you define --replace, you must also define --with_, and vice versa
    if args.replace is not None and args.with_ is None:
        raise Exception("ERROR: You must pass --with_ when passing --replace")
    if args.with_ is not None and args.replace is None:
        raise Exception("ERROR: You must pass --replace when passing --with_")

    seqs = utils.fasta_dict_from_file(args.input)

    ## standard characters (depends on the type of sequence)
    standard_residues = dict()
    if args.type == 'n':
        for base in list("ATGCU"):
            standard_residues[base] = 1
    else:
        for base in list("ACDEFGHIKLMNPQRSTVWY"):
            standard_residues[base] = 1

    if args.list is not None:
        list_fh = open(args.list, 'wt')

    ## build the lookup of characters to ignore
    ignore_residues = dict()
    for residue in list(args.ignore):
        ignore_residues[residue.upper()] = None

    ## process the sequences
    seqs_with_bad_chars = dict()
    
    for seq_id in seqs:
        i = 0
        seq = seqs[seq_id]
        bad_chars = dict()

        for base in list(seq['s']):
            i += 1
            ubase = base.upper()
            if ubase not in standard_residues and ubase not in ignore_residues:
                if ubase in bad_chars:
                    bad_chars[ubase] += 1
                else:
                    bad_chars[ubase] = 1

                if args.print_locations == True:
                    print("Molecule {0} contains residue {1} at position {2}".format(seq_id, ubase, i), file=sys.stderr)
        
        if args.list is not None and len(bad_chars) > 0:
            list_fh.write("{0}".format(seq_id))
            for base in bad_chars:
                list_fh.write( "\t{0}:{1}".format(base, bad_chars[base]) )

            list_fh.write("\n")
        
        if args.replace is not None:
            seq['s'] = seq['s'].replace(args.replace, args.with_)
            out_fh.write( ">{0} {1}\n".format(seq_id, seq['h']) )
            
            for i in range(0, len(seq['s']), 60):
                out_fh.write(seq['s'][i : i + 60] + "\n")

if __name__ == '__main__':
    main()







