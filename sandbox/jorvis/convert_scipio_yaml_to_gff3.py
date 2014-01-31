#!/usr/bin/env python3

"""
Yes, Scipio has an option to output GFF (but not quite GFF3).  This script can be useful if
you want GFF but forgot to use that option and have already used up 100s/1000s of hours of grid
time running the alignments already.

INPUT

The expected input is the default YAML format of Scipio, which looks like:

jgi|Rhior3|11401|RO3G_01724:
  - ID: 1557
    status: incomplete
    reason: mismatches
    prot_len: 197
    prot_start: 0
    prot_end: 155
    prot_seq: MQHRLEWAKKHQDWTVEEWRKVVFSDETKINVWGSDGCKYYWKRPNDPLQPHHLELTVKHGAGKLMMWGCITSEGPGYACQIYDGNMNANVYQHILDTTLRDTMQYYNFNWSNVYFQHDNDPKHKAKSTIDWLDNHQVRYISDWPAQSPDLNPIE
    target: jcf7180000272905
    target_len: 1303
    strand: +
    dna_start: 643
    dna_end: 1109
    matches: 148
    mismatches: 7
    undetermined: 0
    unmatched: 0
    additional: 0
    score: 0.716
    upstream: gggcatttattcgtcgattatctgatatgggcaatagttattggtggaatttgatattgttaaagggaattttataaacatccggcattttttacttgtttacatatgtttggaataacaaattttaatacgacccaaccggttgctatcataatttttcacattcttgacatgattagaatgttgccagtatataaaagcgcttacaattcctttattctaaattttaaattgtcttatttaaattatttctcttttttaccttattaccaatcgcatgcctaaacgactgccagaagaaaaagtaaacgatctcaagcaagcattaactggatccacatcaacgtacgatattgctaaggagactggagtacatgaatcaactgttagtagatactctaggcggttgttttcaaataggacaatgggctcagctggtagacctacattagtctccaaagtgacgaagcggctaatcaaaagaaaagtcattgaaggggttctaaagactgccaaagaggtacaccgggagctggttcaacttggctatgatatctcataccagtcagcaattaacgtactgaaatcaatggagttccattcttcaataaaaaaaaagaagccgtttctgacaaaagctcat
    upstream_gap: ~
    matchings:
      - type: exon
        nucl_start: 0
        nucl_end: 466
        dna_start: 643
        dna_end: 1109
        seq: atgcagcatcgattagaatgggccaagaagcaccaggattggacagtggaagaatggcgaaaagtcgtcttttctgatgaaaccaagatcaacgtatggggctctgatggttgtaagtattattggaaaagaccaagtgatcctcttcaaccgcatcatctagagtttacagtcaagcatggaggcggcaaattgatgatgtggggctgcataacatctgaaggaccaggatatgcttgtcaaatttacgatgggaacatgaacgcggatgtttatcagcatatattagacactacacttcgtgatactatggagtactacaatttcaactggtccaacatttactttcagcatgataatgatccaaagcacaaggcaaaatccacaatagcctggttagataatcatcaagtgcgctatattagcgattggccagctcagagtccggacttaaatccaatagaac
        prot_start: 0
        prot_end: 155
        seqshifts: []
        mismatchlist: [46, 56, 62, 90, 105, 114, 131]
        undeterminedlist: []
        inframe_stopcodons: []
        translation: MQHRLEWAKKHQDWTVEEWRKVVFSDETKINVWGSDGCKYYWKRPSDPLQPHHLEFTVKHGGGKLMMWGCITSEGPGYACQIYDGNMNADVYQHILDTTLRDTMEYYNFNWSNIYFQHDNDPKHKAKSTIAWLDNHQVRYISDWPAQSPDLNPIE
        overlap: ~
    stopcodon: ~
    downstream: ~
    downstream_gap: ~

OUTPUT

This was originally written to support Scipio's inclusion in evidence for EVM, which describes the expected
GFF3 format for alignments here:  http://evidencemodeler.sourceforge.net/#Preparing_inputs

Example:



"""

import argparse
import os
import sys
import yaml


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    args = parser.parse_args()

    print("WARNING: This script is incomplete and is replaced by yaml2gff.1.4.pl from the Scipio distribution")

    ## output will either be a file or STDOUT
    ofh = sys.stdout
    if args.output_file is not None:
        ofh = open(args.output_file, 'wt')

    ofh.write("##gff-version 3\n")
    
    stream = open(args.input_file, 'r')
    yamldata = yaml.load(stream)

    for prot_id in yamldata:
        matches = yamldata[prot_id]
        #print("Match: ({0})".format(prot_id))
        #print(repr(matches))

        for match in matches:
            #print(repr(match))
            for matching in match['matchings']:
                if matching['type'] == 'exon':
                    ofh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t.\t{7}\n".format( \
                            match['target'], 'Scipio', 'nucleotide_to_protein_match', \
                            matching['dna_start'], matching['dna_end'], match['score'], \
                            match['strand'], \
                            "ID={0};Target={1}".format(match['ID'], prot_id)
                            ))
        

if __name__ == '__main__':
    main()







