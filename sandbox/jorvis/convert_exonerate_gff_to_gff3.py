#!/usr/bin/env python3

"""
The default genewise output contains GFF2 regions as long as you pass the --showquerygff or
--showtargetgff options.

The output of this script currently reduces this to cDNA_match features, but a future update
could add the ability to support gene model conventions (gene, mRNA, CDS, exon)

Note: Exonerate has a lot of different runtime modes, and I haven't tested all of them to know if
the output format changes significantly such that this script will break.  So far, it has been
tested with the following --model arguments: est2genome


Example input (only GFF regions are show, the others are ignored):

# --- START OF GFF DUMP ---
#
#
##gff-version 2
##source-version exonerate:est2genome 2.2.0
##date 2013-11-07
##type DNA
#
#
# seqname source feature start end score strand frame attributes
#
jcf7180000085759        exonerate:est2genome    gene    829164  1013996 178     -       .       gene_id 0 ; sequence comp1000_c0_seq1 ;
 gene_orientation +
jcf7180000085759        exonerate:est2genome    utr5    1013927 1013996 .       -       .       
jcf7180000085759        exonerate:est2genome    exon    1013927 1013996 .       -       .       insertions 7 ; deletions 5
jcf7180000085759        exonerate:est2genome    splice5 1013925 1013926 .       -       .       intron_id 1 ; splice_site "AT"
jcf7180000085759        exonerate:est2genome    intron  829487  1013926 .       -       .       intron_id 1
jcf7180000085759        exonerate:est2genome    splice3 829487  829488  .       -       .       intron_id 0 ; splice_site "AG"
jcf7180000085759        exonerate:est2genome    exon    829164  829486  .       -       .       insertions 17 ; deletions 33
jcf7180000085759        exonerate:est2genome    similarity      829164  1013996 178     -       .       alignment_id 0 ; Query comp1000
_c0_seq1 ; Align 1013997 114 8 ; Align 1013988 122 11 ; Align 1013977 138 16 ; Align 1013960 154 3 ; Align 1013955 157 1 ; Align 101395
3 158 15 ; Align 1013936 173 9 ; Align 829487 182 29 ; Align 829457 211 39 ; Align 829418 255 23 ; Align 829395 279 1 ; Align 829394 28
2 9 ; Align 829385 292 9 ; Align 829375 301 11 ; Align 829364 313 4 ; Align 829360 318 1 ; Align 829359 320 14 ; Align 829345 335 32 ; 
Align 829313 374 2 ; Align 829311 378 5 ; Align 829306 384 2 ; Align 829304 387 14 ; Align 829290 405 14 ; Align 829275 419 2 ; Align 8
29272 421 5 ; Align 829266 426 6 ; Align 829260 433 4 ; Align 829253 437 13 ; Align 829238 450 8 ; Align 829229 458 3 ; Align 829225 46
1 12 ; Align 829212 473 16 ; Align 829196 491 3 ; Align 829193 495 7 ; Align 829186 503 7 ; Align 829175 510 11
# --- END OF GFF DUMP ---

Example output:

jcf7180000085759        exonerate:est2genome    cDNA_match      1013927 1013996 .       -       .       ID=exonerate.align.1;Target=comp1000_c0_seq1
jcf7180000085759        exonerate:est2genome    cDNA_match      829164  829486  .       -       .       ID=exonerate.align.2;Target=comp1000_c0_seq1
jcf7180000085727        exonerate:est2genome    cDNA_match      445471  445584  .       -       .       ID=exonerate.align.3;Target=comp1000_c0_seq1
jcf7180000085727        exonerate:est2genome    cDNA_match      444121  444161  .       -       .       ID=exonerate.align.4;Target=comp1000_c0_seq1
jcf7180000085705        exonerate:est2genome    cDNA_match      6038    6142    .       +       .       ID=exonerate.align.5;Target=comp10025_c0_seq1
jcf7180000085705        exonerate:est2genome    cDNA_match      6248    6461    .       +       .       ID=exonerate.align.6;Target=comp10025_c0_seq1


Author: Joshua Orvis
Contact: jorvis AT gmail
"""

import argparse
import os
import re


def main():
    parser = argparse.ArgumentParser( description='Converts exonerate GFF output to GFF3')

    # output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to parse' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    
    args = parser.parse_args()

    fout = open(args.output_file, 'w')
    fout.write("##gff-version 3\n")

    next_id_num = 1
    in_gff3_region = False
    current_sequence = None

    for line in open(args.input_file, 'r'):
        if line.startswith('# --- START OF GFF DUMP ---'):
            in_gff3_region = True
    
        elif line.startswith('# --- END OF GFF DUMP ---'):
            in_gff3_region = False
            current_sequence = None
        
        if line.startswith('#'):
            continue
        
        cols = line.split("\t")

        if len(cols) != 9 or in_gff3_region == False:
            continue

        if cols[2] == 'gene':
            m = re.search( 'sequence (\S+)', cols[8] )

            if m:
                current_sequence = m.group(1)
            else:
                raise Exception("ERROR: expected to find a sequence identifier in the 9th column but didn't.")

        elif cols[2] == 'exon':
            id = "exonerate.align.{0}".format(next_id_num)
            next_id_num += 1

            cols[2] = 'cDNA_match'
            cols[8] = "ID={0};Target={1}".format(id, current_sequence)

            fout.write("\t".join(cols) + "\n")


if __name__ == '__main__':
    main()





















