#!/usr/bin/env python3

"""
Prodigal's 'GFF' output is just a collection of CDS features, which doesn't follow the GFF3
spec's canonical gene model.  This script corrects that.

INPUT example:

gi|169887498|gb|CP000948.1|	Prodigal_v2.6.3	CDS	3	98	1.7	+	0	ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.427;conf=59.47;score=1.67;cscore=-1.05;sscore=2.72;rscore=0.00;uscore=0.00;tscore=3.22;
gi|169887498|gb|CP000948.1|	Prodigal_v2.6.3	CDS	337	2799	342.5	+	0	ID=1_2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.531;conf=99.99;score=342.55;cscore=326.54;sscore=16.01;rscore=10.59;uscore=2.26;tscore=3.82;

OUTPUT example:

gi|169887498|gb|CP000948.1|     Prodigal_v2.6.3 gene    3       98      1.7     +       0       ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.427;conf=59.47;score=1.67;cscore=-1.05;sscore=2.72;rscore=0.00;uscore=0.00;tscore=3.22;
gi|169887498|gb|CP000948.1|     Prodigal_v2.6.3 mRNA    3       98      1.7     +       0       ID=1_1_mRNA;Parent=1_1
gi|169887498|gb|CP000948.1|     Prodigal_v2.6.3 CDS     3       98      1.7     +       0       ID=1_1_cds;Parent=1_1_mRNA
gi|169887498|gb|CP000948.1|     Prodigal_v2.6.3 exon    3       98      1.7     +       0       ID=1_1_exon;Parent=1_1_mRNA
gi|169887498|gb|CP000948.1|     Prodigal_v2.6.3 polypeptide     3       98      1.7     +       0       ID=1_1_polypeptide;Parent=1_1_mRNA
gi|169887498|gb|CP000948.1|     Prodigal_v2.6.3 gene    337     2799    342.5   +       0       ID=1_2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.531;conf=99.99;score=342.55;cscore=326.54;sscore=16.01;rscore=10.59;uscore=2.26;tscore=3.82;
gi|169887498|gb|CP000948.1|     Prodigal_v2.6.3 mRNA    337     2799    342.5   +       0       ID=1_2_mRNA;Parent=1_2
gi|169887498|gb|CP000948.1|     Prodigal_v2.6.3 CDS     337     2799    342.5   +       0       ID=1_2_cds;Parent=1_2_mRNA
gi|169887498|gb|CP000948.1|     Prodigal_v2.6.3 exon    337     2799    342.5   +       0       ID=1_2_exon;Parent=1_2_mRNA
gi|169887498|gb|CP000948.1|     Prodigal_v2.6.3 polypeptide     337     2799    342.5   +       0       ID=1_2_polypeptide;Parent=1_2_mRNA

"""

import argparse
import os
import re


def main():
    parser = argparse.ArgumentParser( description='Converts prodigal GFF to canonical GFF3')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    ofh = open(args.output_file, 'wt')
    ofh.write("##gff-version 3\n")

    for line in open(args.input_file):
        if line.startswith('#'):
            m = re.search("gff-version", line)
            if not m:
                ofh.write(line)
            continue
        
        cols = line.split("\t")
        if len(cols) == 9:
            m = re.match("ID=(.+?);(.+)", cols[8])
            id, annotation = m.groups()

            print_row(ofh, cols, 'gene', "ID={0};{1}".format(id, annotation))
            print_row(ofh, cols, 'mRNA', "ID={0}_mRNA;Parent={0}".format(id))
            print_row(ofh, cols, 'CDS', "ID={0}_cds;Parent={0}_mRNA".format(id))
            print_row(ofh, cols, 'exon', "ID={0}_exon;Parent={0}_mRNA".format(id))
            print_row(ofh, cols, 'polypeptide', "ID={0}_polypeptide;Parent={0}_mRNA".format(id))
            
def print_row(fh, cols, ftype, col9):
    cols[2] = ftype
    fh.write( "\t".join(cols[:-1]) + "\t{0}\n".format(col9))

if __name__ == '__main__':
    main()







