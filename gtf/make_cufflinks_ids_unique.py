#!/usr/bin/env python3

"""
The default GTF format from cufflinks numbers IDs sequentially on each molecule
so there will be duplications.  For example:

jcf7180000645906	Cufflinks	transcript	5663	7071	308	-	.	gene_id "CUFF.2"; transcript_id "CUFF.2.1"; FPKM "3.3112927014"; frac "0.226898"; conf_lo "1.089987"; conf_hi "5.532598"; cov "3.114901";
jcf7180000645906	Cufflinks	exon	5663	5812	308	-	.	gene_id "CUFF.2"; transcript_id "CUFF.2.1"; exon_number "1"; FPKM "3.3112927014"; frac "0.226898"; conf_lo "1.089987"; conf_hi "5.532598"; cov "3.114901";
jcf7180000645906	Cufflinks	exon	5874	6162	308	-	.	gene_id "CUFF.2"; transcript_id "CUFF.2.1"; exon_number "2"; FPKM "3.3112927014"; frac "0.226898"; conf_lo "1.089987"; conf_hi "5.532598"; cov "3.114901";
jcf7180000645906	Cufflinks	exon	6219	6353	308	-	.	gene_id "CUFF.2"; transcript_id "CUFF.2.1"; exon_number "3"; FPKM "3.3112927014"; frac "0.226898"; conf_lo "1.089987"; conf_hi "5.532598"; cov "3.114901";
jcf7180000645906	Cufflinks	exon	6418	6858	308	-	.	gene_id "CUFF.2"; transcript_id "CUFF.2.1"; exon_number "4"; FPKM "3.3112927014"; frac "0.226898"; conf_lo "1.089987"; conf_hi "5.532598"; cov "3.114901";
jcf7180000645906	Cufflinks	exon	6909	7071	308	-	.	gene_id "CUFF.2"; transcript_id "CUFF.2.1"; exon_number "5"; FPKM "3.3112927014"; frac "0.226898"; conf_lo "1.089987"; conf_hi "5.532598"; cov "3.114901";
jcf7180000645901	Cufflinks	transcript	1303	1815	1000	+	.	gene_id "CUFF.2"; transcript_id "CUFF.2.1"; FPKM "20.6408485852"; frac "1.000000"; conf_lo "18.345059"; conf_hi "22.936638"; cov "304.191554";
jcf7180000645901	Cufflinks	exon	1303	1419	1000	+	.	gene_id "CUFF.2"; transcript_id "CUFF.2.1"; exon_number "1"; FPKM "20.6408485852"; frac "1.000000"; conf_lo "18.345059"; conf_hi "22.936638"; cov "304.191554";
jcf7180000645901	Cufflinks	exon	1476	1573	1000	+	.	gene_id "CUFF.2"; transcript_id "CUFF.2.1"; exon_number "2"; FPKM "20.6408485852"; frac "1.000000"; conf_lo "18.345059"; conf_hi "22.936638"; cov "304.191554";
jcf7180000645901	Cufflinks	exon	1633	1815	1000	+	.	gene_id "CUFF.2"; transcript_id "CUFF.2.1"; exon_number "3"; FPKM "20.6408485852"; frac "1.000000"; conf_lo "18.345059"; conf_hi "22.936638"; cov "304.191554";

These are made unique by preprending the molecule identifier onto the IDs.  So:

    gene_id "CUFF.2"; transcript_id "CUFF.2.1

becomes:

    gene_id "jcf7180000645901.CUFF.2"; transcript_id "jcf7180000645901.CUFF.2.1"

Author: Joshua Orvis (jorvis AT gmail)
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser( description='Uniquifies IDs within a GTF file by prepending molecule identifiers')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input GTF file' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output GTF file to be created' )
    args = parser.parse_args()

    ## output will either be a file or STDOUT
    ofh = sys.stdout
    if args.output_file is not None:
        ofh = open(args.output_file, 'wt')

    assemblies = dict()

    for line in open(args.input_file, "r"):
        cols = line.split("\t")

        if len(cols) != 9:
            fout.write(line)
            continue
        
        mol_id = cols[0]

        if mol_id not in assemblies:
            assemblies[mol_id] = dict()

        cols[8] = cols[8].replace(' "CUFF.', " \"{0}.CUFF.".format(mol_id) )
        ofh.write( "{0}\n".format("\t".join(cols)) )
        
        

if __name__ == '__main__':
    main()







