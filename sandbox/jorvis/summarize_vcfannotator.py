#!/usr/bin/env python3.2

'''
The goal here is to take any VCF file output by VCFannotator and summarize the SNPs
contained there by type.  

INPUT

The input expected is from the VCFannotator documentation:

An example output line transposed to a column format would look like so (taken from 
the sample data):

0   TY-2482_chromosome
1   5080
2   .
3   T
4   C
5   8101.55
6   PASS
7   AC=2;AF=1.00;AN=2;DP=212;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=0.0000;MQ=59.52;MQ0=0;QD=38.21;SB=-4130.96
8   GT:AD:DP:GQ:PL
9   1/1:0,212:212:99:8101,605,0
10  CDS,tmpTY2482_00008,tmpTY2482_00008T0,microcin H47 immunity protein mchI,trans_orient:+,loc_in_cds:46,codon_pos:1,codon:Tct-Cct,pep:S->P,Ser-16-Pro,(NSY)

The final field is bundled with the individual feature annotation data, comma-delimited. 
This includes the feature type (eg. CDS, intron, UTR), gene and transcript identifiers, 
name of the gene, and the transcribed orientation of the gene. If the SNP is localized 
to a coding region, then the relative position within that CDS sequence is provided in 
addition to the codon change. The types of coding mutations are provided as synonomous 
(SYN), non-synonomous (NSY), read-thru (RTH), and nonsence (STP). SNPs that are localized 
to intergenic regions are reported as such along with the identifiers of the neighboring 
genes and distance to each.

This isn't directly followed by Brian's VCFannotator, but for my own notes these are the SNP types:

- Non-coding region
    - Intergenic
    - 5' UTR
    - 3' UTR
- Coding region
    - Synonymous
    - Nonsynonymous
        - Missense (results in a different amino acid, but not one of the three below)
        - Nonsense (results in a premature stop codon)
        - Read-through (releases a stop codon)
        - Initiating (translation-initiating fMet is changed)


'''


import argparse
import gzip


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-v', '--vcf_file', type=str, required=True, help='Input VCFannotated file' )
    args = parser.parse_args()

    file = args.vcf_file
    file_is_encoded = False

    fh = None
    if file.endswith(".gz"):
        fh = gzip.GzipFile(file, "r")
        file_is_encoded = True
    else:
        fh = open(file)

    total_snp_c = 0
    cds_snp_c = 0
    intron_snp_c = 0
    intergenic_snp_c = 0
    p3UTR_snp_c = 0
    p5UTR_snp_c = 0

    ## these are each subcategories of SNPs within a CDS
    syn_c = 0
    nonsyn_c = 0
    readthru_c = 0
    nonsense_c = 0
    other_c = 0

    for line in fh:
        if file_is_encoded:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
            
        cols = line.split("\t")

        if line.startswith("#"):
            continue
        
        cols = line.split("\t")
        
        ## VCF annotated files should all have 10 columns.
        if len(cols) < 10:
            raise Exception("ERROR: expected all non-comment lines to have 10 columns")

        annot_col = cols[10]
        snp_loc = None
        snp_type = None
        alt = cols[4]

        if alt != ".":
            total_snp_c += 1
            if annot_col.startswith("intergenic"):
                intergenic_snp_c += 1
            elif annot_col.startswith("intron"):
                intron_snp_c += 1
            elif annot_col.startswith("p5UTR"):
                p5UTR_snp_c += 1
            elif annot_col.startswith("p3UTR"):
                p3UTR_snp_c += 1
            elif annot_col.startswith("CDS"):
                cds_snp_c += 1
                cds_type = annot_col.split(",")[-1]
                if cds_type == "(SYN)":
                    syn_c += 1
                elif cds_type == "(NSY)":
                    nonsyn_c += 1
                elif cds_type == "(RTH)":
                    readthru_c += 1
                elif cds_type == "(STP)":
                    nonsense_c += 1
                else:
                    other_c += 1
            
            
            else:
                raise Exception("ERROR: Unexpected SNP type at beginning of column: {0}".format(annot_col) )
        
            
    print("Total SNPs: {}".format(total_snp_c) )
    print("Intergenic: {}".format(intergenic_snp_c) )
    print("Intronic  : {}".format(intron_snp_c) )
    print("p5UTR  : {}".format(p5UTR_snp_c) )
    print("p3UTR  : {}".format(p3UTR_snp_c) )
    print("Within CDS: {}".format(cds_snp_c) )
    print("\tSynonymous    : {}".format(syn_c) )
    print("\tNon-synonymous: {}".format(nonsyn_c) )
    print("\tRead-through  : {}".format(readthru_c) )
    print("\tNonsense      : {}".format(nonsense_c) )



if __name__ == '__main__':
    main()












    
