#!/usr/bin/env python3.2

'''
The goal here is to take any VCF file containing SNP data and generate a figure like
that in Figure 1a. of PMID: 22863733

In the paper, the X-axis is "Genotype quality" and has a range of 0-100, while the
Y-axis is "divergence from reference (SNPs per kb)"

Important link explaining expected VCF format:
http://www.broadinstitute.org/gatk/guide/article?id=1268

Assumptions
-----------
The input VCF file is filtered on the following conditions:
    - The 5th column must not be '.'
    - The 7th column must be 'PASS'
    - The 9th column are keys like: GT:AD:DP:GQ:PL with 'GQ' in 4th position
    - The sample name is value in the 10th column (of the line starting with '#CHROM')

Will NOT currently work with merged sample files.

Plot improvements:
http://matplotlib.org/users/pyplot_tutorial.html
http://matplotlib.org/api/artist_api.html#matplotlib.lines.Line2D
http://jira.igs.umaryland.edu/browse/EUK-279
'''

import argparse
import gzip
import re

import matplotlib.pyplot as plt
from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-v', '--vcf_list', type=str, required=True, help='Input list of VCF files' )
    parser.add_argument('-o', '--output', type=str, required=True, help='Base output path to be created.' )
    args = parser.parse_args()

    vcf_files = utils.read_list_file(args.vcf_list)

    for vcf_file in vcf_files:
        plot_vcf_file( vcf_file )

    plt.legend(loc='best')
    plt.xlabel("Genotype quality")
    plt.ylabel("Divergence from reference (SNPs per kb)")
    plt.savefig("{0}.png".format(args.output))


def plot_vcf_file( file ):
    fh = None
    if file.endswith(".gz"):
        fh = gzip.GzipFile(file, "r")
    else:
        fh = open(file)

    total_seq_length = 0
    snp_scores = list()
    sample_id = None
        
    for line in fh:
        line = line.decode().rstrip()
        cols = line.split("\t")

        if line.startswith("##contig="):
             m = re.match( ".+length\=(\d+)\>", line )
             if m:
                 total_seq_length += int(m.group(1))
                 continue
             else:
                 raise Exception("ERROR: found a contig line with no length defined: {0}".format(line) )

        elif line.startswith("#CHROM"):
            sample_id = cols[9]
            continue
        
        elif line.startswith("#"):
            continue
        
        cols = line.split("\t")
        
        if len(cols) < 7 or cols[4] == '.' or cols[6] != 'PASS':
            continue

        phred_score = get_scaled_phred_quality_score(cols[8], cols[9])
        snp_scores.append(phred_score)

    ## sanity check
    if total_seq_length < 1:
        raise Exception("ERROR: no sequence length lines found (like ##contig=) in file {0}".format(file))
    
    snps_per_kb = len(snp_scores) / (total_seq_length / 1000)
    print("INFO: Sample {3}, refsequence length is {0} and SNP count is {1} ({2} SNPs/kb)".format(total_seq_length, len(snp_scores), snps_per_kb, sample_id) )

    pscores = list()
    pcounts = list()
    last_score = 0

    ## counts the number of SNPs that score at least as high as each unique score
    for n,score in enumerate(sorted(snp_scores)):
        if score != last_score:
            pscores.append(score)

            snps_scoring_this_high = len(snp_scores) - n + 1
            pcounts.append( snps_scoring_this_high / (total_seq_length / 1000) )
            last_score = score

    if ( snps_per_kb < 1 ):
        this_line = plt.plot(pscores, pcounts, label=sample_id, linewidth=5.0)
    else:
        this_line = plt.plot(pscores, pcounts, label=sample_id, linewidth=2.0)
 #   plt.gca().get_frame().set_linewidth(10).draw()
    

    

def get_scaled_phred_quality_score(keys, values):
    keylist = keys.split(':')
    vallist = values.split(':')

    if keylist[3] == 'GQ':
        return float(vallist[3])
    else:
        raise Exception('ERROR: expected GQ in 4th column of this string: {0}'.format(keys))
    



if __name__ == '__main__':
    main()












    
