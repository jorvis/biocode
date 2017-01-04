#!/usr/bin/env python3

"""

INPUT

Expected input file format (pileup):

Each line consists of 5 (or optionally 6) tab-separated columns:

  1.  Sequence identifier
  2.  Position in sequence (starting from 1)
  3.  Reference nucleotide at that position
  4.  Number of aligned reads covering that position (depth of coverage)
  5.  Bases at that position from aligned reads
  6.  Quality of those bases (OPTIONAL)

In my testing, this was generated like this:

  samtools mpileup -f Trinity.fasta XUMTA_20131112.bowtie.sorted.mappedonly.bam > XUMTA_20131112.bowtie.sorted.mappedonly.mpileup

  
OUTPUT

 The X and Y size of the resulting image is going to be a product of the --mol_size_limit and --mol_bin_size
 parameters.

 If you pass a value of 'plot' to the -o parameter it will invoke the interactive plot viewer rather
 than writing an output file.  (You can still save a file from within the viewer)
  
"""

import argparse
import numpy as np
from collections import defaultdict

import matplotlib.pyplot as plt
from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Generates a figure showing coverage/abundance vs. molecule size.')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input pileup file' )
    parser.add_argument('-f', '--fasta_file', type=str, required=True, help='Path to the FASTA file of reference molecules' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-s', '--mol_size_limit', type=int, required=False, default=5000, help='Results for molecules over this size will be grouped together' )
    parser.add_argument('-b', '--mol_bin_size', type=int, required=False, default=10, help='Set the binning resolution of the transcript size axis')
    args = parser.parse_args()

    ## first, we need a collection of the FASTA data and the molecule lengths
    molecules = utils.fasta_dict_from_file(args.fasta_file)

    ## data points for plotting
    #  structure like this:
    #    500 = { 30 => 2 }
    #  which means: There were 2 transcripts with median coverage of 30 and length between 500 and 500+mol_bin_size
    data_bins = defaultdict(lambda: defaultdict(int))

    current_molecule_id = None
    current_molecule_coverages = list()

    ## These files are usually huge.  For scalability, operations performed within this
    #  loop should be limited.
    for line in open(args.input_file):
        cols = line.split("\t")

        if current_molecule_id is None:
            current_molecule_id = cols[0]
            current_molecule_coverages = [0] * len(molecules[cols[0]]['s'])

        if cols[0] != current_molecule_id:
            mol_length_bin = int(len(molecules[current_molecule_id]['s']) / args.mol_bin_size)
            median_size = np.median(current_molecule_coverages)
            data_bins[mol_length_bin][median_size] += 1
            print("DEBUG: molecule {0} appeared to be {1} bp in length with median coverage of {2}".format(current_molecule_id, len(molecules[current_molecule_id]['s']), median_size))

            # reset
            current_molecule_id = cols[0]
            current_molecule_coverages = [0] * len(molecules[cols[0]]['s'])

        try:
            current_molecule_coverages[int(cols[1]) - 1] = int(cols[3]) 
        except IndexError:
            print("ERROR: pileup file reports position {0} coverage but transcript {1} is only {2} bp in length".format(cols[1], current_molecule_id, len(molecules[cols[0]]['s'])) )

    # don't forget the last one
    mol_length_bin = int(len(molecules[cols[0]]['s']) / args.mol_bin_size)
    median_size = np.median(current_molecule_coverages)
    data_bins[mol_length_bin][median_size] += 1

    ## now generate the plot data - x,y positions and radii
    x = list()
    y = list()
    r = list()

    for bin_size in data_bins:
        for cov in data_bins[bin_size]:
            x.append(bin_size)
            y.append(cov)
            r.append(data_bins[bin_size][cov])

    plt.xlabel('Molecule length')
    plt.ylabel('Median depth of coverage')
    #plt.xlim(0,2000)
    #plt.ylim(0,500)
    plt.scatter(x, y, s=r, alpha=0.5)

    if args.output_file == 'plot':
        plt.show()
    else:
        plt.savefig(args.output_file)
    


if __name__ == '__main__':
    main()







