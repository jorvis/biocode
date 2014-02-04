#!/usr/bin/env python3

"""
I was given data generated using this protocol:

    http://trinityrnaseq.sourceforge.net/analysis/abundance_estimation.html

Specifically, this file:

    RSEM.isoforms.results

OUTPUT

 If you pass a value of 'plot' to the -o parameter it will invoke the interactive plot viewer rather
 than writing an output file.  (You can still save a file from within the viewer)
  
"""

import argparse
import biocodeutils
import os
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser( description='Generates a figure showing coverage/abundance vs. molecule size.')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input pileup file' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-t', '--title', type=str, required=False, help='Plot title' )
    parser.add_argument('-g', '--graph_using', type=str, required=False, default='TPM', help='TPM or FPKM' )
    args = parser.parse_args()

    x = list()
    y = list()
    skipped_datapoints = 0
    total_datapoints = 0

    if args.graph_using not in ['TPM', 'FPKM']:
        raise Exception("ERROR: --graph_using value must be either TPM or FPKM")

    for line in open(args.input_file):
        cols = line.split("\t")

        if cols[0] == 'transcript_id' and cols[1] == 'gene_id':
            continue

        tpm = float(cols[5])
        fpkm = float(cols[6])
        total_datapoints += 1
        
        if args.graph_using == 'TPM':
            if tpm > 0:
                x.append(int(cols[2]))
                y.append(tpm)
            else:
                skipped_datapoints += 1
        else:
            if fpkm > 0:
                x.append(int(cols[2]))
                y.append(fpkm)
            else:
                skipped_datapoints += 1

    if args.graph_using == 'TPM':
        print("LOG: {0}/{1} data points were skipped because the TPM value was 0.0".format(skipped_datapoints, total_datapoints))
    else:
        print("LOG: {0}/{1} data points were skipped because the FPKM value was 0.0".format(skipped_datapoints, total_datapoints))
        
    fig = plt.figure()

    if args.title is not None:
        fig.suptitle(args.title)
                        
    plt.xlabel('Molecule length')

    if args.graph_using == 'TPM':
        plt.ylabel('log TPM')
    else:
        plt.ylabel('log FPKM')
        
    ax = plt.gca()
    ax.plot(x, y, 'o', c='blue', alpha=0.05, markeredgecolor='none')
    ax.set_yscale('log')

    if args.output_file == 'plot':
        plt.show()
    else:
        plt.savefig(args.output_file)
    


if __name__ == '__main__':
    main()







