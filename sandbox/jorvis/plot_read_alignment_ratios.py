#!/usr/bin/env python3

"""

./plot_read_alignment_ratios.py -i /usr/local/projects/aplysia/A1_CNS_reorientation.log -o /usr/local/projects/aplysia/A1_CNS_reorientation.png
Count of ratios < 0.05: 47414
Count of ratios > 0.95: 53006

./plot_read_alignment_ratios.py -i /usr/local/projects/aplysia/A3_digestive_reorientation.log -o /usr/local/projects/aplysia/A3_digestive_reorientation.png
Count of ratios < 0.05: 44087
Count of ratios > 0.95: 49084

./plot_read_alignment_ratios.py -i /usr/local/projects/aplysia/A7_heart_reorientation.log -o /usr/local/projects/aplysia/A7_heart_reorientation.png
Count of ratios < 0.05: 45995
Count of ratios > 0.95: 52188

./plot_read_alignment_ratios.py -i /usr/local/projects/aplysia/A10_gills_reorientation.log -o /usr/local/projects/aplysia/A10_gills_reorientation.png
Count of ratios < 0.05: 49941
Count of ratios > 0.95: 55683

Ranges wanted:
0.125 - 0.875
0.25 - 0.75
0.475 - 0.525

"""

import argparse
import re

import matplotlib
import matplotlib.pyplot as plt
from biocode import utils

matplotlib.use('Agg')
import sys


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here' )

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input reorientation log file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output image file to be created' )
    parser.add_argument('-f', '--fasta_file', type=str, required=True, help='Path to an input FASTA file' )
    parser.add_argument('-rmin', '--ratio_min', type=float, required=False, default=0.05, help='Left bounds of ratio with correctly mapped reads to export' )
    parser.add_argument('-rmax', '--ratio_max', type=float, required=False, default=0.95, help='Right bounds of ratio with correctly mapped reads to export' )
    args = parser.parse_args()

    ratios = list()

    # if set to true, the IDs in the mid-range will be printed to STDOUT
    print_ids = True
    
    LENGTH_CUTOFF = 350
    ratio_min_count = 0
    ratio_bet_count = 0
    ratio_max_count = 0

    fasta = utils.fasta_dict_from_file(args.fasta_file)

    for line in open(args.input_file):
        # lines are like: comp0_c0_seq1   1-T:6   1-F:0   2-T:0   2-F:5
        m = re.search('(.+)\t1-T:(\d+)\t1-F:(\d+)\t2-T:(\d+)\t2-F:(\d+)', line)
        if m:
            seq_id = m.group(1)

            if seq_id in fasta:
                if len(fasta[seq_id]['s']) < LENGTH_CUTOFF:
                    continue
            else:
                raise Exception("Expected but filed to find seq ID {0} in FASTA file".format(seq_id))
            
            f_reads_correctly_mapped = int(m.group(2))
            f_reads_incorrectly_mapped = int(m.group(3))
            r_reads_correctly_mapped = int(m.group(5))
            r_reads_incorrectly_mapped = int(m.group(4))
            f_read_count = f_reads_correctly_mapped + f_reads_incorrectly_mapped

            if f_read_count > 0:
                correct_ratio = f_reads_correctly_mapped / f_read_count
                ratios.append(correct_ratio)

                if correct_ratio < args.ratio_min:
                    ratio_min_count += 1
                elif correct_ratio > args.ratio_max:
                    ratio_max_count += 1
                else:
                    ratio_bet_count += 1
                    if print_ids == True:
                        print(seq_id)

                #print("LOG: Fcorrect:{0} Fwrong:{1} Ftotal:{2} ratio:{3}".format(f_reads_correctly_mapped, f_reads_incorrectly_mapped, f_read_count, correct_ratio))

    plt.hist(ratios, bins=100)
    plt.xlabel("Correct read orientation alignment ratio")
    plt.ylabel("Log of transcript count")
    plt.grid(True)
    #plt.ylim(0,5000)
    plt.gca().set_yscale("log")
    plt.savefig(args.output_file)

    sys.stderr.write("Length cutoff: {0} bp\n".format(LENGTH_CUTOFF))
    sys.stderr.write("Count of ratios < {0}: {1}\n".format(args.ratio_min, ratio_min_count))
    sys.stderr.write("Count where {0} > ratio < {1}: {2}\n".format(args.ratio_min, args.ratio_max, ratio_bet_count))
    sys.stderr.write("Count of ratios > {0}: {1}\n".format(args.ratio_max, ratio_max_count))

    

if __name__ == '__main__':
    main()







