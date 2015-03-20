#!/usr/bin/env python3

import argparse
import matplotlib

# back-end options are here: http://matplotlib.sourceforge.net/faq/usage_faq.html#what-is-a-backend
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import os
import re


def fasta_entry_sizes(file):
    seq_lengths = []
    
    ## these are reset as each seq entry is encountered
    seq_lines = []
    seq_count = 0

    for line in open(file, 'r'):
        if re.match('^\>', line):
            seq_count += 1
            seq = ''.join(seq_lines)
            seq = re.sub(r'\s', '', seq)
            seq_lengths.append( len(seq) )
            seq_lines = []
        else:
            seq_lines.append(line)

    return seq_lengths
        

def get_legend_labels(label_arg, file_count):
    labels = []

    if label_arg is not None:
        labels = label_arg.split(',')

        if len(labels) != file_count:
            raise Exception("Error: number of input files doesn't match number of labels specified in --legend_names")

    return labels


def main():
    parser = argparse.ArgumentParser( description='Generate FASTA file(s) size distribution plot')

    parser.add_argument('fasta_files', metavar='N', type=str, nargs='+', help='Pass one or more FASTA files')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-t', '--title', type=str, required=False, default='FASTA size distribution', \
                        help='Pass a title for the graph')
    parser.add_argument('-b', '--bin_count', type=int, required=False, default=30, \
                        help='Data will be placed into this many bins.  This is the default behavior. ' + \
                        'Alternatively, use --bin_size and --bin_max')
    parser.add_argument('-s', '--bin_size', type=int, required=False, default=0, \
                        help='Instead of a --bin_count, use this to specify the size of your bins.')
    parser.add_argument('-m', '--bin_max', type=int, required=False, default=0, \
                        help='If specifying --bin_size, you can optionally use this to limit the ' + \
                        'maximum bound of your bins (prevents long tails in plots)')
    parser.add_argument('-l', '--legend_names', type=str, required=False, help='For a legend with labels ' + \
                        'of each of your datasets, pass a comma-separated list with no spaces.')
    parser.add_argument('-g', '--log_scale', type=bool, required=False, default=False, help='Set to true for log10 Y scale')
    args = parser.parse_args()

    data_ranges = []
    fasta_files = args.fasta_files

    size_max = 0
    seqs_above_size_max = 0

    for fasta_file in fasta_files:
        print("INFO: parsing seq lengths in file: {0}".format(fasta_file))
        sizes = fasta_entry_sizes(fasta_file)
        print("INFO: {0} sequences found in {1}".format(len(sizes), fasta_file))
        data_ranges.append(sizes)

        this_max_size = max(sizes)

        if this_max_size > size_max:
            size_max = this_max_size

    
    ## calculate the bins.  default is to use the bin_count
    bin_opt = args.bin_count

    if args.bin_size != 0:
        bin_opt = []
        for bin_min in range(args.bin_size, size_max, args.bin_size):
            if args.bin_max == 0 or args.bin_max > bin_min:
                bin_opt.append(bin_min)

                if args.bin_max > bin_min:
                    seqs_above_size_max += 1
    
    plot.xlabel('Sequence size (bins)')
    plot.ylabel('Sequence counts')
    plot.title(args.title)

    if args.log_scale == True:
        n, bins, patches = plot.hist(data_ranges, bin_opt, normed=0, histtype='bar', log=True)
    else:
        n, bins, patches = plot.hist(data_ranges, bin_opt, normed=0, histtype='bar')

    plot.grid(True)

    legend_labels = get_legend_labels( args.legend_names, len(fasta_files) )
    if len(legend_labels) > 0:
        plot.legend(legend_labels)
    
    plot.savefig(args.output_file)

    print("INFO: there were {0} data points above the range defined in the histogram".format(seqs_above_size_max))


if __name__ == '__main__':
    main()







