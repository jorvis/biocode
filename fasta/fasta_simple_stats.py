#!/usr/bin/env python3

"""
Reads one or more (optionally gz compressed) FASTA files and reports some basic statistics.  Example:

$ ./fasta_simple_stats.py Trinity.fasta

File: Trinity.fasta
	Total sequence entries: 583596
	Total bases: 425253191
	Shortest sequence length: 177.0 (TRINITY_DN32084_c0_g1_i3)
	Average sequence length : 728.7
	Longest sequence length : 42580.0 (TRINITY_DN1050_c0_g1_i9)
	GC percentage: 43.2

If you pass multiple files via *.fasta or something similar this will report stats
on these as an aggregate unless you pass the --individual flag, which will cause
reporting on each file individually.

Contact: jorvis@gmail.com
"""

import argparse
import codecs
import gzip
import numpy
import os
import sys
import re

def main():
    parser = argparse.ArgumentParser( description='Provides simple quantitative statistics for one or more FASTA files')

    ## output file to be written
    parser.add_argument('input_files', metavar='N', type=str, nargs='+', help='Path to one or more input files, separated by spaces' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Optional path to an output file to be created, else prints on STDOUT' )
    parser.add_argument('--individual', action='store_true', help='Report stats on each file individually' )
    args = parser.parse_args()

    ## open the output file
    fout = None
    if args.output_file is None:
        fout = codecs.getwriter('utf8')(sys.stdout.buffer)
    else:
        fout = open(args.output_file, "w")

    ## values that will be reported
    data = {
        'entry_count': 0,
        'total_bases': 0,
        'gc_count': 0,
        'smallest_length': None,
        'smallest_id': None,
        'largest_length': 0,
        'largest_id': None,
        'lengths': list()
    }

    for input_file in args.input_files:
        if input_file.endswith('.gz'):
            fh = gzip.open(input_file, 'rb')
            is_compressed = True
        else:
            fh = open(input_file, 'r')
            is_compressed = False

        if args.individual:
           data = {
               'entry_count': 0,
               'total_bases': 0,
               'gc_count': 0,
               'smallest_length': None,
               'smallest_id': None,
               'largest_length': 0,
               'largest_id': None,
               'lengths': list()
           }

        line_number = 0
        this_seq_length = 0
    
        for line in fh:
            line_number += 1
            if is_compressed:
                line = line.decode()
            
            line = line.rstrip()

            ## every line starting with > is the start of a sequence entry
            if line and line[0] == '>':
                data['entry_count'] += 1

                if line_number != 1:
                    if data['smallest_id'] is None:
                        # this is the first record, so it sets all records.
                        data['smallest_id'] = seq_id
                        data['smallest_length'] = this_seq_length
                    elif this_seq_length < data['smallest_length']:
                        data['smallest_id'] = seq_id
                        data['smallest_length'] = this_seq_length

                    if this_seq_length > data['largest_length']:
                        data['largest_length'] = this_seq_length
                        data['largest_id'] = seq_id

                    data['lengths'].append(this_seq_length)

                m = re.match(">(\S+)", line)
                seq_id = m.group(1)
                
                this_seq_length = 0
            else:
                data['total_bases'] += len(line)
                this_seq_length += len(line)
                data['gc_count'] += line.count('G')
                data['gc_count'] += line.count('C')

        if args.individual:
            report_stats(fout, input_file, data)

    if not args.individual:
        if len(args.input_files) == 1:
            report_stats(fout, args.input_files[0], data)
        else:
            report_stats(fout, 'multiple', data)

def report_stats(fout, input_file, data):
    avg_entry_length = data['total_bases'] / data['entry_count']

    # calculate N50
    data['lengths'] = sorted(data['lengths'], reverse=True)
    length_sum = numpy.cumsum(data['lengths'])
    n_2 = int(sum(data['lengths'])/2)
    csum_n2 = min(length_sum[length_sum >= n_2])
    n50_idx = numpy.where(length_sum == csum_n2)
    n_50 = data['lengths'][int(n50_idx[0])]

    val90 = int(sum(data['lengths']) * 0.90)
    csum_n90 = min(length_sum[length_sum >= val90])
    n90_idx = numpy.where(length_sum == csum_n90)
    n_90 = data['lengths'][int(n90_idx[0])]

    fout.write("File: {0}\n".format(input_file))
    fout.write("\tTotal sequence entries: {0}\n".format(data['entry_count']))
    fout.write("\tTotal bases: {0}\n".format(data['total_bases']))
    fout.write("\tShortest sequence length: {0:.1f} ({1})\n".format(data['smallest_length'], data['smallest_id']))
    fout.write("\tAverage sequence length : {0:.1f}\n".format(avg_entry_length))
    fout.write("\tN50 : {0}\n".format(n_50))
    fout.write("\tN90 : {0}\n".format(n_90))
    fout.write("\tLongest sequence length : {0:.1f} ({1})\n".format(data['largest_length'], data['largest_id']))
    fout.write("\tGC percentage: {0:.1f}\n".format( (data['gc_count'] / data['total_bases']) * 100))

if __name__ == '__main__':
    main()







