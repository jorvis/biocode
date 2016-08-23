#!/usr/bin/env python3

"""
If you have an mpileup file this script scans it to report any gaps (or dips)
in coverage levels.

It can report any positions where coverage drops below a certain level, for example, using
the --min_coverage_depth (-mcd) option to report any spans where the coverage drops below 5x:

  $ ./report_coverage_gaps.py -i foo.pileup -mcd 5

If you want to use a sliding window for coverage calculation, using the --min_coverage_span (-mcs)
option.  Here, it will report any windows of at least 10 bp where the average coverage is below
6 bp.

  $ ./report_coverage_gaps.py -i foo.pileup -mcd 6 -mcs 10

Example command before using this script:

  $ samtools mpileup -f toi_20150925.fna alignment.mapped.sorted.bam > alignment.mapped.sorted.mpileup

NOTE: This script accommodates for the fact that the samtools mpileup output does not contain
coverage information for every base position.  If a position is omitted, it is assumed to
be a zero-depth locus.  This is also why the reference fasta file (-f) is required, in order to
be able to report 3' end sections.
"""

import argparse
import biocodeutils
import os
import sys

def main():
    parser = argparse.ArgumentParser( description='Report coverage gaps/dips from a samtools mpileup file')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    parser.add_argument('-f', '--fasta_file', type=str, required=True, help='Reference fasta file, against which reads were aligned.  Needed for low 3-prime end coverage' )
    parser.add_argument('-mcd', '--min_coverage_depth', type=int, required=True, help='Min coverage depth, below which is reported' )
    parser.add_argument('-mcs', '--min_coverage_span', type=int, required=False, default=0, help='Coverage window size, the avg of which is calculated for depth cutoff' )
    args = parser.parse_args()

    # Check, this isn't ready yet:
    if args.min_coverage_span is not None:
        raise Exception("ERROR: Sorry, --min_coverage_span not yet implemented.")

    if args.output_file is None:
        out_fh = sys.stdout
    else:
        out_fh = open(args.output_file, 'wt')

    lengths = biocodeutils.fasta_sizes_from_file(args.fasta_file)
        
    current_seq_id = None
    next_seq_pos = 1
    low_cov_start = None
    low_cov_end = None

    stats_total_molecules = 0
    stats_total_bases = 0
    stats_depth_sum = 0
    
    for line in open(args.input_file):
        stats_total_bases += 1
        contig, coord, base, depth = line.split("\t")[0:4]
        stats_total_bases += 1
        stats_depth_sum += int(depth)
        coord = int(coord)
        print("DEBUG: contig:{0} coord:{1} depth:{2} == current_seq_id:{3} low_cov_start:{5} low_cov_end:{4} next_seq_pos:{6}".format(
            contig, coord, depth, current_seq_id, low_cov_end, low_cov_start, next_seq_pos))

        if contig != current_seq_id:
            # if this isn't the end of the contig, report until the end
            if current_seq_id is not None and next_seq_pos - 1 != lengths[current_seq_id]:
                print("DEBUG: Found an extended end (last_coord:{2}) at the end of {0} which is actually {1} bp".format(current_seq_id, lengths[current_seq_id], next_seq_pos - 1))
            
            # purge and reset
            if low_cov_start is not None:
                out_fh.write("{2}\t{0}\t{1}\n".format(low_cov_start, low_cov_end, current_seq_id))
                low_cov_start = None
                
            current_seq_id = contig
            next_seq_pos = 1
            stats_total_molecules += 1

        # The mpileup output has gaps, so we need to check each coordinate and compare with the current
        #  line to see if we've found one.
        if coord < next_seq_pos:
            if low_cov_start is None:
                low_cov_start = next_seq_pos

            next_seq_pos = coord + 1
        else:
            next_seq_pos += 1

        if int(depth) < args.min_coverage_depth:
            # We're starting a new low-cov region
            if low_cov_start is None:
                low_cov_start = coord
            else:
                # We're just extending upon an open low-cov region.  Nothing to see here.
                pass

            low_cov_end = int(coord)
        else:
            # did we finish a low_cov_region?:
            if low_cov_start is not None:
                out_fh.write("Low coverage region: {2}: {0} - {1}\n".format(low_cov_start, low_cov_end, contig))
                low_cov_start = None

    # Handle possible spans at the end of the last contig
    if low_cov_start is not None:
        print("{2}\t{0}\t{1}".format(low_cov_start, low_cov_end, current_seq_id))

    print("INFO: Total molecules: {0}".format(stats_total_molecules), file=sys.stderr)
    print("INFO: Total bases    : {0}".format(stats_total_bases), file=sys.stderr)
    print("INFO: Avg cov depth  : {0}x".format(int(stats_depth_sum / stats_total_bases)), file=sys.stderr)
        

if __name__ == '__main__':
    main()







