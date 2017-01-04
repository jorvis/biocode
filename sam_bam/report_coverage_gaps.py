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

Test file:
/local/projects-t3/aplysia/read_alignment/A1_to_TOI_matches/alignment.mapped.sorted.mpileup

"""

import argparse
import sys

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Report coverage gaps/dips from a samtools mpileup file')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input mpileup file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    parser.add_argument('-f', '--fasta_file', type=str, required=True, help='Reference fasta file, against which reads were aligned.  Needed for low 3-prime end coverage' )
    parser.add_argument('-mcd', '--min_coverage_depth', type=int, required=True, help='Min coverage depth, below which is reported' )
    parser.add_argument('-mcs', '--min_coverage_span', type=int, required=False, help='Coverage window size, the avg of which is calculated for depth cutoff' )
    parser.add_argument('-eb', '--end_buffer', type=int, required=False, default=0, help='If specified, gaps this far from either end of the molecule will not be reported.' )
    args = parser.parse_args()

    # Check, this isn't ready yet:
    if args.min_coverage_span is not None:
        raise Exception("ERROR: Sorry, --min_coverage_span not yet implemented.")

    if args.output_file is None:
        out_fh = sys.stdout
    else:
        out_fh = open(args.output_file, 'wt')

    lengths = utils.fasta_sizes_from_file(args.fasta_file)

    stats = {'total_molecules': 0, 'total_bases': 0, 'depth_sum': 0}
    
    # In mpileup gaps are reported either with positions of coverage 0 OR omitted rows.
    depths = list()
    current_seq_id = None

    for line in open(args.input_file):
        contig, this_coord, base, depth = line.split("\t")[0:4]
        this_coord = int(this_coord)

        # new molecule
        if contig != current_seq_id:
            stats['total_molecules'] += 1
            
            # purge the last one
            if current_seq_id != None:
                print_spans(current_seq_id, depths, args.min_coverage_depth, out_fh, stats, lengths, args.end_buffer)

            depths = [0] * lengths[contig]
            current_seq_id = contig

        depths[this_coord - 1] = depth

    print_spans(current_seq_id, depths, args.min_coverage_depth, out_fh, stats, lengths, args.end_buffer)

    print("INFO: Total molecules: {0}".format(stats['total_molecules']), file=sys.stderr)
    print("INFO: Total bases    : {0}".format(stats['total_bases']), file=sys.stderr)
    print("INFO: Avg cov depth  : {0}x".format(int(stats['depth_sum'] / stats['total_bases'])), file=sys.stderr)
     
def print_spans(id, depths, cutoff, fh, stats, lengths, skip_end_dist):
    print("Skip end dist is: {0}".format(skip_end_dist))
    i = 0
    span_start = None
    in_span = False

    for depth in depths:
        stats['total_bases'] += 1
        stats['depth_sum'] += int(depth)

        if int(depth) < cutoff:
            if in_span == False:
                in_span = True
                span_start = i
        else:
            if in_span == True:
                if skip_end_dist is not None and skip_end_dist > 0:
                    print("DEBUG: ended first span within {0}, skipped end needs evaluating".format(id))
                    if i > skip_end_dist and i < (lengths[id] - skip_end_dist):
                        print("\tDEBUG: span with i:{0} falls within range to include".format(id))
                        fh.write("{0}\t{1}\t{2}\t{3}\n".format(id, span_start + 1, i + 1, i - span_start))
                    else:
                        print("\tDEBUG: span with i:{0} falls within range to exclude".format(id))
                else:
                    print("DEBUG: ended first span within {0}, skipped end doesn't need evaluating".format(id))
                    fh.write("{0}\t{1}\t{2}\t{3}\n".format(id, span_start + 1, i + 1, i - span_start))

            in_span = False
            span_start = None

        i += 1

if __name__ == '__main__':
    main()







