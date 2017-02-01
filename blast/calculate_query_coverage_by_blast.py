#!/usr/bin/env python3

"""

Written to see how well a set of assembled transcripts covers reference transcripts,
this script accepts a BLAST file (with -m 8 or -m 9 options) and reports the coverage
of each query.

Assumes BLAST was run like this:

  blastall -p blastn -d assembly.fasta -i query.fasta -m 9 -e 1 -o assembly.tab -a 8

Coverage of the entries within query.fasta will be reported.

# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score



"""

import argparse

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Parse BLAST output and report coverage of queries')

    ## output file to be written
    parser.add_argument('-f', '--fasta_file', type=str, required=True, help='FASTA file of query.fasta (for lengths)' )
    parser.add_argument('-b', '--blast_file', type=str, required=True, help='BLAST output using -m 8 or -m 9 options' )
    #parser.add_argument('-e', '--evalue_cutoff', type=str, required=False, help='E-value cutoff' )
    parser.add_argument('-slpc', '--subject_length_percentage_cutoff', type=int, required=False, help='Ignore hits from transcripts of length > N% relative to the query transcript' )
    parser.add_argument('-sf', '--subject_fasta', type=str, required=False, help='Only required if -slpc is passed.  FASTA file of the sub' )
    parser.add_argument('-o', '--output_base', type=str, required=True, help='base name/path for output.  Two files starting with this will be created' )
    args = parser.parse_args()

    qsizes = utils.fasta_sizes_from_file(args.fasta_file)
    
    if args.subject_fasta is not None:
        ssizes = utils.fasta_sizes_from_file(args.subject_fasta)

    all_fh = open("{0}.cov.all.perc.txt".format(args.output_base), 'wt')
    longest_fh = open("{0}.cov.longest.perc.txt".format(args.output_base), 'wt')

    current_query_id = None
    covs_all = None

    match_segments = dict()
    longest_match_length = None

    for line in open(args.blast_file):
        if line.startswith('#'): continue

        line = line.rstrip()
        cols = line.split("\t")
        query_id = cols[0]
        subj_id = cols[1]
        qstart = int(cols[6])
        qend = int(cols[7])

        if args.subject_length_percentage_cutoff is not None:
            perc_length_diff = (ssizes[subj_id] / qsizes[query_id]) * 100
            if perc_length_diff > args.subject_length_percentage_cutoff:
                continue
        
        if query_id != current_query_id:
            if current_query_id is not None:
                # report last coverage
                # first the 'all' coverage
                all_cov_perc = (qsizes[current_query_id] - covs_all.count(0)) / qsizes[current_query_id]
                all_fh.write("{0}\t{1:.1f}\n".format(current_query_id, all_cov_perc * 100))

                # now the 'longest' coverage
                longest_cov_transcript_id = None
                longest_cov_perc = None
                for tid in match_segments:
                    # calculate coverage of this tid
                    this_cov = [0] * qsizes[current_query_id]
                    for seg in match_segments[tid]:
                        for i in range(seg[0] - 1, seg[1]):
                            this_cov[i] += 1

                    this_cov_perc = (len(this_cov) - this_cov.count(0)) / len(this_cov)
                    if longest_cov_perc is None or this_cov_perc > longest_cov_perc:
                        longest_cov_perc = this_cov_perc
                        longest_cov_transcript_id = tid
                        
                print("LOG: transcript {0} covers {1} (len:{3}) best at {2:.1f}%".format(longest_cov_transcript_id, current_query_id, longest_cov_perc * 100, len(this_cov)))
                longest_fh.write("{0}\t{1:.1f}\t{2}\n".format(current_query_id, longest_cov_perc * 100, longest_cov_transcript_id))

            # now reset and init this transcript
            current_query_id = query_id
            covs_all = [0] * qsizes[query_id]
            longest_match_length = None
            match_segments = dict()

        # now handle this row
        if subj_id not in match_segments:
            match_segments[subj_id] = list()

        match_segments[subj_id].append([qstart, qend])

        for i in range(qstart - 1, qend):
            covs_all[i] += 1

    # handle last ones
    # first the 'all' coverage
    all_cov_perc = (qsizes[current_query_id] - covs_all.count(0)) / qsizes[current_query_id]
    all_fh.write("{0}\t{1:.1f}\n".format(current_query_id, all_cov_perc * 100))

    # now the 'longest' coverage
    longest_cov_transcript_id = None
    longest_cov_perc = None
    for tid in match_segments:
        # calculate coverage of this tid
        this_cov = [0] * qsizes[current_query_id]
        for seg in match_segments[tid]:
            for i in range(seg[0] - 1, seg[1]):
                this_cov[i] += 1

        this_cov_perc = (len(this_cov) - this_cov.count(0)) / len(this_cov)
        if longest_cov_perc is None or this_cov_perc > longest_cov_perc:
            longest_cov_perc = this_cov_perc
            longest_cov_transcript_id = tid

    print("LOG: transcript {0} covers {1} (len:{3}) best at {2:.1f}%".format(longest_cov_transcript_id, current_query_id, longest_cov_perc * 100, len(this_cov)))
    longest_fh.write("{0}\t{1:.1f}\n".format(current_query_id, longest_cov_perc * 100))  

    all_fh.close()
    longest_fh.close()
            
if __name__ == '__main__':
    main()







