#!/usr/bin/env python3

"""
Converts the btab output of Ergatis ncbi-blast to GFF3

Example input:
GluR_4          3300    BLASTN  Trinity.fasta   comp93586_c0_seq1       26      2208    2180    1       99.6    0.0     2139    4240            len=2337 path=[1:0-1874 1876:1875-1882 1884:1883-2336]  1       Plus    2337    0.0     0.0
GluR_4          3300    BLASTN  Trinity.fasta   comp103174_c0_seq1      2415    3277    2105    1245    98.7    0.0     819     1624            len=2105 path=[2313:0-204 2518:205-302 2616:303-366 2680:367-367 2681:368-390 2704:391-1395 3710:1396-1752 4067:1753-1763 4078:1764-2001 4316:2002-2081 5675:2082-2104] 1       Plus    2105    0.0     0.0
GluR_4          3300    BLASTN  Trinity.fasta   comp103174_c0_seq2      2186    2439    1       254     100.0   100.0   254     504             len=254 path=[1:0-168 170:169-200 202:201-253]  1       Plus    254     1e-140  1e-140
PICK1           1242    BLASTN  Trinity.fasta   comp3011803_c0_seq1     45      170     159     34      90.5    0.0     78      155             len=517 path=[495:0-516]        1       Plus    517     6e-36   6e-36
RICTOR          5676    BLASTN  Trinity.fasta   comp102759_c0_seq1      182     4349    6303    2136    99.1    0.0     4048    8025            len=6303 path=[1:0-1699 1701:1700-1706 1708:1707-4017 4019:4018-4024 4026:4025-4070 4072:4071-4075 4077:4076-6302]      1       Plus    6303    0.0     0.0
RICTOR          5676    BLASTN  Trinity.fasta   comp102759_c0_seq1      4428    5676    2057    809     99.9    0.0     1245    2468            len=6303 path=[1:0-1699 1701:1700-1706 1708:1707-4017 4019:4018-4024 4026:4025-4070 4072:4071-4075 4077:4076-6302]      1

The columns are:

    1   query_name
    2   date
    3   query_length
    4   algorithm
    5   database_name
    6   hit_name
    7   qry_start
    8   qry_end
    9   hit_start
    10  hit_end
    11  percent_identity
    12  percent_similarity
    13  raw_score
    14  bit_score
    15  NULL
    16  hit_description
    17  blast_frame
    18  qry_strand (Plus | Minus)
    19  hit_length
    20  e_value
    21  p_value


Example output:



WARNING:
Currently only tested with blastn

"""

import argparse
import os
from operator import itemgetter

next_ids = {'match':1, 'match_part':1 }

def main():
    parser = argparse.ArgumentParser( description='NCBI-BLAST (btab) converter to GFF3 format')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to parse' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-p', '--perc_identity_cutoff', type=float, required=False, help='Filters on the perc identity of each HSP' )
    
    args = parser.parse_args()
    algorithm = None
    current_qry_id = None
    current_hit_id = None
    current_match_parts = list()

    ofh = open(args.output_file, 'w')
    ofh.write("##gff-version 3\n")
    
    for line in open(args.input_file, 'r'):
        cols = line.split("\t")
        qry_id = cols[0]
        hit_id = cols[5]
        
        segment = { 'qry_id':qry_id,  'qry_start':int(cols[6]), 'qry_end':int(cols[7]), \
                    'hit_id':hit_id, 'hit_start':int(cols[8]), 'hit_end':int(cols[9]), \
                    'pct_id':float(cols[10]), 'strand':cols[17], 'eval':float(cols[19]) }

        if algorithm is None:
            algorithm = cols[3]

        if qry_id != current_qry_id or hit_id != current_hit_id:
            # this is a new chain, export the previous one
            if current_qry_id is not None:
                export_match( current_match_parts, ofh, algorithm, args.perc_identity_cutoff )
                
            current_match_parts = list()
            current_qry_id = qry_id
            current_hit_id = hit_id

        current_match_parts.append( segment )
    
    ## make sure to do the last one
    if current_qry_id is not None:
        export_match( current_match_parts, ofh, algorithm, args.perc_identity_cutoff )

    



def export_match( segments, out, source, perc_id_cutoff ):
    contig_min = None
    contig_max = None
    contig_id  = None
    hit_id     = None
    hit_min    = None
    hit_max    = None
    segment_strand = '+'

    for segment in segments:
        segment_min = min( segment['qry_start'], segment['qry_end'] )
        segment_max = max( segment['qry_start'], segment['qry_end'] )
        segment_hit_min = min( segment['hit_start'], segment['hit_end'] )
        segment_hit_max = max( segment['hit_start'], segment['hit_end'] )

        if segment['strand'] == 'Minus':
            segment_strand = '-'

        if contig_min is None or segment_min < contig_min:
            contig_min = segment_min
            hit_min    = segment_hit_min
        
        if contig_max is None or segment_max > contig_max:
            contig_max = segment_max
            hit_max = segment_hit_max

        if contig_id is None:
            contig_id = segment['qry_id']

        if hit_id is None:
            hit_id = segment['hit_id']

    ## does this score above any user-defined cutoffs.
    if perc_id_cutoff is not None:
        pct_id = global_pct_id( segments )

        if pct_id < perc_id_cutoff:
            return False

    
    ## write the match feature
    match_id = "match.{0}".format(next_ids['match'])
    next_ids['match'] += 1
    data_column = "ID={0};Target={1} {2} {3}".format(match_id, hit_id, hit_min, hit_max)
    out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
            contig_id, source, 'match', contig_min, contig_max, '.', segment_strand, '.', data_column
            ) )

    ## write the match_parts
    for segment in sorted(segments, key=itemgetter('qry_start', 'qry_end')):
        match_part_id = "match_part.{0}".format(next_ids['match_part'])
        next_ids['match_part'] += 1
        mp_start = min( segment['qry_start'], segment['qry_end'] )
        mp_hit_start = min(segment['hit_start'], segment['hit_end'])
        mp_end   = max( segment['qry_start'], segment['qry_end'] )
        mp_hit_end = max( segment['hit_start'], segment['hit_end'] )

        data_column = "ID={0};Parent={1};Target={2} {3} {4}".format(match_part_id, match_id, hit_id, \
                mp_hit_start, mp_hit_end)
        out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
                contig_id, source, 'match_part', mp_start, mp_end, \
                '.', segment_strand, '.', data_column 
                ) )

    return True

if __name__ == '__main__':
    main()





















