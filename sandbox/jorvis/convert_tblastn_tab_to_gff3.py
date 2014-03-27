#!/usr/bin/env python3

"""
Converts the tabular output of NCBI's tblastn output when using the -m 9
or -m 9 options to GFF3 with match and match_part features.

This script assumes you're using the protein file as the query and the
DNA one as the subject.

Note, this is NOT for 'btab' output, which is the tabular output from
and Ergatis pipeline.  There's another script for that.

INPUT FIELDS:

1.  Query id
2.  Subject id
3.  % identity
4.  alignment length
5.  mismatches
6.  gap openings
7.  q. start
8.  q. end
9.  s. start
10. s. end
11. e-value
12. bit score

EXAMPLE INPUT

# TBLASTN 2.2.21 [Jun-14-2009]
# Query: CMU_015380 | organism=Cryptosporidium_muris_RN66 | product=zinc knuckle family protein | location=DS989730:711568-713079(-) | length=503 | sequence_SO=supercontig | SO=protein_coding
# Database: B_microti_4_chrom_plus_mito.fasta
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
CMU_015380      ChromosomeIII_BmicrotiR1        72.06   68      19      0       164     231     1637106 1637309 2e-51    111
CMU_015380      ChromosomeIII_BmicrotiR1        42.19   128     47      2       232     332     1637333 1637689 2e-51    108
CMU_015380      ChromosomeIII_BmicrotiR1        53.28   122     57      0       42      163     1636378 1636743 9e-31    129
CMU_015380      ChromosomeIV_BmicrotiR1 25.76   66      40      2       265     321     697665  697477  1.6     29.3
CMU_015380      ChromosomeII_BmicrotiR1 24.73   93      68      3       151     241     1393304 1393558 3.8     28.1

(Note that any lines starting with '#' are ignored.  These are present in the -m 9 output only.)

EXAMPLE OUTPUT




"""

import argparse
import os
from operator import itemgetter

next_ids = {'match':1, 'match_part':1 }

def main():
    parser = argparse.ArgumentParser( description='NCBI-BLAST (tblastn/tab) converter to GFF3 format')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to parse' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-p', '--perc_identity_cutoff', type=float, required=False, help='Filters on the perc identity of each HSP' )
    parser.add_argument('-cf', '--custom_filter', action='store_true', help='Applies a custom set of filters defined by a certain PI.  These are not the droids you are looking for.')
    
    args = parser.parse_args()
    algorithm = 'TBLASTN'
    current_qry_id = None
    current_sbj_id = None
    current_match_parts = list()

    ofh = open(args.output_file, 'w')
    ofh.write("##gff-version 3\n")
    
    for line in open(args.input_file, 'r'):
        if line.startswith('#'):
            continue
        
        cols = line.split("\t")
        qry_id = cols[1]
        sbj_id = cols[0]
        strand = '+'

        if int(cols[9]) < int(cols[8]):
            strand = '-'
        
        segment = { 'qry_id':qry_id,  'qry_start':int(cols[8]), 'qry_end':int(cols[9]), \
                    'hit_id':sbj_id, 'hit_start':int(cols[6]), 'hit_end':int(cols[7]), \
                    'pct_id':float(cols[2]), 'strand':strand, 'eval':float(cols[10]) }

        if qry_id != current_qry_id or sbj_id != current_sbj_id:
            # this is a new chain, export the previous one
            if current_qry_id is not None:
                export_match( current_match_parts, ofh, algorithm, args.perc_identity_cutoff, args.custom_filter )
                
            current_match_parts = list()
            current_qry_id = qry_id
            current_sbj_id = sbj_id

        current_match_parts.append( segment )
    
    ## make sure to do the last one
    if current_qry_id is not None:
        export_match( current_match_parts, ofh, algorithm, args.perc_identity_cutoff, args.custom_filter )

    



def export_match( segments, out, source, perc_id_cutoff, custom_filter_requested ):
    contig_min = None
    contig_max = None
    contig_id  = None
    hit_id     = None
    hit_min    = None
    hit_max    = None

    for segment in segments:
        segment_min = min( segment['qry_start'], segment['qry_end'] )
        segment_max = max( segment['qry_start'], segment['qry_end'] )
        segment_hit_min = min( segment['hit_start'], segment['hit_end'] )
        segment_hit_max = max( segment['hit_start'], segment['hit_end'] )

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

    match_length = contig_max - contig_min
    pct_id = global_pct_id( segments )

    if custom_filter_requested == True:
        criteria_met = False
        if (match_length <= 20 and pct_id >= 80) or \
           (match_length > 20 and match_length <= 50 and pct_id >= 60) or \
           (match_length > 50 and pct_id > 40):
            criteria_met = True

        if criteria_met == False:
            return False

    ## does this score above any user-defined cutoffs.
    if perc_id_cutoff is not None:
        if pct_id < perc_id_cutoff:
            return False

    
    ## write the match feature
    match_id = "match.{0}".format(next_ids['match'])
    next_ids['match'] += 1
    data_column = "ID={0};Target={1} {2} {3}".format(match_id, hit_id, hit_min, hit_max)
    out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
            contig_id, source, 'match', contig_min, contig_max, '.', segments[0]['strand'], '.', data_column
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
                '.', segment['strand'], '.', data_column 
                ) )

    return True


def global_pct_id( segments ):
    """
    Calculated like this:
    
    10bp @ 50% id = 5 matching residues, 10 total residues
    10bp @ 80% id = 8 matching residues, 10 total residues
                   13 matching residues, 20 total residues
                   ---------------------------------------
                   13 / 20 * 100 = 65%

    """
    match_length = 0
    identical_residues = 0
    
    for segment in segments:
        segment_length = abs(segment['qry_start'] - segment['qry_end']) + 1
        match_length += segment_length
        matched_residues = segment_length * (segment['pct_id'] / 100)
        identical_residues += matched_residues

    return (identical_residues / match_length) * 100


def global_pct_sim( segments ):
    """
    Calculated as in global_pct_id(), but with the percent similarity field
    """
    match_length = 0
    similar_residues = 0
    
    for segment in segments:
        segment_length = abs(segment['contig_start'] - segment['contig_end']) + 1
        match_length += segment_length
        matched_residues = segment_length * (segment['pct_sim'] / 100)
        similar_residues += matched_residues

    return (similar_residues / match_length) * 100


if __name__ == '__main__':
    main()





















