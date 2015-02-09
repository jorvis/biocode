#!/usr/bin/env python3

"""

We have an array of files containing the IDs of transcripts whose orientations were
ambiguous based on ratios of read mapping direction.  The point of this script is to
find which IDs fell into this category across different tissue sites.

Input:

The input directory should contain files like:

A10_gills_reorientation.midrange.0.05-0.95.id.list
A10_gills_reorientation.midrange.0.125-0.875.id.list
A10_gills_reorientation.midrange.0.25-0.75.id.list
A10_gills_reorientation.midrange.0.475-0.525.id.list
A1_CNS_reorientation.midrange.0.05-0.95.id.list
A1_CNS_reorientation.midrange.0.125-0.875.id.list
A1_CNS_reorientation.midrange.0.25-0.75.id.list
A1_CNS_reorientation.midrange.0.475-0.525.id.list
A3_digestive_reorientation.midrange.0.05-0.95.id.list
A3_digestive_reorientation.midrange.0.125-0.875.id.list
A3_digestive_reorientation.midrange.0.25-0.75.id.list
A3_digestive_reorientation.midrange.0.475-0.525.id.list
A7_heart_reorientation.midrange.0.05-0.95.id.list
A7_heart_reorientation.midrange.0.125-0.875.id.list
A7_heart_reorientation.midrange.0.25-0.75.id.list
A7_heart_reorientation.midrange.0.475-0.525.id.list
all.prepended.sorted.uc
x

PROCESS:
- Read through cluster file store IDs of all members
    IL100038439_comp81521_c0_seq4
    IL100038709_IL100038708_comp128043_c3_seq141

- Look up this ID in its list file to see if its in or out of the midrange
    IL100038435 -> A3
    A3 ->
        A3_reorientation.midrange.0.05-0.95.id.list
        A3_reorientation.midrange.0.125-0.875.id.list
        A3_reorientation.midrange.0.25-0.75.id.list
        A3_reorientation.midrange.0.475-0.525.id.list

counts[cluster_id][$range_label] = N

"""

import argparse
import os
import re

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-id', '--input_directory', type=str, required=True, help='Directory where all the input files are found' )
    parser.add_argument('-c', '--cluster_file', type=str, required=True, help='' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    ranges = {'A': [0.05, 0.95], 'B': [0.125, 0.875], 'C': [0.25, 0.75], 'D': [0.475, 0.525]}
    targets = ['A1', 'A3', 'A7', 'A10']
    
    cluster_members = dict()

    ofh = open(args.output_file, 'wt')

    ## cluster file IDs:  IL100038709_IL100038708_comp128043_c3_seq141
    ## mapping file: A3_reorientation.midrange.0.05-0.95.id.list : comp1893_c0_seq1
    sample_labels = {
        'IL100038432': 'A1',
        'IL100038432_IL100038433': 'A1',
        'IL100038434': 'A2',
        'IL100038435': 'A3',
        'IL100038436': 'A5',
        'IL100038438': 'A7',
        'IL100038439': 'A8',
        'IL100038440': 'A10',
        'IL100038441': 'A10',
        'IL100038709': 'A1',
        'IL100038709_IL100038708': 'A1',
        'IL100038721': 'B6'
    }

    # structure like
    # dict['A3']['A'] = list()
    midrange_members = {}
    for sample_label in ['A1', 'A3', 'A7', 'A10']:
        if sample_label not in midrange_members:
            midrange_members[sample_label] = dict()

        for range_label in ranges:
            range_min = ranges[range_label][0]
            range_max = ranges[range_label][1]
            range_string = "{0}-{1}".format(range_min, range_max)

            if range_label not in midrange_members[sample_label]:
                midrange_members[sample_label][range_label] = list()

            for line in open("{0}_reorientation.midrange.{1}.id.list".format(sample_label, range_string)):
                line = line.rstrip()
                midrange_members[sample_label][range_label].append(line)

    # Load up the cluster membership
    for line in open(args.cluster_file):
        if line.startswith('#'): continue
        cols = line.split()
        cluster_num = cols[1]

        # have we seen this cluster number?
        if cluster_num not in cluster_members:
            cluster_members[cluster_num] = list()

        # We only want S (seed) and H (hit) lines
        if cols[0] == 'H' or cols[1] == 'S':
            m = re.search("^(\S+)", cols[8])
            transcript_id = None
            if m:
                transcript_id = m.group(1)
                cluster_members[cluster_num].append(transcript_id)

    # Now check the membership of each cluster, then print the counts
    for cluster_id in cluster_members:
        ofh.write("{0}\t".format(cluster_id))
        
        for range_label in ranges:
            transcript_count = 0

            for transcript_id in cluster_members[cluster_id]:
                # need to get the sample label from the prefix of the transcript ID
                m = re.match("^(IL.+?)_c", transcript_id)
                if m:
                    sample_id = m.group(1)
                else:
                    raise Exception("ERROR: transcript ID didn't follow expected convention: {0}".format(transcript_id))

                sample_name = sample_labels[sample_id]

                if sample_name not in targets:
                    continue

                if transcript_id not in midrange_members[sample_name][range_label]:
                    transcript_count += 1

            ofh.write("{0}: {1}\t".format(range_label, transcript_count))

        ofh.write("\n")

                    

                


        

if __name__ == '__main__':
    main()







