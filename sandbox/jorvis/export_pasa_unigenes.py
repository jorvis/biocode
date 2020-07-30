#!/usr/bin/env python3

"""

After running PASA (https://github.com/PASApipeline/PASApipeline) we want to create a 'unigene'
set, which includes the transcript from each PASA cluster with the longest ORF.  This approach
was published in ??

Input expectations:

In the input_pasta_fasta file the IDs are expected to look like this:

   >asmbl_1
   >asmbl_2
   >asmbl_3

Within the GTF, the last column should look like:

   gene_id "PASA_cluster_1"; transcript_id "align_id:907554|asmbl_1";
   gene_id "PASA_cluster_2"; transcript_id "align_id:907555|asmbl_2";
   gene_id "PASA_cluster_19"; transcript_id "align_id:907573|asmbl_20";
   gene_id "PASA_cluster_19"; transcript_id "align_id:907574|asmbl_21";

Since transdecoder is run on the PASA FASTA output, its IDs should match those PASA provided.

Assumes transdecoder output looks like this ('longest_orfs.cds'):

   >asmbl_6.p1 type:5prime_partial len:126 asmbl_6:2055-1678(-)
   >asmbl_6.p2 type:complete len:104 asmbl_6:374-685(+)

"""

import argparse
import os
import re
import sys

from biocode import utils

def main():
    parser = argparse.ArgumentParser( description='Export the transcript from each PASA cluster with the longest ORF')
    parser.add_argument('-if', '--input_pasa_fasta', type=str, required=True, help='Path to PASAs predicted FASTA file' )
    parser.add_argument('-itsc', '--input_pasa_transcript_size_cutoff', type=str, required=False, default=300, help='Nucleotide size cutoff of input transcripts' )
    parser.add_argument('-itc', '--input_pasa_tpm_cutoff', type=str, required=False, default=0.1, help='TPM cutoff of input transcripts' )
    parser.add_argument('-ig', '--input_pasa_gtf', type=str, required=True, help='Path to PASAs predicted GFT file' )
    parser.add_argument('-it', '--input_transdecoder', type=str, required=True, help='Path to the longest_orfs.cds from transdecoder' )
    parser.add_argument('-is', '--input_salmon', type=str, required=True, help='Path to quant.sf file from Salmon' )
    parser.add_argument('-ol', '--output_file_longest_orf', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-oa', '--output_file_most_abundant', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    seqs = utils.fasta_dict_from_file(args.input_pasa_fasta)
    pasa_clusters = load_pasa_clusters(args.input_pasa_gtf)
    orf_lengths = longest_orf_lengths(args.input_transdecoder)
    load_abundances(seqs, args.input_salmon)

    pasa_cluster_count = len(pasa_clusters)
    pasa_transcript_count = len(seqs)
    size_filtered_pasa_transcript_count = 0
    size_filtered_pasa_transcript_count_with_cds = 0
    size_filtered_pasa_cluster_count = 0
    tpm_filtered_pasa_transcript_count = 0
    tpm_filtered_pasa_cluster_count = 0

    print("Stats:\n")
    print("Initial PASA cluster count: {0}".format(pasa_cluster_count))
    print("Initial PASA transcript count: {0}".format(pasa_transcript_count))
    
    seqs, pasa_clusters = apply_transcript_size_cutoff(seqs, pasa_clusters, args.input_pasa_transcript_size_cutoff)
    seqs, pasa_clusters = apply_abundance_filter(seqs, pasa_clusters, args.input_pasa_tpm_cutoff)
    seqs, pasa_clusters = apply_orf_filter(seqs, pasa_clusters, orf_lengths, 100)

def write_transcripts(seqs, filename):
    with open(filename, 'wt') as fh:
        for asmbl_id in seqs:
            fh.write(">{0} {1}\n{2}\n".format(asmbl_id, seqs[asmbl_id]['h'], 
                                              utils.wrapped_fasta(''.join(seqs[asmbl_id]['s']))))
            
def apply_orf_filter(seqs, clusters, orf_lengths, orf_length_cutoff):
    clusters_to_delete = list()
    
    for cluster_id in clusters:
        asmbl_ids_to_delete = list()
        
        for asmbl_id in clusters[cluster_id]:
            if asmbl_id not in orf_lengths or orf_lengths[asmbl_id] < orf_length_cutoff:
                del seqs[asmbl_id]
                asmbl_ids_to_delete.append(asmbl_id)

        for id in asmbl_ids_to_delete:
            del clusters[cluster_id][id]

            if len(clusters[cluster_id]) == 0:
                clusters_to_delete.append(cluster_id)

    for id in clusters_to_delete:
        del clusters[id]

    print("After ORF length filtering:")
    print("\tTranscripts: {0}".format(len(seqs)))
    print("\tPASA clusters: {0}".format(len(clusters)))

    write_transcripts(seqs, 'pasa.orf_filtered.fasta')
    
    return seqs, clusters
    
def apply_abundance_filter(seqs, clusters, tpm_cutoff):
    clusters_to_delete = list()
    
    for cluster_id in clusters:
        asmbl_ids_to_delete = list()
        
        for asmbl_id in clusters[cluster_id]:
            if 'tpm' in seqs[asmbl_id]:
                if seqs[asmbl_id]['tpm'] < tpm_cutoff:
                    del seqs[asmbl_id]
                    asmbl_ids_to_delete.append(asmbl_id)

        for id in asmbl_ids_to_delete:
            del clusters[cluster_id][id]

            if len(clusters[cluster_id]) == 0:
                clusters_to_delete.append(cluster_id)

    for id in clusters_to_delete:
        del clusters[id]

    print("After TPM filtering:")
    print("\tTranscripts: {0}".format(len(seqs)))
    print("\tPASA clusters: {0}".format(len(clusters)))

    write_transcripts(seqs, 'pasa.tpm_filtered.fasta')
                    
    return seqs, clusters

def apply_transcript_size_cutoff(seqs, clusters, size_cutoff):
    clusters_to_delete = list()
    
    for cluster_id in clusters:
        asmbl_ids_to_delete = list()
        
        for asmbl_id in clusters[cluster_id]:
            if len(seqs[asmbl_id]['s']) < size_cutoff:
                del seqs[asmbl_id]
                asmbl_ids_to_delete.append(asmbl_id)

        for id in asmbl_ids_to_delete:
            del clusters[cluster_id][id]

            if len(clusters[cluster_id]) == 0:
                clusters_to_delete.append(cluster_id)

    for id in clusters_to_delete:
        del clusters[id]

    print("After size filtering:")
    print("\tTranscripts: {0}".format(len(seqs)))
    print("\tPASA clusters: {0}".format(len(clusters)))

    write_transcripts(seqs, 'pasa.size_filtered.fasta')
                    
    return seqs, clusters
    
def load_abundances(seqs, salmon_file):
    line_num = 0
    for line in open(salmon_file):
        line_num += 1
        if line_num == 1:
            continue
        
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) < 5:
            continue

        seq_id = cols[0]
        tpm = float(cols[3])

        if seq_id not in seqs:
            raise Exception("ERROR: Got seq ID ({0}) in salmon input that wasn't in FASTA input".format(seq_id))

        seqs[seq_id]['tpm'] = tpm

    
def load_pasa_clusters(pasa_gtf):
    clusters = dict()
    for line in open(pasa_gtf):
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) < 9:
            continue

        if cols[2] != 'transcript':
            continue

        m = re.match("gene_id \"(.+?)\"; transcript_id \"align_id:\d+\|(.+?)\";", cols[8])

        if m:
            cluster_id = m.group(1)
            asmbl_id = m.group(2)

            if cluster_id not in clusters:
                clusters[cluster_id] = dict()

            # initialize all their longest-ORF lengths to zero
            clusters[cluster_id][asmbl_id] = 0

    return clusters
        
def longest_orf_lengths(transdecoder_file):
    orf_lengths = dict()
    
    for line in open (transdecoder_file):
        line = line.rstrip()
        m = re.match(">.+? type:\S+ len:(\d+) (asmbl_\d+)\:.+", line)
        if m:
            orf_length = int(m.group(1))
            asmbl_id = m.group(2)

            if asmbl_id in orf_lengths:
                if orf_length > orf_lengths[asmbl_id]:
                    orf_lengths[asmbl_id] = orf_length
            else:
                orf_lengths[asmbl_id] = orf_length

    return orf_lengths

if __name__ == '__main__':
    main()







