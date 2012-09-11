#!/usr/bin/env python3.2

import argparse
import os

next_ids = {'gene':1, 'mRNA':1, 'exon':1, 'CDS':1 }

def main():
    parser = argparse.ArgumentParser( description='AAT (nap/gap2) converter to GFF3 format')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to parse' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()
    algorithm = None

    ofh = open(args.output_file, 'w')
    current_chain_number = 1
    gene_segments = []

    for line in open(args.input_file, 'r'):
        cols = line.split("\t")
        chain_num = int(cols[13]);
        
        segment = { 'contig_id':cols[0], 'contig_start':int(cols[6]), 'contig_end':int(cols[7]), \
                    'hit_start':int(cols[8]), 'hit_end':int(cols[9]), \
                    'strand':cols[17], 'product':cols[15], 'score':cols[18] }

        if algorithm is None:
            algorithm = os.path.basename(cols[3])

        if chain_num != current_chain_number:
            # this is a new chain
            export_gene( gene_segments, ofh, current_chain_number, algorithm )
            current_chain_number = chain_num
            gene_segments = []

        gene_segments.append( segment )
    
    ## make sure to do the last one
    export_gene( gene_segments, ofh, current_chain_number, algorithm )
        

def get_next_id(type):
    id = type + "{0:06d}".format( next_ids[type] )
    next_ids[type] += 1

    return id


def export_gene( segments, out, chn_num, source ):
    contig_min = None
    contig_max = None
    contig_id  = None
    segment_strand = '+'

    for segment in segments:
        segment_min = min( segment['contig_start'], segment['contig_end'] )
        segment_max = max( segment['contig_start'], segment['contig_end'] )

        if segment['strand'] == 'Minus':
            segment_strand = '-'

        if contig_min is None or segment_min < contig_min:
            contig_min = segment_min
        
        if contig_max is None or segment_max > contig_max:
            contig_max = segment_max

        if contig_id is None:
            contig_id = segment['contig_id']

    print("Exported chain: ({0}) min:({1}) max({2})".format(chn_num, contig_min, contig_max))
    
    ## write the gene feature
    gene_id = get_next_id('gene')
    data_column = "ID={0}".format(gene_id)
    out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
            contig_id, source, 'gene', contig_min, contig_max, '.', segment_strand, '.', data_column 
            ) )

    ## write the mRNA
    mRNA_id = get_next_id('mRNA')
    data_column = "ID={0};Parent={1}".format(mRNA_id, gene_id)
    out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
            contig_id, source, 'mRNA', contig_min, contig_max, '.', segment_strand, '.', data_column 
            ) )

    ## write the exons and CDS
    for segment in segments:
        exon_id = get_next_id('exon')
        exon_start = min( segment['contig_start'], segment['contig_end'] )
        exon_end   = max( segment['contig_start'], segment['contig_end'] )
        data_column = "ID={0};Parent={1}".format(exon_id, mRNA_id)
        out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
                contig_id, source, 'exon', exon_start, exon_end, \
                '.', segment_strand, '.', data_column 
                ) )

        CDS_id = get_next_id('CDS')
        data_column = "ID={0};Parent={1}".format(exon_id, mRNA_id)
        out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
                contig_id, source, 'CDS', exon_start, exon_end, \
                '.', segment_strand, 0, data_column 
                ) )

if __name__ == '__main__':
    main()





















