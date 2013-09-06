#!/usr/bin/env python3.2

"""
Converts the btab output of AAT (nap/gap2) to GFF3 format with proper gene models.

Example input:
jcf7180000787896        Aug 28 2013     16242   /usr/local/packages/aat/nap     /usr/local/projects/mucormycosis/protein_alignments/fungi_jgi/fungi_jgi.faa     jgi|Rhior3|10928|RO3G_01207     1       103     1162    1198    25.000000       41.666667       15      1       1               -1      Plus    1484
jcf7180000787896        Aug 28 2013     16242   /usr/local/packages/aat/nap     /usr/local/projects/mucormycosis/protein_alignments/fungi_jgi/fungi_jgi.faa     jgi|Rhior3|10928|RO3G_01207     870     1729    1198    1484    91.958042       93.356643       1384    1       2               -1      Plus    1484

Example output:
jcf7180000787896        nap     gene    1       1729    .       +       .       ID=gene.000001
jcf7180000787896        nap     mRNA    1       1729    .       +       .       ID=mRNA.000001;Parent=gene.000001
jcf7180000787896        nap     exon    1       103     .       +       .       ID=exon.000001;Parent=mRNA.000001
jcf7180000787896        nap     CDS     1       103     .       +       0       ID=exon.000001;Parent=mRNA.000001
jcf7180000787896        nap     exon    870     1729    .       +       .       ID=exon.000002;Parent=mRNA.000001
jcf7180000787896        nap     CDS     870     1729    .       +       0       ID=exon.000002;Parent=mRNA.000010

Future development possibilities:
- Add an output mode for exporting alignments as 'nucleotide_to_protein_match' SO features
  rather than gene models.
"""

import argparse
import os

next_ids = {'gene':1, 'mRNA':1, 'exon':1, 'CDS':1 }

def main():
    parser = argparse.ArgumentParser( description='AAT (nap/gap2) converter to GFF3 format')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to parse' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-p', '--name_prefix', type=str, required=False, help='So that resulting GFF3 files from multiple runs can be merged, you can supply a prefix for the IDs generated here' )
    args = parser.parse_args()
    algorithm = None

    ofh = open(args.output_file, 'w')
    current_chain_number = 1
    gene_segments = []

    for line in open(args.input_file, 'r'):
        cols = line.split("\t")
        chain_num = int(cols[13]);
        
        segment = { 'contig_id':cols[0], 'contig_start':int(cols[6]), 'contig_end':int(cols[7]), \
                    'hit_id':cols[5], 'hit_start':int(cols[8]), 'hit_end':int(cols[9]), \
                    'strand':cols[17], 'product':cols[15], 'score':cols[18] }

        if algorithm is None:
            algorithm = os.path.basename(cols[3])

        if chain_num != current_chain_number:
            # this is a new chain
            export_gene( gene_segments, ofh, current_chain_number, algorithm, args.name_prefix )
            current_chain_number = chain_num
            gene_segments = []

        gene_segments.append( segment )
    
    ## make sure to do the last one
    export_gene( gene_segments, ofh, current_chain_number, algorithm, args.name_prefix )
        

def get_next_id(type, prefix):
    if prefix == None:
        id = type + ".{0:06d}".format( next_ids[type] )
    else:
        id = prefix + type + ".{0:06d}".format( next_ids[type] )
    
    next_ids[type] += 1

    return id


def export_gene( segments, out, chn_num, source, prefix ):
    contig_min = None
    contig_max = None
    contig_id  = None
    hit_id     = None
    hit_min    = None
    hit_max    = None
    segment_strand = '+'

    for segment in segments:
        segment_min = min( segment['contig_start'], segment['contig_end'] )
        segment_max = max( segment['contig_start'], segment['contig_end'] )
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
            contig_id = segment['contig_id']

        if hit_id is None:
            hit_id = segment['hit_id']

    print("Exported chain: ({0}) min:({1}) max({2})".format(chn_num, contig_min, contig_max))
    
    ## write the gene feature
    gene_id = get_next_id('gene', prefix)
    data_column = "ID={0};Target={1} {2} {3}".format(gene_id, hit_id, hit_min, hit_max)
    out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
            contig_id, source, 'gene', contig_min, contig_max, '.', segment_strand, '.', data_column
            ) )

    ## write the mRNA
    mRNA_id = get_next_id('mRNA', prefix)
    data_column = "ID={0};Parent={1};Target={2} {3} {4}".format(mRNA_id, gene_id, hit_id, hit_min, hit_max)
    out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
            contig_id, source, 'mRNA', contig_min, contig_max, '.', segment_strand, '.', data_column
            ) )

    ## write the exons and CDS
    for segment in segments:
        exon_id = get_next_id('exon', prefix)
        exon_start = min( segment['contig_start'], segment['contig_end'] )
        exon_hit_start = min(segment['hit_start'], segment['hit_end'])
        exon_end   = max( segment['contig_start'], segment['contig_end'] )
        exon_hit_end = max( segment['hit_start'], segment['hit_end'] )

        data_column = "ID={0};Parent={1};Target={2} {3} {4}".format(exon_id, mRNA_id, hit_id, \
                exon_hit_start, exon_hit_end)
        out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
                contig_id, source, 'exon', exon_start, exon_end, \
                '.', segment_strand, '.', data_column 
                ) )

        CDS_id = get_next_id('CDS', prefix)
        data_column = "ID={0};Parent={1}Target={2} {3} {4}".format(exon_id, mRNA_id, hit_id, \
                exon_hit_start, exon_hit_end)
        out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
                contig_id, source, 'CDS', exon_start, exon_end, \
                '.', segment_strand, 0, data_column 
                ) )

if __name__ == '__main__':
    main()





















