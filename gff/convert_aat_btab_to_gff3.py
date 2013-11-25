#!/usr/bin/env python3

"""
Converts the btab output of AAT (nap/gap2) to GFF3 format with option-dependent modeling.

Example input:
jcf7180000787896        Aug 28 2013     16242   /usr/local/packages/aat/nap     /usr/local/projects/mucormycosis/protein_alignments/fungi_jgi/fungi_jgi.faa     jgi|Rhior3|10928|RO3G_01207     1       103     1162    1198    25.000000       41.666667       15      1       1               -1      Plus    1484
jcf7180000787896        Aug 28 2013     16242   /usr/local/packages/aat/nap     /usr/local/projects/mucormycosis/protein_alignments/fungi_jgi/fungi_jgi.faa     jgi|Rhior3|10928|RO3G_01207     870     1729    1198    1484    91.958042       93.356643       1384    1       2               -1      Plus    1484

The output here depends on the mode requested.  If --export_mode=model (default), then it is:

jcf7180000787896        nap     gene    1       1729    .       +       .       ID=gene.000001
jcf7180000787896        nap     mRNA    1       1729    .       +       .       ID=mRNA.000001;Parent=gene.000001
jcf7180000787896        nap     exon    1       103     .       +       .       ID=exon.000001;Parent=mRNA.000001
jcf7180000787896        nap     CDS     1       103     .       +       0       ID=exon.000001;Parent=mRNA.000001
jcf7180000787896        nap     exon    870     1729    .       +       .       ID=exon.000002;Parent=mRNA.000001
jcf7180000787896        nap     CDS     870     1729    .       +       0       ID=exon.000002;Parent=mRNA.000010

On the other hand, if you specify --export_mode=match you'll get:


Note that chain membership is (intentionally) lost in this mode, and all segments appear as
individual, distinct matches.

"""

import argparse
import os
import re

next_ids = {'gene':1, 'mRNA':1, 'exon':1, 'CDS':1, 'match':1 }
MATCH_TERM = 'nucleotide_to_protein_match'

def main():
    parser = argparse.ArgumentParser( description='AAT (nap/gap2) converter to GFF3 format')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to parse' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-p', '--name_prefix', type=str, required=False, help='So that resulting GFF3 files from multiple runs can be merged, you can supply a prefix for the IDs generated here' )
    parser.add_argument('-d', '--perc_identity_cutoff', type=float, required=False, help='Filters on the percent identity of over the length of the alignment' )
    parser.add_argument('-s', '--perc_similarity_cutoff', type=float, required=False, help='Filters on the percent similarity of over the length of the alignment' )
    parser.add_argument('-m', '--max_intron_cutoff', type=int, required=False, help='Excludes alignment chains which propose an intron greater than this value' )
    parser.add_argument('-e', '--export_mode', type=str, required=False, default='model', help='Controls whether chains are modeled as genes or match regions' )
    
    args = parser.parse_args()
    algorithm = None
    total_chains = 0
    chains_exported = 0

    ofh = open(args.output_file, 'w')
    ofh.write("##gff-version 3\n")
    
    current_chain_number = 1
    gene_segments = []

    for line in open(args.input_file, 'r'):
        cols = line.split("\t")
        chain_num = int(cols[13]);
        
        segment = { 'contig_id':cols[0], 'contig_start':int(cols[6]), 'contig_end':int(cols[7]), \
                    'hit_id':cols[5], 'hit_start':int(cols[8]), 'hit_end':int(cols[9]), \
                    'pct_id':float(cols[10]), 'pct_sim':float(cols[11]), \
                    'strand':cols[17], 'product':cols[15], 'score':cols[18] }

        ## aat product fields are a bit of mess, so lets tease a few parts out:
        #   | organism=Neospora_caninum_Liverpool | product=Iron regulatory protein-like protein, related | location=FR823391:1231501-1237765(-) | length=986 | sequence_SO=chromosome | SO=protein_coding
        #   | organism=Toxoplasma_gondii_ME49 | product=aconitate hydratase ACN/IRP | location=TGME49_chrX:1400716-1407887(-) | length=1055 | sequence_SO=chromosome | SO=protein_coding
        m = re.search("product=(.+?) \|", segment['product'])
        if m:
            hit_product = m.group(1)
        else:
            hit_product = None
        
        if algorithm is None:
            algorithm = os.path.basename(cols[3])

        if args.export_mode == 'match':
            segment_min = min( segment['contig_start'], segment['contig_end'] )
            segment_max = max( segment['contig_start'], segment['contig_end'] )
            hit_min = min( segment['hit_start'], segment['hit_end'] )
            hit_max = max( segment['hit_start'], segment['hit_end'] )
            
            if segment['strand'] == 'Minus':
                segment_strand = '-'
            else:
                segment_strand = '+'

            match_id = get_next_id('match', args.name_prefix)
            data_column = "ID={0};Target={1} {2} {3};Name={4}".format(match_id, segment['hit_id'], hit_min, hit_max, hit_product)
            ofh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
                    segment['contig_id'], algorithm, MATCH_TERM, segment_min, segment_max, '.', segment_strand, \
                    '.', data_column ))
        else:
            if chain_num != current_chain_number:
                chain_pct_id = global_pct_id( gene_segments )

                # this is a new chain
                total_chains += 1
                exported = export_gene( gene_segments, ofh, current_chain_number, algorithm, args.name_prefix, args.perc_identity_cutoff, \
                                        args.perc_similarity_cutoff, args.max_intron_cutoff )
                if exported is True:
                    chains_exported += 1

                current_chain_number = chain_num
                gene_segments = []

            gene_segments.append( segment )
    
    ## make sure to do the last one
    if args.export_mode == 'match':
        pass
    else:
        total_chains += 1
        exported = export_gene( gene_segments, ofh, current_chain_number, algorithm, args.name_prefix, args.perc_identity_cutoff, \
                                args.perc_similarity_cutoff, args.max_intron_cutoff )
        if exported is True:
            chains_exported += 1

        print("\nTotal alignment chains found: {0}".format(total_chains) )
        print("Total chains exported (after cutoffs applied): {0}\n".format(chains_exported) )
    

def get_next_id(type, prefix):
    if prefix == None:
        id = type + ".{0:06d}".format( next_ids[type] )
    else:
        id = prefix + "." + type + ".{0:06d}".format( next_ids[type] )
    
    next_ids[type] += 1

    return id


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
        segment_length = abs(segment['contig_start'] - segment['contig_end']) + 1
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


def export_match():
    pass


def export_gene( segments, out, chn_num, source, prefix, perc_id_cutoff, perc_sim_cutoff, max_intron_cutoff ):
    contig_min = None
    contig_max = None
    contig_id  = None
    hit_id     = None
    hit_min    = None
    hit_max    = None
    segment_strand = '+'

    segment_last_max = None

    # these are always returned in the file in order, so there's no need to sort segments
    for segment in segments:
        segment_min = min( segment['contig_start'], segment['contig_end'] )
        segment_max = max( segment['contig_start'], segment['contig_end'] )
        segment_hit_min = min( segment['hit_start'], segment['hit_end'] )
        segment_hit_max = max( segment['hit_start'], segment['hit_end'] )

        if segment_last_max is not None:
            intron_size = segment_min - segment_last_max - 1
            if max_intron_cutoff is not None and intron_size > max_intron_cutoff:
                return False

        segment_last_max = segment_max
        
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

    ## does this score above any user-defined cutoffs.
    if perc_id_cutoff is not None:
        pct_id = global_pct_id( segments )

        if pct_id < perc_id_cutoff:
            return False

    if perc_sim_cutoff is not None:
        pct_sim = global_pct_sim( segments )

        if pct_sim < perc_sim_cutoff:
            return False
    
    #print("Exported chain: ({0}) min:({1}) max({2})".format(chn_num, contig_min, contig_max))
    
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
        data_column = "ID={0};Parent={1};Target={2} {3} {4}".format(exon_id, mRNA_id, hit_id, \
                exon_hit_start, exon_hit_end)
        out.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
                contig_id, source, 'CDS', exon_start, exon_end, \
                '.', segment_strand, 0, data_column 
                ) )

    return True

if __name__ == '__main__':
    main()





















