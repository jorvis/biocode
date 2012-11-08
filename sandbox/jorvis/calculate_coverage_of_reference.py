#!/usr/bin/env python3.2

import argparse
import math
import os
import re

from operator import itemgetter

"""
========
OVERVIEW
========

This reads nucmer output and an annotated GFF3 to provide a lot of statistics and lists
related to coverage of a reference genome.

=====
INPUT
=====

Options
-------
-c, --coords_file         required
    Path to a nucmer coords file with non-overlapping results (requires -l -r -T options
    of show-coords)

-o, --output_prefix       required
    This script creates several output files.  The string specified here will be the first
    part of the name for each of them.  This is especially helpful when pooling the results
    of runs on different alignment files.  See the OUTPUT section for more details.

-a, --annotation_file     required
    Path to a coordinate-sorted GFF3 annotation file.  Follow the spec.

-k, --annotation_key      optional
    There is no standard place in the 9th column of GFF3 where an annotation string like
    a gene product name should be kept.  The most consistent is probably 'Note'.  This
    is a key=value pair column, so if you specify the key this script will find the values
    for each gene and include those in the analysis output.

======
OUTPUT
======

The script creates 5 files, named using the prefix specified by the -o option:

    $prefix.list.genes_missing
        A list of the genes completely missing from the mapping process.  This will not
        include those partically covered.  The columns of the file are:

            reference_id | gene_id | gene_fmin | gene_fmax | gene product
    
    $prefix.stats.gene_coverage
        Summary statistics for gene coverage in stat\tvalue format.  Example:

            Total molecule count            9
            Total gene count                4172
            Genes not covered at all        35
            Genes covered (partially)       169
            Genes covered (completely)      3968
    
    $prefix.stats.refmol_coverage
        Summary statistics for reference molecule coverage in stat\tvalue format.  Example:

            Total bases in reference molecules	        8347606
            Ref bases covered by query fragments	8181412
            Ref % covered by query fragments	        98.01
            Ref % identity by query fragments    	98.49

    $prefix.tab.refmol_coverage
        Here each reference molecule is listed along with its coverage statistics. The
        columns are:

            reference_id | reference_length | reference_bases_covered | percent_coverage
    
    $prefix.tab.extensions
        If the fragment alignments show evidence of extension on the ends of the reference
        molecule, these data are listed here.  The columns of the file are:

            reference_id | ref_fmin | ref_fmax | ref_strand | qry_id | qry_fmin | qry_fmax | qry_strand | qry_length

"""


## If a query molecule aligns to the reference at the end, this defines
#   the cutoff for the percentage of that fragment (that isn't hanging
#   off the end) that must align with the reference.  This should be
#   set pretty high.
fragment_overlap_perc_cutoff = 90

def parse_annotation(gff3_file, product_key):
    annotation = {}

    for line in open(gff3_file, 'r'):
        line.rstrip()
        cols = line.split("\t")

        if len(cols) != 9 or cols[2] != 'gene':
            continue

        asm_id = cols[0]

        if asm_id not in annotation:
            annotation[asm_id] = []

        gene_id = get_last_column_value( 'ID', cols[8] )

        if gene_id is None:
            print("ERROR: failed to get a gene ID for the following line: " + line)

        gene = {}
        gene['id'] = gene_id
        gene['fmin'] = int(cols[3]) - 1
        gene['fmax'] = int(cols[4])
        gene['cov'] = None

        ## grab the annotation, if a key was specified
        if product_key:
            gene['prod'] = get_last_column_value( product_key, cols[8] )
        else:
            gene['prod'] = None
        
        annotation[asm_id].append( gene )
    
    ## if debugging, report how many genes were found on each assembly
    #for asm_id in annotation:
        #print("DEBUG: from molecule ({0}), parsed ({1}) genes".format(asm_id, len(annotation[asm_id]) ) )

    return annotation


def get_last_column_value( key, colstring ):
    ## first check if there's a version that ends with the ';' symbol
    re_str = "{0}=(.+?)\;".format(key)
    m = re.match( re_str, colstring )
    key_value = None

    if m:
        key_value = m.group(1)
    else:
        re_str = "{0}=(.+?)".format(key)
        m = re.match( re_str, colstring )

        if m:
            key_value = m.group(1)

    return key_value


def calculate_gene_coverage( annot, frags ):
    """
    Iterates through the fragments aligned to a reference molecule and tracks the
    overlaps of each fragment with the genes that are annotated on that reference
    """
    ## this can be reduced to a sliding analysis window if this performs unreasonably
    for frag in frags:
        for gene in annot:
            #if gene['id'] == 'TP05_0045' and frag['id'] == 'ctg7180000023240':
            #    print("DEBUG: comparing the target gene and fragment")
            #    print("DEBUG: comparing r({0}-{1}) with q({2}-{3})".format(gene['fmin'],gene['fmax'],frag['rfmin'],frag['rfmax']) )

            ## don't check this one if it's already reported as completely covered by another fragment
            if gene['cov'] == 'complete':
                continue

            ## if the gene fmin falls within range of this fragment, it has at least partial coverage
            if gene['fmin'] >= frag['rfmin'] and gene['fmin'] <= frag['rfmax']:
                gene['cov'] = 'partial'

                ## if the gene fmax also falls within range, the gene is fully covered
                if gene['fmax'] <= frag['rfmax']:
                    gene['cov'] = 'complete'

            ## also check for fmax-only coverage of the gene
            elif gene['fmax'] >= frag['rfmin'] and gene['fmax'] <= frag['rfmax']:
                gene['cov'] = 'partial'


def report_gene_coverage_results( annot, stats_ofh, missing_ofh ):
    total_molecule_count = 0
    total_gene_count = 0
    partial_covered_gene_count = 0
    completely_covered_gene_count = 0
    completely_missed_gene_count = 0

    for assembly_id in annot:
        total_molecule_count += 1

        for gene in annot[assembly_id]:
            total_gene_count += 1

            if gene['cov'] == 'partial':
                partial_covered_gene_count += 1
            elif gene['cov'] == 'complete':
                completely_covered_gene_count += 1
            else:
                completely_missed_gene_count += 1
                missing_ofh.write("assembly:{0}\t{1}\t{2}\t{3}\t{4}\n".format( \
                        assembly_id, gene['id'], gene['fmin'], gene['fmax'], gene['prod']) )

    stats_ofh.write("Total molecule count\t{0}\n".format(total_molecule_count) )
    stats_ofh.write("Total gene count\t{0}\n".format(total_gene_count) )
    stats_ofh.write("Genes not covered at all\t{0}\n".format(completely_missed_gene_count) )
    stats_ofh.write("Genes covered (partially)\t{0}\n".format(partial_covered_gene_count) )
    stats_ofh.write("Genes covered (completely)\t{0}\n".format(completely_covered_gene_count) )
                    

def calculate_fragment_coverage( ref_id, frags, ref_length, cov_stats, covinfo_ofh, ext_ofh ):
    last_ref_coord = 0
    last_frag_rfmax = 0
    bases_covered = 0
    bases_uncovered = 0
    
    ## these are the possible number of bases that might be extended off the end of
    #   the reference sequence based on the fragments aligned.
    base_extension_fmin = 0
    base_extension_fmax = 0
    ## for calculating summary percent identity, taking into account overlaps
    last_frag_pctid = 0
    pctid_base_score = 0

    ## iterate through the fragments sorted by fmin coordinate on the reference
    for frag in sorted(frags, key=itemgetter('rfmin')):

        ## does the aligned query fragment overlap either end of the reference?  (suggesting extension?)
        possible_rfmin_ext = 0
        possible_rfmax_ext = 0
        if frag['qstrand'] == 1:
            possible_rfmin_ext = frag['qfmin'] - frag['rfmin']
            possible_rfmax_ext = (frag['rlen'] - frag['rfmax']) + (frag['qlen'] - frag['qfmax'])
        else:
            possible_rfmin_ext = frag['qlen'] - frag['qfmax'] - frag['rfmin']
            possible_rfmax_ext = frag['qfmax'] - (frag['rlen'] - frag['rfmin'])

        if possible_rfmin_ext > 0:
            #print("DEBUG: apparent 5' overhang between {0} and {1} at frag position {2}-{3}".format(ref_id, frag['id'], frag['qfmin'], frag['qfmax']))
            ## does it pass our alignment cutoff(s)? fragment_overlap_perc_cutoff
            bp_aligned = frag['qfmax'] - frag['qfmin']
            bp_unaligned_and_not_overhang = frag['qlen'] - bp_aligned - possible_rfmin_ext
            perc_unaligned_and_not_overhang = (bp_unaligned_and_not_overhang / frag['qlen']) * 100;
            
            if perc_unaligned_and_not_overhang <= (100 - fragment_overlap_perc_cutoff):
                ext_ofh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
                        ref_id, frag['rfmin'], frag['rfmax'], 1,\
                        frag['id'], frag['qfmin'], frag['qfmax'], frag['qstrand'], frag['qlen']))
                print("DEBUG: apparent 5' overhang ({4}bp) between {0} and {1} at frag position {2}-{3}".format(ref_id, frag['id'], frag['qfmin'], frag['qfmax'], possible_rfmin_ext))
        
        if possible_rfmax_ext > 0:
            bp_aligned = frag['qfmax'] - frag['qfmin']
            bp_unaligned_and_not_overhang = frag['qfmin'] + (frag['rlen'] - frag['rfmax'])
            perc_unaligned_and_not_overhang = (bp_unaligned_and_not_overhang / frag['qlen']) * 100;
        
            if perc_unaligned_and_not_overhang <= (100 - fragment_overlap_perc_cutoff):
                ext_ofh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( \
                        ref_id, frag['rfmin'], frag['rfmax'], 1,\
                        frag['id'], frag['qfmin'], frag['qfmax'], frag['qstrand'], frag['qlen']))
                print("DEBUG: apparent 3' overhang ({4}bp) between {0} and {1} at frag position {2}-{3}".format(ref_id, frag['id'], frag['qfmin'], frag['qfmax'], possible_rfmax_ext))
        


        ## for fragments contained within a previous fragment (ex: 0-10, then 3-5)
        #   i'm actually not sure if nucmer lets this happen, but it's caught anyway
        if frag['rfmax'] <= last_ref_coord:
            pass
        
        ## for fragments starting after or even with the last fragment end (ex: 1-3, then 5-10 or 3-10)
        elif frag['rfmin'] >= last_ref_coord:
            bases_covered_by_this_frag = frag['rfmax'] - frag['rfmin']
            bases_covered += bases_covered_by_this_frag
            bases_uncovered += (frag['rfmin'] - last_ref_coord)
            last_ref_coord = frag['rfmax']
            pctid_base_score += (bases_covered_by_this_frag * frag['pctid'])
            last_frag_rfmax = frag['rfmax']
            last_frag_pctid = frag['pctid']
            
        ## for a new fragment overlapping the previous (ex: 1-5, then 3-8)
        elif frag['rfmin'] < last_ref_coord:
            bases_covered_by_this_frag = frag['rfmax'] - last_ref_coord
            bases_covered += bases_covered_by_this_frag
            pctid_base_score += (bases_covered_by_this_frag * frag['pctid'])

            ## if this fragment has a higher percent identity than the previous one
            #  pct_id score needs to be corrected/increased
            if frag['pctid'] > last_frag_pctid:
                adjustment = (last_frag_rfmax - frag['rfmin']) * (frag['pctid'] - last_frag_pctid)
                print("adjustment:{0} because {4} < {5}: ({3} - {4}) * ({1} - {2})".format( \
                       adjustment, frag['pctid'], last_frag_pctid, last_frag_rfmax, frag['rfmin'], \
                       last_ref_coord))
                pctid_base_score += adjustment
            
            last_ref_coord = frag['rfmax']
            last_frag_rfmax = frag['rfmax']
            last_frag_pctid = frag['pctid']
        else:
            raise Exception("ERROR: unhandled fragment overlap: [{0}-{1} and {2}-{3}]".format(frag['rfmin'], frag['rfmax'], frag['qfmin'], frag['qfmax']) )

    ## then, finally, see if there's a region from the end of the last fragment to the end
    #   of the reference
    if last_frag_rfmax < ref_length:
        bases_uncovered += ref_length - last_frag_rfmax
        
    perc_cov = (bases_covered / ref_length) * 100
    ## output: reference mol, ref mol size, # bases covered, percent coverage
    covinfo_ofh.write("{0}\t{1}\t{2}\t{3}\n".format(ref_id, ref_length, bases_covered, perc_cov))
    #print("DEBUG: {0} size: {1} has {2}bp covered, ({3:.2f}%)".format(ref_id, ref_length, bases_covered, perc_cov))

    cov_stats['n_total'] += ref_length
    cov_stats['n_cov'] += bases_covered
    cov_stats['n_uncov'] += bases_uncovered
    cov_stats['sum_pct_id'] = pctid_base_score / ref_length
    

def main():
    parser = argparse.ArgumentParser( description='Parses nucmer coords output to provide an overall coverage report')

    ## coords file generated with: show-coords -l -r -T out.delta
    parser.add_argument('-c', '--coords_file', type=str, required=True, \
                            help='Path to a nucmer coords file with non-overlapping results (requires -l -r -T options of show-coords)' )
    parser.add_argument('-o', '--output_prefix', type=str, required=True, help='Several output files will be created with this prefix.' )
    parser.add_argument('-a', '--annotation_file', type=str, required=True, help='Path to a sorted GFF3 annotation file' )
    parser.add_argument('-k', '--annotation_key', type=str, required=False, help='Optional.  Key string to look for in the 9th column of the GFF3 file for an annotation string.' )
    args = parser.parse_args()

    ## like: h[$assem] = [ {id=?,fmin=?,fmax=?}, ...  ]
    annot = parse_annotation( args.annotation_file, args.annotation_key )

    ## like: [ {id=?,qfmin=?,qfmax=?,rfmin=?,rfmax=?} ]
    query_fragments = []

    ## open the output files
    genecov_stats_ofh     = open(args.output_prefix + ".stats.gene_coverage", "wt")
    genesmissing_list_ofh = open(args.output_prefix + ".list.genes_missing", "wt")
    refmol_stats_ofh      = open(args.output_prefix + ".stats.refmol_coverage", "wt")
    refcov_stats_ofh      = open(args.output_prefix + ".tab.refmol_coverage", "wt")
    refext_list_ofh       = open(args.output_prefix + ".tab.extensions", "wt")
    refext_list_ofh.write("# {0}\n".format(args.output_prefix) )
    refext_list_ofh.write("# reference_id\tref_fmin\tref_fmax\tref_strand\tqry_id\tqry_fmin\tqry_fmax\tqry_strand\tqry_length\n");
    
    ref_cov_stats = { 'n_total': 0, 'n_cov': 0, 'n_uncov': 0 }

    alignment_lines_found = 0
    current_ref_id = None
    current_ref_length = None
    
    for line in open(args.coords_file, 'r'):
        cols = line.split()

        if len(cols) == 11:
            alignment_lines_found += 1
        else:
            continue

        cols[0] = int(cols[0])
        cols[1] = int(cols[1])
        cols[2] = int(cols[2])
        cols[3] = int(cols[3])
        
        if cols[9] != current_ref_id:
            if current_ref_id is not None:
                calculate_gene_coverage( annot[current_ref_id], query_fragments )
                calculate_fragment_coverage( current_ref_id, query_fragments, current_ref_length, ref_cov_stats, refcov_stats_ofh, refext_list_ofh )
                
            ## reste
            current_ref_id = cols[9]
            current_ref_length = int(cols[7])
            query_fragments = []
        
        qstrand = 1

        if cols[2] > cols[3]:
            qstrand = -1

        fragment = {}
        fragment['id'] = cols[10]
        fragment['qfmin']   = min(cols[2], cols[3]) - 1
        fragment['qfmax']   = max(cols[2], cols[3])
        fragment['qlen']    = int(cols[8])
        fragment['qstrand'] = qstrand
        fragment['rfmin']   = min(cols[0], cols[1]) - 1
        fragment['rfmax']   = max(cols[0], cols[1])
        fragment['rlen']    = int(cols[7])
        fragment['pctid']   = float(cols[6])
        query_fragments.append(fragment)

    ## don't forget the last one
    calculate_gene_coverage( annot[current_ref_id], query_fragments )
    calculate_fragment_coverage( current_ref_id, query_fragments, current_ref_length, ref_cov_stats, refcov_stats_ofh, refext_list_ofh )
    
    if alignment_lines_found == 0:
        raise Exception("ERROR: failed to find any 11-column alignment lines")
    else:
        print("INFO: {0} alignment lines found".format(alignment_lines_found) )

    report_gene_coverage_results( annot, genecov_stats_ofh, genesmissing_list_ofh )
    
    cov_perc = (ref_cov_stats['n_cov'] / ref_cov_stats['n_total']) * 100
    refmol_stats_ofh.write("Total bases in reference molecules\t{0}\n".format(ref_cov_stats['n_total']) )
    refmol_stats_ofh.write("Ref bases covered by query fragments\t{0}\n".format(ref_cov_stats['n_cov']) )
    refmol_stats_ofh.write("Ref % covered by query fragments\t{0:.2f}\n".format(cov_perc))
    refmol_stats_ofh.write("Ref % identity by query fragments\t{0:.2f}\n".format(ref_cov_stats['sum_pct_id']))


if __name__ == '__main__':
    main()







