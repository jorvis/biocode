#!/usr/bin/env python3.2

import argparse
import os
import re
import sqlite3
from collections import OrderedDict

from biocode.utils import read_list_file


def main():
    """This is the second script I've written in Python.  I'm sure it shows."""
    parser = argparse.ArgumentParser( description='Reads a BLAST m8 file and taxonomy DB to produce a taxonomic profile at any user-specified ranking level.')

    ## input formats: btab, blast_m8
    parser.add_argument('-f', '--input_format', type=str, required=True, help='Blast format: current options are btab or blast_m8' )

    ## The SQLite3 file that will be read for taxonomy information
    parser.add_argument('-t', '--taxonomy_db', type=str, required=True, help='Path to a taxonomy.db file created by "create_taxonomy_db.py"' )

    ## BLAST list file
    parser.add_argument('-b', '--blast_list_file', type=str, required=True, help='List of BLAST files (m8 format)' )

    ## output file to be written
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path where the result file should written' )

    ## E-value cutoff to use
    parser.add_argument('-e', '--eval_cutoff', type=float, required=False, help='Optional E-value cutoff to use.' )

    ## Top N hits per query to score.  Only counts those where the taxon could be looked up in the indexes
    parser.add_argument('-n', '--top_n', type=int, required=False, default=1, help=' Top N hits per query to score.  Only counts unique taxon matches which could be looked up in the indexes' )

    ## rank on which matches will be grouped and reported.  values like: species, genus, order, family, etc.
    parser.add_argument('-r', '--rank', type=str, required=True, help='Taxonomy rank on which to group all matches, such as: species, genus, order, family, etc.' )
    args = parser.parse_args()

    conn = sqlite3.connect( args.taxonomy_db )
    c = conn.cursor()

    blast_files = read_list_file( args.blast_list_file )

    taxon_counts = {}
    processed_file_count = 0
    stats = {}
    stats['gi_lookup_success_count'] = 0
    stats['gi_lookup_fail_count'] = 0
    stats['taxon_lookup_success_count'] = 0
    stats['taxon_lookup_failure_count'] = 0
     
    for file in blast_files:
        print("Processing file: ", file)
        if args.input_format == 'blast_m8' or args.input_format == 'btab':
            parse_blast_file( file, c, taxon_counts, args.eval_cutoff, args.input_format, stats, args.top_n )
        else:
            raise Exception("Unsupported input format passed: {0}".format(args.input_format) )

        processed_file_count += 1
        #if processed_file_count == 50:
            #break
   
    ## process the taxon counts, conforming them to the user-specified rank
    result_table = group_taxa_by_rank( args.rank, taxon_counts, c )
    node_names = get_selected_node_names( result_table, c )

    c.close()

    fout = open(args.output_file, mode='w')

    ## write the results to the output file in order of most-found clade first
    for tax_id in OrderedDict(sorted(result_table.items(), reverse=True, key=lambda t: t[1])):
        sci_name = ''
        if tax_id in node_names:
            sci_name = node_names[tax_id]
            
        fout.write( "{0}\t{1}\t{2}\n".format(tax_id, int(result_table[tax_id]), sci_name ) )

    fout.close()

    print("INFO: successful GI lookups: {0}/{1}".format(stats['gi_lookup_success_count'], \
                                                        (stats['gi_lookup_fail_count'] + stats['gi_lookup_success_count'])) )
    print("INFO: successful taxon lookups: {0}/{1}".format( stats['taxon_lookup_success_count'], \
                                                        (stats['taxon_lookup_success_count'] + stats['taxon_lookup_failure_count']) ) )
    

def get_selected_node_names( res_table, cursor):
    node_names = {}

    for taxon_id in res_table:
        cursor.execute('''SELECT scientific_name FROM orgs WHERE tax_id=?''', (taxon_id,) )
        row = cursor.fetchone()
        
        if row:
            node_names[taxon_id] = row[0]
        else:
            print("WARN: failed to get scientific name for tax_id:", taxon_id)

    return node_names


def get_ranked_taxon( tax_id, c, rank, rec_depth ):
    c.execute("""SELECT parent_tax_id, rank FROM nodes WHERE tax_id = ?""", (tax_id,) )
    row = c.fetchone()
    
    if rec_depth > 20:
        print("WARN: deep recursion detected for tax ID:", tax_id)
        return None

    if row:
        if row[1] == rank:
            return tax_id
        else:
            return get_ranked_taxon( row[0], c, rank, rec_depth + 1 )
    else:
        print("WARN: unable to find ranked taxon for tax_id:", tax_id)
        return None


def group_taxa_by_rank(rank, counts, cursor):
    """Given a taxonomic rank, the input count table is regrouped by walking
       up the taxonomy tree until all nodes are at the level of the passed
       rank
    """
    ranked_counts = dict()
    unranked_taxon_count = 0

    for taxon_id in counts:
        ranked_taxon_id = get_ranked_taxon(taxon_id, cursor, rank, 0)
        if ranked_taxon_id:
            if ranked_taxon_id in ranked_counts:
                #print("DEBUG: increasing count for taxon: {0} by: {1}".format(taxon_id, counts[taxon_id]['n']) )
                ranked_counts[ranked_taxon_id] += counts[taxon_id]['n']
            else:
                #print("DEBUG: initializing a count for ranked_taxon_id {0} from taxon_id {1}".format(ranked_taxon_id, taxon_id) )
                ranked_counts[ranked_taxon_id] = counts[taxon_id]['n']
        else:
            unranked_taxon_count += 1

    return ranked_counts


def parse_blast_file( file, cursor, tax, eval_cutoff, format, stats, hits_per_query ):
    """ For each query sequence find the top match above the E-val cutoff (if any)
        which has an NCBI taxonomy assignment.
    """
    if ( not os.path.isfile(file) ):
        raise Exception("Couldn't find file: " + file)

    ## presets are for ncbi_m8, for which lines should have 12 columns.  they are
    # 1: Query - The query sequence id
    # 2: Subject - The matching subject sequence id
    # 3: Percent identity
    # 4: alignment length
    # 5: mismatches
    # 6: gap openings
    # 7: q.start
    # 8: q.end
    # 9: s.start
    # 10: s.end
    # 11: e-value
    # 12: bit score
    ID_COLUMN_NUM = 0
    SUBJECT_LABEL_COLUMN_NUM = 1
    EVAL_COLUMN_NUM = 10
    ALIGN_LEN_COLUMN_NUM = 3
    BIT_SCORE_COLUMN_NUM = 11

    if format == 'btab':
        ID_COLUMN_NUM = 0
        SUBJECT_LABEL_COLUMN_NUM = 5
        EVAL_COLUMN_NUM = 19

    current_id = ""
    current_id_classified = False
    current_id_match_count = 0
    current_match_ids = dict()

    for line in open(file, "r"):
        cols = line.split("\t")
        if len(cols) >= 10:
            this_id = cols[ID_COLUMN_NUM]

            ## this controls that we only look at the top hit for query
            if this_id != current_id:
                current_id_classified = False
                current_id_match_ids = dict()
                current_id_match_count = 0
                current_id = this_id
            
            if current_id_match_count < hits_per_query:
                if eval_cutoff is None or eval_cutoff >= float(cols[EVAL_COLUMN_NUM]):
                    #print("DEBUG: attempting to parse a GI for header: ({0})".format(cols[SUBJECT_LABEL_COLUMN_NUM]) )
                    gi = parse_gi(cols[SUBJECT_LABEL_COLUMN_NUM], cursor)
                                        
                    if gi:
                        #print("DEBUG: Got a GI ({0}) for hit with id this_id".format(gi))
                        stats['gi_lookup_success_count'] += 1
                        taxon_id = get_taxon_id_by_gi(gi, cursor)

                        if taxon_id:
                            stats['taxon_lookup_success_count'] += 1
                            if taxon_id not in current_match_ids:
                                #print("DEBUG: adding match to taxon_id: {0}".format(taxon_id) )
                                match_score = int(cols[BIT_SCORE_COLUMN_NUM])/int(cols[ALIGN_LEN_COLUMN_NUM])
                                add_taxon_match( taxon_id, tax, cursor, match_score )
                                current_id_match_count += 1
                                current_match_ids[taxon_id] = True
                        else:
                            stats['taxon_lookup_failure_count'] += 1
                            print("WARN: failed to find a taxon_id for gi: {0}".format(gi))
                    else:
                        stats['gi_lookup_fail_count'] += 1


def add_taxon_match( id, tax, c, score ):
    if id in tax:
        tax[id]['n'] += score
        #print("DEBUG: match count for taxon id {0} increased to {1}".format(id, tax[id]['n']) )
    else:
        tax[id] = {}
        tax[id]['n'] = score
        tax[id]['l'] = get_clade_name_by_taxon_id( c, id )


def get_clade_name_by_taxon_id( cursor, taxon_id ):
    cursor.execute("""SELECT scientific_name FROM orgs WHERE tax_id = ?""", (taxon_id,) )
    row = cursor.fetchone()

    if row:
        return row[0]
    else:
        return None


def get_taxon_id_by_gi( gi, cursor ):
    """Attempts to fetch a taxon ID using the GI from first the nucl_acc and then prot_acc tables"""
    taxon_id = None
    
    cursor.execute( 'SELECT tax_id FROM nucl_acc WHERE gi = ?', (gi,) )
    row = cursor.fetchone()

    if row:
        taxon_id = row[0]
    else:
        cursor.execute( 'SELECT tax_id FROM prot_acc WHERE gi = ?', (gi,) )
        row = cursor.fetchone()

        if row:
            taxon_id = row[0]

    return taxon_id



def get_gi_by_accession( acc, cursor ):
    ## The table we query depends on the source
    #   Protein queries have a three-letter prefix, followed by five digits.
    m = re.match( '^[A-Z]{3}[0-9]{5}\.\d', acc )

    ## CURRENT HACK, because of a db index bug the version number was indexed with
    #   the period character removed.
    acc = acc.replace('.', '')

    if m:
        cursor.execute("""SELECT gi FROM prot_acc WHERE version = ?""", (acc,) )
        row = cursor.fetchone()
    else:
        cursor.execute("""SELECT gi FROM nucl_acc WHERE version = ?""", (acc,) )
        row = cursor.fetchone()

    if row:
        return row[0]
    else:
        return None



def parse_gi( match_header, cursor ):
    """Parses the GI number out of a NCBI-formatted BLAST match header string.  If the
       header has an accession instead of a GI, it performs a lookup
    """
    ## look for a GI first
    #print("DEBUG: looking for a gi in header: ({0})".format(match_header))
    m = re.match( 'gi\|(\d+)', match_header  )
    gi = None
    accession = None

    if m:
        gi = m.group(1)
    else:
        ## a direct GI pull failed, so look for something that resembles an accession instead
        m = re.match( '[a-z]+\|([A-Z]{2,}\d{5,}\.\d)', match_header )

        if m:
            accession = m.group(1)
            gi = get_gi_by_accession( accession, cursor )

    #print("DEBUG:\treturning GI ({0}) for header ({1}), accession ({2})".format( gi, match_header, accession ) )
    return gi

    


if __name__ == '__main__':
    main()







