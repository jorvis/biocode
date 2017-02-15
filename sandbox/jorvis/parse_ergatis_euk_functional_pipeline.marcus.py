#!/usr/bin/env python3

import argparse
import os
import re
import sqlite3
import sys
from collections import defaultdict

from biocode import utils, annotation, gff, things

## constants
DEFAULT_PRODUCT_NAME = "hypothetical protein"

next_ids = defaultdict(int)

def main():
    parser = argparse.ArgumentParser( description='Parses multiple sources of evidence to generate a consensus functional annotation')

    ## output file to be written
    parser.add_argument('-f', '--input_fasta', type=str, required=True, help='Protein FASTA file of source molecules' )
    parser.add_argument('-m', '--hmm_htab_list', type=str, required=True, help='List of htab files from hmmpfam3' )
    parser.add_argument('-bs', '--blast_sprot_btab_list', type=str, required=True, help='List of btab files from BLAST against UniProtKB/SWISS-PROT' )
    parser.add_argument('-bt', '--blast_trembl_btab_list', type=str, required=False, help='List of btab files from BLAST against UniProtKB/Trembl' )
    parser.add_argument('-d', '--hmm_db', type=str, required=True, help='SQLite3 db with HMM information' )
    parser.add_argument('-u', '--uniprot_sprot_db', type=str, required=True, help='SQLite3 db with UNIPROT/SWISSPROT information' )
    parser.add_argument('-a', '--format', type=str, required=False, default='tab', help='Output format.  Current options are: "tab", "fasta", "gff3"' )
    parser.add_argument('-s', '--source_gff', type=str, required=False, help='Source GFF file from which proteins were derived.  Required if you want to export any format other than tab-delimited.' )
    parser.add_argument('-e', '--blast_eval_cutoff', type=float, required=False, default=1e-5, help='Skip BLAST hits unless they have an E-value at least as low as this' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Optional output file path (else STDOUT)' )
    parser.add_argument('-r', '--organism_table', type=str, required=False, help='Optional table with counts of organism frequency based on top BLAST match for each protein' )
    parser.add_argument('-g', '--genomic_fasta', type=str, required=False, help='If passed, the genomic FASTA sequence will be included in the exported GFF3')
    args = parser.parse_args()

    check_arguments(args)

    # connection to the HMM-associated SQLite3 database
    hmm_db_conn = sqlite3.connect(args.hmm_db)
    hmm_db_curs = hmm_db_conn.cursor()

    # connection to the UniProt_Sprot SQLite3 database
    usp_db_conn = sqlite3.connect(args.uniprot_sprot_db)
    usp_db_curs = usp_db_conn.cursor()

    sources_log_fh = open("{0}.sources.log".format(args.output_file), 'wt')
    
    # this is a dict of biothings.Polypeptide objects
    polypeptides = initialize_polypeptides( sources_log_fh, args.input_fasta )

    # keyed on polypeptide ID, the values here are the organism name for the top BLAST match of each
    polypeptide_blast_org = dict()

    # get source structural annotation, if necessary:
    if args.source_gff is not None:
        print("INFO: parsing source GFF")
        (assemblies, features) = gff.get_gff3_features(args.source_gff)

    print("INFO: parsing HMM evidence")
    parse_hmm_evidence( sources_log_fh, polypeptides, args.hmm_htab_list, hmm_db_curs )
    hmm_db_curs.close()

    print("INFO: parsing BLAST (SWISS-PROT) evidence")
    parse_sprot_blast_evidence( sources_log_fh, polypeptides, polypeptide_blast_org, args.blast_sprot_btab_list, usp_db_curs, args.blast_eval_cutoff )
    usp_db_curs.close()

    if args.blast_trembl_btab_list is not None:
        print("INFO: parsing BLAST (TrEMBL) evidence")
        parse_trembl_blast_evidence(polypeptides, args.blast_trembl_btab_list, args.blast_eval_cutoff)

    ## output will either be a file or STDOUT
    print("INFO: writing output")
    fout = sys.stdout
    if args.output_file is not None:
        fout = open(args.output_file, 'wt')

    if args.format == 'tab':
        write_tab_results(fout, polypeptides)
    elif args.format == 'fasta':
        write_fasta_results(fout, polypeptides)
    elif args.format == 'gff3':
        write_gff3_results(fout, polypeptides, assemblies, features, args.genomic_fasta)
    
    fout.close()

    if args.organism_table is not None:
        create_organism_table(args.organism_table, polypeptide_blast_org)

def check_arguments( args ):
    # Check the acceptable values for format
    if args.format not in ('tab', 'gff3', 'fasta'):
        raise Exception("ERROR: The format provided ({0}) isn't supported.  Please check the documentation again.".format(args.format))

    # The --source_gff argument is required if --format=gff3
    if args.source_gff is None and args.format == 'gff3':
        raise Exception("ERROR: You must pass the --source_gff argument if --format=gff3")

    # It doesn't make sense to pass --genomic_fasta if the --format isn't GFF3
    if args.genomic_fasta is not None and args.format != 'gff3':
        raise Exception("ERROR: Argument --genomic_fasta not valid unless --format=gff3")

    # make sure all passed files are found before we do anything else:
    for path in ( args.input_fasta, args.hmm_htab_list, args.blast_sprot_btab_list, args.hmm_db, args.uniprot_sprot_db ):
        if not os.path.isfile( path ):
            raise Exception("ERROR: You passed this file but the script failed to find it: {0}".format(path))


def write_gff3_results( f, polypeptides, assemblies, features, genomic_fasta ):
    f.write("##gff-version 3\n")
    
    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            for mRNA in gene.mRNAs():
                # the polypeptide ID needs to be the same as the mRNA one but with the type changed
                #polypeptide_id = mRNA.id.replace('mRNA', 'polypeptide')
                polypeptide_id = "{0}-Protein".format(mRNA.id)
                
                # add polypeptides to the mRNAs
                if polypeptide_id in polypeptides:
                    # the polypeptide ID needs to be the same as the mRNA one but with the type changed
                    new_polypeptide_id = mRNA.id.replace('mRNA', 'polypeptide').rstrip('-Protein')
                    polypeptides[polypeptide_id].id = new_polypeptide_id
                    mRNA.add_polypeptide( polypeptides[polypeptide_id] )
            
            gene.print_as(fh=f, source='IGS', format='gff3')

    if genomic_fasta is not None:
        f.write("##FASTA\n")

        for line in open(genomic_fasta):
            f.write(line)


def write_fasta_results( f, polypeptides ):
    """
    Produces headers like:
    >ID PRODUCT_NAME gene::GENE_SYMBOL ec::EC_NUMBERS go::GO_TERMS

    Example:
    
    """
    for polypeptide_id in polypeptides:
        polypeptide = polypeptides[polypeptide_id]
        go_string = ""
        ec_string = ""

        for go_annot in polypeptide.annotation.go_annotations:
            go_string += "GO:{0},".format(go_annot.go_id)
        
        go_string = go_string.rstrip(',')

        for ec_annot in polypeptide.annotation.ec_numbers:
            ec_string += "{0},".format(ec_annot.number)
        
        ec_string = ec_string.rstrip(',')

        header = "{0} {1}".format(polypeptide_id, polypeptide.annotation.product_name)

        if polypeptide.annotation.gene_symbol is not None:
            header = "{0} gene::{1}".format(header, polypeptide.annotation.gene_symbol)

        if ec_string != "":
            header = "{0} ec::{1}".format(header, ec_string)
            
        if go_string != "":
            header = "{0} go::{1}".format(header, go_string)
            
        f.write( ">{0}\n".format( header ) )
        f.write( "{0}\n".format(utils.wrapped_fasta(polypeptide.residues)))
        
        


def write_tab_results( f, polypeptides ):
    f.write("# polypeptide ID\tpolypeptide_length\tgene_product_name\tGO_terms\tEC_nums\tgene_symbol\n")

    for polypeptide_id in polypeptides:
        polypeptide = polypeptides[polypeptide_id]
        go_string = ""
        ec_string = ""

        for go_annot in polypeptide.annotation.go_annotations:
            go_string += "GO:{0},".format(go_annot.go_id)
        
        go_string = go_string.rstrip(',')

        for ec_annot in polypeptide.annotation.ec_numbers:
            ec_string += "{0},".format(ec_annot.number)
        
        ec_string = ec_string.rstrip(',')
                
        f.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format( \
                polypeptide_id, polypeptide.length, polypeptide.annotation.product_name, go_string, ec_string, \
                polypeptide.annotation.gene_symbol) )


def initialize_polypeptides( log_fh, fasta_file ):
    '''
    Reads a FASTA file of (presumably) polypeptide sequences and creates a dict of Polypeptide
    objects, keyed by ID.  No bioannotation.FunctionalAnnotation objects will be attached yet.
    '''
    seqs = utils.fasta_dict_from_file(fasta_file)

    polypeptides = dict()

    for seq_id in seqs:
        polypeptide = things.Polypeptide(id=seq_id, length=len(seqs[seq_id]['s']), residues=seqs[seq_id]['s'])
        annotation = annotation.FunctionalAnnotation(product_name=DEFAULT_PRODUCT_NAME)
        log_fh.write("INFO: {0}: Set initial product name to '{1}'\n".format(seq_id, DEFAULT_PRODUCT_NAME))
        polypeptide.annotation = annotation
        
        polypeptides[seq_id] = polypeptide
    
    return polypeptides

def parse_trembl_blast_evidence(polypeptides, blast_list, eval_cutoff):
    '''
    Reads a list file of NCBI BLAST evidence against TrEMBL and a dict of polypeptides,
    populating each with Annotation evidence where appropriate.  Only attaches evidence if
    the product name is the default.

    Currently only considers the top BLAST hit for each query.
    '''
    for file in utils.read_list_file(blast_list):
        last_qry_id = None
        
        for line in open(file):
            line = line.rstrip()
            cols = line.split("\t")

            # We're going to ignore any lines which have 'uncharacterized' in the name
            if 'ncharacterized' in cols[15]:
                continue
            
            this_qry_id = cols[0]

            # skip this line if it doesn't meet the cutoff
            if float(cols[19]) > eval_cutoff:
                continue

            # the BLAST hits are sorted already with the top hit for each query first
            if last_qry_id != this_qry_id:
                annot = polypeptides[this_qry_id].annotation

                # get the accession from the cols[5]
                #  then process for known accession types
                accession = cols[5]

                # save it, unless the gene product name has already changed from the default
                if annot.product_name == DEFAULT_PRODUCT_NAME:
                    # current hack until DB is updated:
                    # some products look like this:
                    #    Coatomer subunit gamma-2 OS=Bos taurus GN=COPG2 PE=2 SV=1
                    # take off everything after the OS=
                    m = re.search("(.+) OS=", cols[15])

                    if m:
                        annot.product_name = m.group(1)
                    else:
                        annot.product_name = cols[15]

                # remember the ID we just saw
                last_qry_id = this_qry_id

def parse_sprot_blast_evidence( log_fh, polypeptides, blast_org, blast_list, cursor, eval_cutoff ):
    '''
    Reads a list file of NCBI BLAST evidence and a dict of polypeptides, populating
    each with Annotation evidence where appropriate.  Only attaches evidence if
    the product name is the default.

    Currently only considers the top BLAST hit for each query.
    '''
    for file in utils.read_list_file(blast_list):
        last_qry_id = None
        
        for line in open(file):
            line = line.rstrip()
            cols = line.split("\t")
            this_qry_id = cols[0]

            # skip this line if it doesn't meet the cutoff
            if float(cols[19]) > eval_cutoff:
                continue

            # the BLAST hits are sorted already with the top hit for each query first
            if last_qry_id != this_qry_id:
                annot = polypeptides[this_qry_id].annotation

                # get the accession from the cols[5]
                #  then process for known accession types
                accession = cols[5]

                if accession.startswith('sp|'):
                    # pluck the second part out of this:
                    #  sp|Q4PEV8|EIF3M_USTMA
                    accession = accession.split('|')[1]

                assertions = get_uspdb_annot( accession, cursor )
                blast_org[this_qry_id] = assertions['organism']
                    
                # save it, unless the gene product name has already changed from the default
                if annot.product_name == DEFAULT_PRODUCT_NAME:
                    # current hack until DB is updated:
                    # some products look like this:
                    #    Coatomer subunit gamma-2 OS=Bos taurus GN=COPG2 PE=2 SV=1
                    # take off everything after the OS=
                    m = re.search("(.+) OS=", cols[15])

                    if m:
                        annot.product_name = m.group(1)
                    else:
                        annot.product_name = cols[15]

                    log_fh.write("INFO: {0}: Updated product name to '{1}' based on BLAST hit to SPROT accession '{2}'".format(this_qry_id, annot.product_name, accession))

                # if no EC numbers have been set, they can inherit from this
                if len(annot.ec_numbers) == 0:
                    for ec_annot in get_uspdb_ec_nums( accession, cursor ):
                        annot.add_ec_number(ec_annot)

                # if no GO IDs have been set, they can inherit from this
                if len(annot.go_annotations) == 0:
                    for go_annot in get_uspdb_go_terms( accession, cursor ):
                        annot.add_go_annotation(go_annot)

                # if no gene symbol has been set, it can inherit from this
                if annot.gene_symbol is None:
                    annot.gene_symbol = assertions['symbol']

                # remember the ID we just saw
                last_qry_id = this_qry_id





def parse_hmm_evidence( log_fh, polypeptides, htab_list, cursor ):
    '''
    Reads a list file of HMM evidence and dict of polypeptides, populating each with
    Annotation evidence where appropriate.  Each file in the list can have results
    for multiple queries, but it's assumed that ALL candidate matches for any given
    query are grouped together.

    Currently only the top hit for any given query polypeptide is used.
    '''
    for file in utils.read_list_file(htab_list):
        last_qry_id = None
        
        for line in open(file):
            line = line.rstrip()
            cols = line.split("\t")
            
            ## only consider the row if the total score is above the total trusted cutoff
            if cols[12] >= cols[17]:
                continue

            this_qry_id = cols[5]
            accession = cols[0]
            version = None

            # if this is a PFAM accession, handle the version
            m = re.match("^(PF\d+)\.\d+", accession)
            if m:
                version = accession
                accession = m.group(1)

            ## the HMM hits are sorted already with the top hit for each query first
            if last_qry_id != this_qry_id:
                ## save it
                annot = polypeptides[this_qry_id].annotation
                annot.product_name = cols[15]
                log_fh.write("INFO: {0}: Updated product name to '{1}' based on HMM hit to accession '{2}'".format(this_qry_id, annot.product_name, accession))
                
                # does our hmm database provide GO terms for this accession?
                for go_annot in get_hmmdb_go_terms( accession, cursor ):
                    annot.add_go_annotation(go_annot)

                ## remember the ID we just saw
                last_qry_id = this_qry_id


def get_uspdb_annot( acc, c ):
    """
    This returns a dict of any annotation assertions for a given accession
    in which there will only be one.  Examples: organism name, gene symbol, etc.

    Other methods pull GO terms and EC numbers
    """
    qry = """
          SELECT us.organism, us.symbol
            FROM uniprot_sprot us
                 JOIN uniprot_sprot_acc us_acc ON us.id = us_acc.id
           WHERE us_acc.accession = ?
          """
    c.execute(qry, (acc,))
    #print("DEBUG: executing annot query where accession = ({0})".format(acc))

    assertions = { 'organism':None, 'symbol':None }

    for row in c:
        assertions['organism'] = row[0]
        assertions['symbol'] = row[1]
        break

    return assertions


def get_uspdb_ec_nums( acc, c ):
    """
    This returns a lsit of bioannotation:ECAnnotation objects
    """
    qry = """
          SELECT us_ec.ec_num
            FROM uniprot_sprot_ec us_ec
                 JOIN uniprot_sprot_acc us_acc ON us_ec.id=us_acc.id
           WHERE us_acc.accession = ?
          """
    c.execute(qry, (acc,))

    ec_annots = list()

    for row in c:
        ec = annotation.ECAnnotation(number=row[0])
        ec_annots.append(ec)

    return ec_annots


def get_uspdb_go_terms( acc, c ):
    """
    This returns a list of bioannotation:GOAnnotation objects
    """
    qry = """
          SELECT us_go.go_id
            FROM uniprot_sprot_go us_go
                 JOIN uniprot_sprot_acc us_acc ON us_go.id=us_acc.id
           WHERE us_acc.accession = ?
          """
    c.execute(qry, (acc,))

    go_annots = list()
    
    for row in c:
        go = annotation.GOAnnotation(go_id=row[0], ev_code='ISA', with_from=acc)
        go_annots.append(go)
    
    return go_annots

    


def get_hmmdb_go_terms( acc, c ):
    """
    This returns a list of bioannotation:GOAnnotation objects
    """
    qry = "SELECT hg.go_id FROM hmm_go hg JOIN hmm ON hg.hmm_id=hmm.id WHERE hmm.accession = ?"
    c.execute(qry, (acc,))

    go_annots = list()
    
    for row in c:
        go = annotation.GOAnnotation(go_id=row[0], ev_code='ISM', with_from=acc)
        go_annots.append(go)
    
    return go_annots


def create_organism_table( file_path, blast ):
    org = defaultdict(int)
    
    for protein_id in blast:
        org[ blast[protein_id] ] += 1

    table = open(file_path, 'wt')

    for name in sorted(org, key=org.get, reverse=True):
        table.write("{0}\t{1}\n".format(name, org[name]))
    
    table.close()


if __name__ == '__main__':
    main()







