#!/usr/bin/env python3

import argparse
import math
import os
import re
import sqlite3
import sys
from collections import defaultdict

from biocode import annotation, gff, things, utils

## constants
DEFAULT_PRODUCT_NAME = "hypothetical protein"

## example invocation:
# parse_ergatis_euk_functional_pipeline.py -f HUZ239124.metagenemark.nt.faa -m hmmpfam3.htab.list -bs ncbi-blastp.btab.list --hmm_db /usr/local/projects/jorvis/dbs/coding_hmm_lib.sqlite3 -u /usr/local/projects/jorvis/dbs/uniprot_sprot.sqlite3 -a gff3 -s HUZ239124.metagenemark.nt.gff3 -e 1e-10 -o HUZ239124.metagenemark.nt.annotated.gff3 -g HUZ239124.fas

# /home/jorvis/git/biocode/sandbox/jorvis/parse_ergatis_euk_functional_pipeline.py -f SRS019986.scaffolds.metagenemark.faa -s SRS019986.scaffolds.metagenemark.gff3 -e 1e-10 -g SRS019986.scaffolds.fa -m hmmpfam3.htab.list --hmm_db /usr/local/projects/jorvis/dbs/coding_hmm_lib.sqlite3 -o SRS019986.scaffolds.metagenemark.annotated.hmm.gff3  --format=gff3 

next_ids = defaultdict(int)

def main():
    parser = argparse.ArgumentParser( description='Parses multiple sources of evidence to generate a consensus functional annotation')

    ## output file to be written
    parser.add_argument('-f', '--input_fasta', type=str, required=True, help='Protein FASTA file of source molecules' )
    parser.add_argument('-m', '--hmm_htab_list', type=str, required=False, help='List of htab files from hmmpfam3' )
    parser.add_argument('-bs', '--blast_sprot_btab_list', type=str, required=False, help='List of btab files from BLAST against UniProtKB/SWISS-PROT' )
    parser.add_argument('-rs', '--rapsearch_sprot_btab_list', type=str, required=False, help='List of m8 files from RAPSEARCH2 against UniProtKB/SWISS-PROT' )
    parser.add_argument('-bt', '--blast_trembl_btab_list', type=str, required=False, help='List of btab files from BLAST against UniProtKB/Trembl' )
    parser.add_argument('-bk', '--blast_kegg_btab_list', type=str, required=False, help='List of btab files from BLAST against KEGG' )
    parser.add_argument('-bu100', '--blast_uniref100_btab_list', type=str, required=False, help='List of btab files from BLAST against UniRef100' )
    parser.add_argument('-ru100', '--rapsearch_uniref100_btab_list', type=str, required=False, help='List of m8 files from RAPSEARCH2 against UniRef100' )
    parser.add_argument('-u100f', '--uniref100_fasta', type=str, required=False, help='Only required if also passing RAPSEARCH2 against UniRef100 evidence' )
    parser.add_argument('-tm', '--tmhmm_raw_list', type=str, required=False, help='List of raw files from a tmhmm search' )
    parser.add_argument('-d', '--hmm_db', type=str, required=False, help='SQLite3 db with HMM information' )
    parser.add_argument('-u', '--uniprot_sprot_db', type=str, required=False, help='SQLite3 db with UNIPROT/SWISSPROT information' )
    parser.add_argument('-ur', '--uniref_db', type=str, required=False, help='SQLite3 db with UNIREF information' )
    parser.add_argument('-a', '--format', type=str, required=False, default='tab', help='Output format.  Current options are: "tab", "fasta", "gff3"' )
    parser.add_argument('-s', '--source_gff', type=str, required=False, help='Source GFF file from which proteins were derived.  Required if you want to export any format other than tab-delimited.' )
    parser.add_argument('-e', '--blast_eval_cutoff', type=float, required=False, default=1e-5, help='Skip BLAST hits unless they have an E-value at least as low as this' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Optional output file path (else STDOUT)' )
    parser.add_argument('-r', '--organism_table', type=str, required=False, help='Optional table with counts of organism frequency based on top BLAST match for each protein' )
    parser.add_argument('-g', '--genomic_fasta', type=str, required=False, help='If passed, the genomic FASTA sequence will be included in the exported GFF3')
    parser.add_argument('-eon', '--export_organism_names', help='If passed, includes organism names from top BLAST hit into 9th column when available.  Mostly useful for metagenomic samples.', action='store_true')
    args = parser.parse_args()

    check_arguments(args)

    # If --rapsearch_uniref100_btab_list passed, --uniref100_fasta is required
    if args.rapsearch_uniref100_btab_list is not None:
        if args.uniref100_fasta is None:
            raise Exception("ERROR: --uniref100_fasta required if --rapsearch_uniref100_btab_list is passed")

    sources_log_fh = open("{0}.sources.log".format(args.output_file), 'wt')
    
    # this is a dict of biothings.Polypeptide objects
    polypeptides = initialize_polypeptides( sources_log_fh, args.input_fasta )

    # Keyed on polypeptide ID (from the FASTA, which is actually the mRNA gff feature ID), the
    #  values here are the organism name for the top BLAST match of each
    polypeptide_blast_org = dict()

    # get source structural annotation, if necessary:
    if args.source_gff is not None:
        print("INFO: parsing source GFF")
        (assemblies, features) = gff.get_gff3_features(args.source_gff)

    if args.hmm_htab_list is not None:
        # connection to the HMM-associated SQLite3 database
        hmm_db_conn = sqlite3.connect(args.hmm_db)
        hmm_db_curs = hmm_db_conn.cursor()
        
        if args.hmm_db is None:
            raise Exception("ERROR: You specified HMM results but not the db with the -d option")
        
        print("INFO: parsing HMM evidence")
        parse_hmm_evidence( sources_log_fh, polypeptides, args.hmm_htab_list, hmm_db_curs )
        hmm_db_curs.close()

    if args.blast_sprot_btab_list is not None:
        if args.uniprot_sprot_db is None:
            raise Exception("ERROR: You specified BLAST evidence vs UnitProt/SwissProt results but not the db with the -u option")
        
        # connection to the UniProt_Sprot SQLite3 database
        usp_db_conn = sqlite3.connect(args.uniprot_sprot_db)
        usp_db_curs = usp_db_conn.cursor()
        print("INFO: parsing BLAST (SWISS-PROT) evidence")
        parse_sprot_blast_evidence( sources_log_fh, polypeptides, polypeptide_blast_org, args.blast_sprot_btab_list, usp_db_curs, args.blast_eval_cutoff, 'blast' )
        usp_db_curs.close()

    if args.rapsearch_sprot_btab_list is not None:
        if args.uniprot_sprot_db is None:
            raise Exception("ERROR: You specified RAPSEARCH2 evidence vs UnitProt/SwissProt results but not the db with the -u option")
        
        # connection to the UniProt_Sprot SQLite3 database
        usp_db_conn = sqlite3.connect(args.uniprot_sprot_db)
        usp_db_curs = usp_db_conn.cursor()
        print("INFO: parsing RAPSEARCH2 (SWISS-PROT) evidence")
        parse_sprot_blast_evidence( sources_log_fh, polypeptides, polypeptide_blast_org, args.rapsearch_sprot_btab_list, usp_db_curs, args.blast_eval_cutoff, 'rapsearch2' )
        usp_db_curs.close()

    if args.blast_trembl_btab_list is not None:
        print("INFO: parsing BLAST (TrEMBL) evidence")
        parse_trembl_blast_evidence(polypeptides, args.blast_trembl_btab_list, args.blast_eval_cutoff)

    if args.blast_kegg_btab_list is not None:
        print("INFO: parsing BLAST (KEGG) evidence")
        parse_kegg_blast_evidence(sources_log_fh, polypeptides, args.blast_kegg_btab_list, args.blast_eval_cutoff)

    if args.blast_uniref100_btab_list is not None:
        print("INFO: parsing BLAST (UniRef100) evidence")
        # connection to the UniRef SQLite3 database
        uniref_db_conn = sqlite3.connect(args.uniref_db)
        uniref_db_curs = uniref_db_conn.cursor()
        parse_uniref100_blast_evidence(sources_log_fh, polypeptides, args.blast_uniref100_btab_list, uniref_db_curs, args.blast_eval_cutoff, 'blast', args.uniref100_fasta)
        uniref_db_curs.close()

    if args.rapsearch_uniref100_btab_list is not None:
        print("INFO: parsing RAPSEARCH2 (UniRef100) evidence")
        # connection to the UniRef SQLite3 database
        uniref_db_conn = sqlite3.connect(args.uniref_db)
        uniref_db_curs = uniref_db_conn.cursor()
        parse_uniref100_blast_evidence(sources_log_fh, polypeptides, args.rapsearch_uniref100_btab_list, uniref_db_curs, args.blast_eval_cutoff, 'rapsearch2', args.uniref100_fasta)
        uniref_db_curs.close()
        
    if args.tmhmm_raw_list is not None:
        print("INFO: parsing TMHMM evidence")
        parse_tmhmm_evidence(sources_log_fh, polypeptides, args.tmhmm_raw_list)

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

    ## There isn't a method in biocodegff3 to add arbitrary key=value pairs.  So we have to cheat here.
    if args.export_organism_names is True:
        if args.output_file:
            append_organism_names_to_gff(args.output_file, polypeptide_blast_org)
        else:
            raise Exception("ERROR: an --output_file must be specified when using the --export_organism_names option.")

    if args.organism_table is not None:
        create_organism_table(args.organism_table, polypeptide_blast_org)


def append_organism_names_to_gff(file_path, poly_orgs):
    # we have to write to a temp file and copy over
    fout = open("{0}.orgtmp".format(file_path), 'wt')
    orgs_found = 0
    last_RNA_id = None
    
    for line in open(file_path):
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) == 9 and cols[2].endswith('RNA'):
            last_RNA_id = gff.column_9_value(cols[8], 'ID')
        if len(cols) == 9 and cols[2] == 'polypeptide':
            if last_RNA_id in poly_orgs:
                cols[8] += ";top_organism_from_blast={0}".format(poly_orgs[last_RNA_id], gff.escape(poly_orgs[last_RNA_id]))
                orgs_found += 1

            fout.write("{0}\n".format("\t".join(cols)) )

        else:
            fout.write("{0}\n".format(line))

    if orgs_found == 0:
        print("WARNING: The --export_organism_names option was passed, but parsing failed to find any organism names at all.  This might be an error.")
            
    ## now move the temp file over the original copy
    fout.close()
    os.rename("{0}.orgtmp".format(file_path), file_path)


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
        if path is not None:
            if not os.path.isfile( path ):
                raise Exception("ERROR: You passed this file but the script failed to find it: {0}".format(path))


def write_gff3_results( f, polypeptides, assemblies, features, genomic_fasta ):
    f.write("##gff-version 3\n")
    
    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            for mRNA in gene.mRNAs():
                # the polypeptide ID needs to be the same as the mRNA one but with the type changed
                polypeptide_id = mRNA.id.replace('mRNA', 'polypeptide')
                #polypeptide_id = mRNA.id
                
                # add polypeptides to the mRNAs
                if polypeptide_id in polypeptides:
                    # the polypeptide ID needs to be the same as the mRNA one but with the type changed
                    new_polypeptide_id = mRNA.id.replace('mRNA', 'polypeptide')
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
    objects, keyed by ID, with bioannotation.FunctionalAnnotation objects attached.
    '''
    seqs = utils.fasta_dict_from_file(fasta_file)

    polypeptides = dict()

    for seq_id in seqs:
        polypeptide = things.Polypeptide(id=seq_id, length=len(seqs[seq_id]['s']), residues=seqs[seq_id]['s'])
        annot = annotation.FunctionalAnnotation(product_name=DEFAULT_PRODUCT_NAME)
        log_fh.write("INFO: {0}: Set initial product name to '{1}'\n".format(seq_id, DEFAULT_PRODUCT_NAME))
        polypeptide.annotation = annot
        
        polypeptides[seq_id] = polypeptide
    
    return polypeptides

def parse_kegg_blast_evidence(log_fh, polypeptides, blast_list, eval_cutoff):
    '''
    Reads a list file of NCBI BLAST evidence against KEGG and a dict of polypeptides,
    populating each with Annotation evidence where appropriate.  Only attaches evidence if
    the product name is the default.

    Currently only considers the top BLAST hit for each query which doesn't have
    'uncharacterized' or hypothetical in the product name.
    '''
    for file in utils.read_list_file(blast_list):
        last_qry_id = None
        
        for line in open(file):
            line = line.rstrip()
            cols = line.split("\t")

            # We're going to ignore any lines which have a few keywords in the name
            # First character left off for initcap reasons
            if 'ncharacterized' in cols[15] or 'ypothetical' in cols[15]:
                continue
            
            this_qry_id = cols[0]

            # skip this line if it doesn't meet the cutoff
            if float(cols[19]) > eval_cutoff:
                continue

            # the BLAST hits are sorted already with the top hit for each query first
            if last_qry_id != this_qry_id:
                annot = polypeptides[this_qry_id].annotation

                # get the accession from the cols[5]
                accession = cols[5]

                # save it, unless the gene product name has already changed from the default
                if annot.product_name == DEFAULT_PRODUCT_NAME:
                    accession = cols[5]

                    # the product field looks like this:
                    # dam; adenine-specific DNA methyltransferase; K06223 DNA adenine methylase [EC:2.1.1.72]
                    # troponin I type 1 (skeletal, slow); K10371 troponin I, slow skeletal muscle
                    if ' [EC' in cols[15] and cols[15].endswith(']'):
                        m = re.search("\; (K\d+)\s+(.+) \[EC\:(.+)\]", cols[15])
                    else:
                        m = re.search("\; (K\d+)\s+(.+)", cols[15])

                    if m:
                        kegg_id = m.group(1)
                        product = m.group(2)
                        
                        if len(m.groups()) == 3:
                            ec_num = m.group(3)
                        else:
                            ec_num = None

                        annot.product_name = product
                        log_fh.write("INFO: {0}: Updated product name to '{1}' based on BLAST hit to KEGG accession '{2}'\n".format(this_qry_id, annot.product_name, accession))

                        if ec_num is not None and ec_num is not '':
                            ec = annotation.ECAnnotation(number=ec_num)
                            annot.add_ec_number(ec)

                        kegg_dbxref = annotation.Dbxref(db='KEGG', identifier=kegg_id)
                        annot.add_dbxref(kegg_dbxref)
                        
                # remember the ID we just saw
                last_qry_id = this_qry_id

def parse_trembl_blast_evidence(polypeptides, blast_list, eval_cutoff):
    '''
    Reads a list file of NCBI BLAST evidence against TrEMBL and a dict of polypeptides,
    populating each with Annotation evidence where appropriate.  Only attaches evidence if
    the product name is the default.

    Currently only considers the top BLAST hit for each query which doesn't have
    'uncharacterized' in the product name.
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


def parse_sprot_blast_evidence( log_fh, polypeptides, blast_org, blast_list, cursor, eval_cutoff, algorithm ):
    '''
    Reads a list file of NCBI BLAST evidence and a dict of polypeptides, populating
    each with Annotation evidence where appropriate.  Only attaches evidence if
    the product name is the default.

    Currently only considers the top BLAST hit for each query.
    '''
    if algorithm not in ['blast', 'rapsearch2']:
        raise Exception("algorithm argument must be either blast or rapsearch2")
    
    for file in utils.read_list_file(blast_list):
        last_qry_id = None
        
        for line in open(file):
            # 0 indexing is faster than startswith()
            if line[0] == '#':
                continue
            
            line = line.rstrip()
            cols = line.split("\t")
            this_qry_id = cols[0]

            if algorithm == 'blast':
                e_value = float(cols[19])
            elif algorithm == 'rapsearch2':
                ## rapsearch2 can actually report values outside of python's double range.  Handle these 
                try:
                    e_value = math.pow(10, float(cols[10]))
                except OverflowError:
                    print("WARN: couldn't handle E-value math on the following line (setting to 0):\n{0}".format(line))
                    e_value = 0

            # skip this line if it doesn't meet the cutoff
            if e_value > eval_cutoff:
                continue

            # the BLAST hits are sorted already with the top hit for each query first
            if last_qry_id != this_qry_id:
                annot = polypeptides[this_qry_id].annotation

                # get the accession from the cols[5]
                #  then process for known accession types
                if algorithm == 'blast':
                    accession = cols[5]
                elif algorithm == 'rapsearch2':
                    accession = cols[1]

                if accession.startswith('sp|'):
                    # pluck the second part out of this:
                    #  sp|Q4PEV8|EIF3M_USTMA
                    accession = accession.split('|')[1]

                assertions = get_uspdb_annot( accession, cursor )
                blast_org[this_qry_id] = assertions['organism']

                # save it, unless the gene product name has already changed from the default
                if annot.product_name == DEFAULT_PRODUCT_NAME:
                    if algorithm == 'blast':
                        # current hack until DB is updated:
                        # some products look like this:
                        #    Coatomer subunit gamma-2 OS=Bos taurus GN=COPG2 PE=2 SV=1
                        # take off everything after the OS=
                        m = re.search("(.+) OS=", cols[15])

                        if m:
                            annot.product_name = m.group(1)
                        else:
                            annot.product_name = cols[15]
                    elif algorithm == 'rapsearch2':
                        annot.product_name = assertions['product']

                    log_fh.write("INFO: {0}: Updated product name to '{1}' based on BLAST hit to SPROT accession '{2}'\n".format(this_qry_id, annot.product_name, accession))

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


def parse_tmhmm_evidence( log_fh, polypeptides, htab_list ):
    '''
    Reads a list of raw TMHMM evidence and a dict of polypeptides, adding annotation
    attributes where possible.

    Notes from the esteemed M Giglio:
    The GO term to use would be GO:0016021 "integral component of membrane"
    Or if you want to be more conservative you could go with GO:0016020 "membrane"
    
    Depends on the evidence. For the prok pipe we are pretty conservative, we require five TMHMM
    domains and then we call it putative integral membrane protein. 

    On ECO - in fact Marcus and I are the developers of ECO.  It is an ontology of evidence types.
    An annotation to an ECO term is used in conjunction with another annotation, like a GO term
    (but many other types of annotation can, and are, used with ECO). It provides additional
    information about the annotation. In fact for GO, the assignment of an evidence term along
    with a GO term is a required part of a GO annotation. (ECO terms are the "evidence codes" in GO.)

    INPUT: Expected TMHMM input (all HTML lines are skipped)
    # CHARM010_V2.mRNA.887 Length: 904
    # CHARM010_V2.mRNA.887 Number of predicted TMHs:  6
    # CHARM010_V2.mRNA.887 Exp number of AAs in TMHs: 133.07638
    # CHARM010_V2.mRNA.887 Exp number, first 60 AAs:  21.83212
    # CHARM010_V2.mRNA.887 Total prob of N-in:        0.99994
    # CHARM010_V2.mRNA.887 POSSIBLE N-term signal sequence
    CHARM010_V2.mRNA.887	TMHMM2.0	inside	     1    11
    CHARM010_V2.mRNA.887	TMHMM2.0	TMhelix	    12    34
    CHARM010_V2.mRNA.887	TMHMM2.0	outside	    35   712
    CHARM010_V2.mRNA.887	TMHMM2.0	TMhelix	   713   735
    CHARM010_V2.mRNA.887	TMHMM2.0	inside	   736   755
    CHARM010_V2.mRNA.887	TMHMM2.0	TMhelix	   756   773
    CHARM010_V2.mRNA.887	TMHMM2.0	outside	   774   782
    CHARM010_V2.mRNA.887	TMHMM2.0	TMhelix	   783   805
    CHARM010_V2.mRNA.887	TMHMM2.0	inside	   806   809
    CHARM010_V2.mRNA.887	TMHMM2.0	TMhelix	   810   832
    CHARM010_V2.mRNA.887	TMHMM2.0	outside	   833   871
    CHARM010_V2.mRNA.887	TMHMM2.0	TMhelix	   872   894
    CHARM010_V2.mRNA.887	TMHMM2.0	inside	   895   904
    '''
    # The number of helices spanning the membrane required before counted as a membrane protein
    MIN_HELICAL_SPANS = 3

    # For successful matches, this is the product name which gets applied
    GENE_PRODUCT_NAME = 'Putative integral membrane protein'
    
    for file in utils.read_list_file(htab_list):
        last_qry_id = None
        current_helix_count = 0
        
        for line in open(file):
            if line.startswith('<'): continue
            m = re.match("# (.+?)\s+Length: \d+", line)

            if m:
                current_id = m.group(1)
                
                # purge previous result
                if current_helix_count >= MIN_HELICAL_SPANS:
                    annot = polypeptides[last_qry_id].annotation

                    if annot.product_name == DEFAULT_PRODUCT_NAME:
                        annot.product_name = GENE_PRODUCT_NAME
                        log_fh.write("INFO: {0}: Updated product name to '{1}' because it had {2} TMHelix domains predicted by TMHMM\n".format(last_qry_id, annot.product_name, current_helix_count))
                    else:
                        log_fh.write("INFO: {0}: TMHMM predicted {1} TMHelix domains but gene product name unchanged because of previous assignment\n".format(last_qry_id, current_helix_count))

                    ## we add the GO terms no matter what
                    annot.add_go_annotation(annotation.GOAnnotation(go_id='0016021'))

                # reset
                last_qry_id = current_id
                current_helix_count = 0
                continue

            cols = line.split()
            if len(cols) == 5 and cols[2] == 'TMhelix':
                current_helix_count += 1

def parse_uniref100_blast_evidence( log_fh, polypeptides, blast_list, cursor, eval_cutoff, algorithm, uniref100_fasta_path ):
    '''
    Reads a list file of NCBI BLAST evidence and a dict of polypeptides, populating
    each with Annotation evidence where appropriate.  Only attaches evidence if
    the product name is the default.

    Currently only considers the top BLAST hit for each query.
    '''
    if algorithm not in ['blast', 'rapsearch2']:
        raise Exception("algorithm argument must be either blast or rapsearch2")

    ## need to load the UniRef100 to TREMBL accession lookup from teh FASTA
    # like UniRef100_K1T359 -> K1T359_9ZZZZ
    uniref2acc = dict()
    print("INFO: parsing UniRef100 FASTA headers for annotation")
    if algorithm == 'rapsearch2':
        for line in open(uniref100_fasta_path):
            if line[0] == '>':
                m = re.match("\>(\S+) (.+) n=.+RepID=(\S+)", line)
                if m:
                    uniref2acc[m.group(1)] = {'acc': m.group(3), 'prod': m.group(2)}
    
    for file in utils.read_list_file(blast_list):
        last_qry_id = None
        
        for line in open(file):
            # 0 indexing is faster than startswith()
            if line[0] == '#':
                continue
            
            line = line.rstrip()
            cols = line.split("\t")
            this_qry_id = cols[0]

            # We're going to ignore any lines which have a few keywords in the name
            # First character left off for initcap reasons
            if algorithm == 'blast':
                skip_products = ['ncharacterized', 'ypothetical', 'enomic scaffold']
                skip = False
                for keyword in skip_products:
                    if keyword in cols[15]:
                        skip = True

                if skip == True:
                    continue

            if algorithm == 'blast':
                e_value = float(cols[19])
            elif algorithm == 'rapsearch2':
                ## rapsearch2 can actually report values outside of python's double range.  Handle these 
                try:
                    e_value = math.pow(10, float(cols[10]))
                except OverflowError:
                    print("WARN: couldn't handle E-value math on the following line (setting to 0):\n{0}".format(line))
                    e_value = 0

            # skip this line if it doesn't meet the cutoff
            if e_value > eval_cutoff:
                continue

            # the BLAST hits are sorted already with the top hit for each query first
            if last_qry_id != this_qry_id:
                annot = polypeptides[this_qry_id].annotation

                # get the accession then process for known accession types
                accession = None

                if algorithm == 'blast':
                    # UniRef100_K1T359 -> K1T359_9ZZZZ
                    m = re.search("RepID\=(\S+)", cols[15])
                    if m:
                        accession = m.group(1)
                    else:
                        raise Exception("ERROR: Unexpected product format in UniRef BLAST results: {0}".format(cols[15]))
                elif algorithm == 'rapsearch2':
                    accession = uniref2acc[cols[1]]['acc']

                assertions = get_uniref_annot( accession, cursor )

                # save it, unless the gene product name has already changed from the default
                if annot.product_name == DEFAULT_PRODUCT_NAME:
                    if algorithm == 'blast':
                        # these hits look like this:
                        #  AD-specific glutamate dehydrogenase n=1 Tax=Ceriporiopsis subvermispora (strain B) RepID=M2RLB9_CERS8
                        m = re.match("(.+) n\=.+", cols[15])
                        if m:
                            annot.product_name = m.group(1)
                        else:
                            raise Exception("ERROR: Unexpected product format in UniRef BLAST results: {0}".format(cols[15]))

                        log_fh.write("INFO: {0}: Updated product name to '{1}' based on BLAST hit to UniRef100 accession '{2}'\n".format(this_qry_id, annot.product_name, accession))
                        
                    elif algorithm == 'rapsearch2':
                        annot.product_name = uniref2acc[cols[1]]['prod']
                        
                # if no EC numbers have been set, they can inherit from this
                if len(annot.ec_numbers) == 0:
                    for ec_annot in get_uniref_ec_nums( accession, cursor ):
                        annot.add_ec_number(ec_annot)

                # if no GO IDs have been set, they can inherit from this
                if len(annot.go_annotations) == 0:
                    for go_annot in get_uniref_go_terms( accession, cursor ):
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

                # do we have a gene symbol for this accession?
                annot.gene_symbol = get_hmmdb_gene_symbol( accession, cursor )

                # do we have an EC number?
                for ec_annot in get_hmmdb_ec_nums( accession, cursor ):
                    annot.add_ec_number(ec_annot)

                ## remember the ID we just saw
                last_qry_id = this_qry_id


def get_uniref_annot( acc, c ):
    """
    This returns a dict of any annotation assertions for a given accession
    in which there will only be one.  Examples: organism name, gene symbol, etc.

    Other methods pull GO terms and EC numbers
    """
    qry = """
          SELECT us.organism, us.symbol
            FROM uniref us
           WHERE us.id = ?
          """
    c.execute(qry, (acc,))
    #print("DEBUG: executing annot query where accession = ({0})".format(acc))

    assertions = { 'organism':None, 'symbol':None }

    for row in c:
        assertions['organism'] = row[0]
        assertions['symbol'] = row[1]
        break

    return assertions


def get_uspdb_annot( acc, c ):
    """
    This returns a dict of any annotation assertions for a given accession
    in which there will only be one.  Examples: organism name, gene symbol, etc.

    Other methods pull GO terms and EC numbers
    """
    qry = """
          SELECT us.organism, us.symbol, us.full_name
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
        assertions['product'] = row[2]
        break

    return assertions


def get_uniref_ec_nums( acc, c ):
    """
    This returns a list of bioannotation:ECAnnotation objects
    """
    qry = """
          SELECT us_ec.ec_num
            FROM uniref_ec us_ec
           WHERE us_ec.id = ?
          """
    c.execute(qry, (acc,))

    ec_annots = list()

    for row in c:
        ec = annotation.ECAnnotation(number=row[0])
        ec_annots.append(ec)

    return ec_annots


def get_uspdb_ec_nums( acc, c ):
    """
    This returns a list of bioannotation:ECAnnotation objects
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

def get_uniref_go_terms( acc, c ):
    """
    This returns a list of bioannotation:GOAnnotation objects
    """
    qry = """
          SELECT us_go.go_id
            FROM uniref_go us_go
           WHERE us_go.id = ?
          """
    c.execute(qry, (acc,))

    go_annots = list()
    
    for row in c:
        go = annotation.GOAnnotation(go_id=row[0], ev_code='ISA', with_from=acc)
        go_annots.append(go)
    
    return go_annots

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

    
def get_hmmdb_ec_nums( acc, c ):
    """
    This returns a list of bioannotation:ECAnnotation objects
    """
    qry = "SELECT he.ec_id FROM hmm_ec he JOIN hmm ON he.hmm_id=hmm.id WHERE hmm.accession = ?"
    c.execute(qry, (acc,))

    ec_annots = list()

    for row in c:
        ec = annotation.ECAnnotation(number=row[0])
        ec_annots.append(ec)
    
    return ec_annots

def get_hmmdb_gene_symbol( acc, c):
    """
    This returns a list of bioannotation:ECAnnotation objects
    """
    qry = "SELECT gene_symbol FROM hmm WHERE accession = ?"
    c.execute(qry, (acc,))

    sym = None

    for row in c:
        sym = row[0]

    return sym


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







