#!/usr/bin/env python3

import argparse
import os
import re
import sqlite3

def main():
    parser = argparse.ArgumentParser( description='Transforms NCBI taxonomy flat-files into a SQLite3 DB.')

    ## categories file columns: 0. top-level category letter, 1. species-level taxid, 2. tax id
    parser.add_argument('-c', '--categories_file', type=str, required=True, help='Path to a categories.dmp file' )

    ## prot taxdump file columns: 0. GenBank GI, 1. tax id
    parser.add_argument('-p', '--prot_taxdump_file', type=str, required=True, help='Path to a gi_taxid_prot.dmp file' )

    ## nucl taxdump file columns: 0. GenBank GI, 1. tax id
    parser.add_argument('-n', '--nucl_taxdump_file', type=str, required=True, help='Path to a gi_taxid_nucl.dmp file' )

    ## names.dmp file columns: 0. tax_id, 1. name_txt, 2. unique name, 3. name class
    #   this script only uses the name class 'scientific name'
    parser.add_argument('-a', '--names_file', type=str, required=True, help='Path to a names.dmp file' )

    ## nodes.dmp file columns: 0. tax_id, 1. parent_tax_id, 2. rank, ... (other columns are ignored)
    parser.add_argument('-d', '--nodes_file', type=str, required=True, help='Path to a nodes.dmp file' )

    ## GbAccList file rows are like: AACY024122252,1,129568417    (acc,version_num,gi)
    #  Example source: ftp://ftp.ncbi.nlm.nih.gov/genbank/livelists/GbAccList.0526.2013.gz
    parser.add_argument('-g', '--gbacclist_file', type=str, required=True, help='Path to a GbAccList.* file' )

    ## RefSeq Accession list file with rows like: 0. tax_id, 1. GI, 2. transcript version, 3. protein version
    #  Example source: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/release59.accession2geneid.gz
    parser.add_argument('-r', '--refseq_acclist_file', type=str, required=True, help='Path to a RefSeq accession2geneid file' )


    ## output file to be written
    parser.add_argument('-o', '--output_sqlite_file', type=str, required=True, help='Path to a the SQLite3 DB to be created' )
    args = parser.parse_args()

    conn = sqlite3.connect( args.output_sqlite_file )
    c = conn.cursor()

    create_tables(c)

    categories = parse_categories( args.categories_file )
    names = parse_names( args.names_file )

    ##-------------------------------------------------------------------------------------
    populate_orgs_table(c, conn, names, categories)

    ##-------------------------------------------------------------------------------------
    populate_nodes_table(c, conn, args.nodes_file)

    ##-------------------------------------------------------------------------------------
    populate_nucl_acc_table(c, conn, args.nucl_taxdump_file)

    ##-------------------------------------------------------------------------------------
    populate_prot_acc_table(c, conn, args.prot_taxdump_file)

    ##-------------------------------------------------------------------------------------
    add_gbk_accessions_to_tables(c, conn, args.gbacclist_file)

    ##-------------------------------------------------------------------------------------
    add_refseq_accessions_to_tables(c, conn, args.refseq_acclist_file)

    c.close()


def create_tables(c):
    c.execute('''
        CREATE TABLE prot_acc
          ( gi INTEGER primary key not null,
            version TEXT,
            tax_id INTEGER )
    ''')

    c.execute('''
        CREATE TABLE nucl_acc
          ( gi INTEGER primary key not null,
            version TEXT,
            tax_id INTEGER )
    ''')

    c.execute('''
        CREATE TABLE orgs
          ( tax_id INTEGER primary key not null,
            cat TEXT, 
            species_tax_id INTEGER,
            scientific_name TEXT )
    ''')

    c.execute('''
        CREATE TABLE nodes
         ( tax_id INTEGER primary key not null,
           parent_tax_id INTEGER,
           rank TEXT )
    ''')



def parse_names( file ):
    names = {}

    if ( not os.path.isfile(file) ):
        raise Exception("Couldn't find file: " + file)

    ofh = open(file, "r")

    for line in ofh:
        columns = line.split('|')

        if ( len(columns) == 5 ):
            sci_name = columns[1].strip()   
            
            if ( columns[3].strip() == 'scientific name' ):
                names[columns[0].strip()] = sci_name
                
    name_count = len(names.keys())
    print("INFO: Loaded", name_count, "taxonomy names")

    return names



def parse_categories( file ):
    categories = {}
    
    if ( not os.path.isfile(file) ):
        raise Exception("Couldn't find file: " + file)

    ofh = open(file, "r")

    for line in ofh:
        columns = line.split()

        ## cat: category letter, spe: species-level tax ID
        if ( len(columns) == 3 ):
            categories[ columns[2] ] = { 'cat': columns[0], 'spe': columns[1] }
    
    category_count = len(categories.keys())
    print("INFO: Loaded", category_count, "categories");

    return categories


def populate_orgs_table( c, conn, names, categories ):
    print("INFO: populating orgs table ...")
    orgs_table_count = 0
    
    for tax_id in names:
        sci_name = ''
        category = ''
        species_id = ''
        
        if tax_id in names:
                sci_name = names[tax_id]
        
        if tax_id in categories:
            category = categories[tax_id]['cat']
            species_id = categories[tax_id]['spe']

        c.execute("""INSERT INTO orgs (tax_id, cat, species_tax_id, scientific_name)
                     VALUES (?,?,?,?)""", \
                  ( tax_id, category, species_id, sci_name ) )
        orgs_table_count += 1

    conn.commit()
    
    print("INFO:\tloaded {0} entries into the orgs table".format(orgs_table_count) )
    print("INFO:\tadding indexes to orgs table")
    c.execute('''CREATE UNIQUE INDEX idx_org_tax_id ON orgs ( tax_id )''')


def populate_nodes_table( c, conn, nodes_file ):
    print("INFO: populating nodes table ...")
    nodes_table_count = 0
    
    if ( not os.path.isfile(nodes_file) ):
        raise Exception("Couldn't find file: " + nodes_file)

    for line in open(nodes_file, "r"):
        columns = line.split('|')
        if len(columns) > 10:
            c.execute("""INSERT INTO nodes (tax_id, parent_tax_id, rank) VALUES (?,?,?)""", \
                      ( columns[0].strip(), columns[1].strip(), columns[2].strip() ) )
            nodes_table_count += 1

    conn.commit()

    print("INFO:\tloaded {0} entries into the nodes table".format(nodes_table_count) )
    print("INFO:\tadding indexes to nodes table")
    c.execute('''CREATE INDEX idx_nodes_parent ON nodes ( parent_tax_id )''')
    c.execute('''CREATE UNIQUE INDEX idx_nodes_tax_id ON nodes ( tax_id )''')



def populate_nucl_acc_table( c, conn, nucl_taxdump_file):
    print("INFO: populating nucl_acc table ...")
    nuc_table_count = 0

    if ( not os.path.isfile(nucl_taxdump_file) ):
        raise Exception("Couldn't find file: " + nucl_taxdump_file)

    for line in open(nucl_taxdump_file, "r"):
        columns = line.split()

        if len(columns) == 2:
            c.execute("""INSERT INTO nucl_acc (gi, tax_id)
                         VALUES (?,?)""", \
                      (columns[0], columns[1]) )
            nuc_table_count += 1

    conn.commit()

    print("INFO:\tloaded {0} entries into the nucl_acc table".format(nuc_table_count) )
    print("INFO:\tadding indexes to nucl_acc table")
    c.execute('''CREATE UNIQUE INDEX idx_nucl_gi ON nucl_acc ( gi )''')



def populate_prot_acc_table( c, conn, prot_taxdump_file):
    print("INFO: populating prot_acc table ...")
    prot_table_count = 0
    
    if ( not os.path.isfile(prot_taxdump_file) ):
        raise Exception("Couldn't find file: " + prot_taxdump_file)

    for line in open(prot_taxdump_file, "r"):
        columns = line.split()

        if len(columns) == 2:
            c.execute("""INSERT INTO prot_acc (gi, tax_id)
                             VALUES (?,?)""", \
                          (columns[0], columns[1]) )
            prot_table_count += 1

    conn.commit()

    print("INFO:\tloaded {0} entries into the prot_acc table".format(prot_table_count) )
    print("INFO:\tadding indexes to prot_acc table")
    c.execute('''CREATE UNIQUE INDEX idx_prot_gi ON prot_acc ( gi )''')



def add_gbk_accessions_to_tables( c, conn, gbacclist_file ):
    print("INFO: adding GenBank accessions to prot_acc and nucl_acc tables ...")
    
    uncommitted_updates = 0

    if ( not os.path.isfile(gbacclist_file) ):
        raise Exception("Couldn't find file: " + gbacclist_file)

    for line in open(gbacclist_file, "r"):
        columns = line.split(',')

        if len(columns) == 3:
            version = columns[0] + '.' + columns[1]
            gi = columns[2].strip()

            ## Protein accessions can be easily distinguished from nucleotide accessions because they have a
            #   three-letter prefix, followed by five digits. The remaining accessions are nucleotide 
            #   accessions, in either a one-letter/five-digit format or a two-letter/six-digit format
            m = re.match( '^[A-Z]{3}[0-9]{5}', version )

            if m:
                #print("INFO: UPDATE prot_acc SET version = {0} WHERE gi = ({1})".format(version, gi))
                c.execute("""UPDATE prot_acc SET version = ? WHERE gi = ?""", (version, gi) )
            else:
                #print("INFO: UPDATE nucl_acc SET version = {0} WHERE gi = ({1})".format(version, gi))
                c.execute("""UPDATE nucl_acc SET version = ? WHERE gi = ?""", (version, gi) )
                
            uncommitted_updates += 1
            if uncommitted_updates == 10000:
                print("INFO: committed 10,000 accession updates ...")
                conn.commit()
                uncommitted_updates = 0

    print("INFO:\tadding indexes for the accessions added for prot and nucs")
    c.execute('''CREATE UNIQUE INDEX idx_prot_version ON prot_acc ( version )''')
    c.execute('''CREATE UNIQUE INDEX idx_nucl_version ON nucl_acc ( version )''')


def add_refseq_accessions_to_tables(c, conn, refseq_acclist_file):
    print("INFO: adding RefSeq accession to prot_acc and nucl_acc tables ...")

    uncommitted_updates = 0

    if ( not os.path.isfile(refseq_acclist_file) ):
        raise Exception("Couldn't find file: " + refseq_acclist_file)

    for line in open(refseq_acclist_file, "r"):
        columns = line.split()

        if len(columns) == 4:
            gi = columns[1].strip()
            nucl_acc = columns[2].strip()
            prot_acc = columns[3].strip()

            if nucl_acc != 'na':
                uncommitted_updates += 1
                #print("DEBUG: UPDATE nucl_acc SET version = {0} WHERE gi = {1}".format(nucl_acc, gi) )
                c.execute("""UPDATE nucl_acc SET version = ? WHERE gi = ?""", (nucl_acc, gi))

            if prot_acc != 'na':
                uncommitted_updates += 1
                #print("DEBUG: UPDATE prot_acc SET version = {0} WHERE gi = {1}".format(prot_acc, gi) )
                c.execute("""UPDATE prot_acc SET version = ? WHERE gi = ?""", (prot_acc, gi))

        if uncommitted_updates == 10000:
            print("INFO: committed 10,000 accession updates ...")
            conn.commit()
            uncommitted_updates = 0



if __name__ == '__main__':
    main()
