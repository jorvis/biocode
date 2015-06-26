#!/usr/bin/env python3

"""

INPUT

The --info_dir is read for .INFO files, each with a format like this:

    ID  purD
    AC  TIGR00877
    DE  phosphoribosylamine--glycine ligase
    AU  Haft DH
    TC  200.00 200.00
    NC  -50.00 -50.00
    AL  clustalw_manual
    IT  equivalog_domain
    EN  phosphoribosylamine--glycine ligase
    GS  purD
    EC  6.3.4.13
    TP  TIGRFAMs
    CC  Alternate name: glycinamide ribonucleotide synthetase (GARS).
    CC  This enzyme appears as a monofunctional protein in prokaryotes but as part of a larger, multidomain protein in eukaryotes.
    DR  EXPERIMENTAL; EGAD|7937|EC4005; Escherichia coli
    DR  ECOCYC; EG10792; purD
    
The only fields I currently look for are these:

    IT  equivalog_domain
    GS  purD
    EC  6.3.4.13
    AC  TIGR00877
    TC  200.00 200.00
    NC  -50.00 -50.00
    CC  This enzyme appears as a monofunctional protein in prokaryotes but as part of a larger, multidomain protein in eukaryotes.  

EC numbers are space-separated in cases with multiple.
    
The --go_link is a tab-delimited text file with this format:

    TIGR00025       GO:0005887      NULL
    TIGR00025       GO:0006810      NULL
    TIGR00025       GO:0043190      NULL
    TIGR00025       GO:0042626      contributes_to
    TIGR00029       GO:0003735      NULL
    TIGR00029       GO:0006412      NULL
    TIGR00029       GO:0000312      NULL

Currently, this script only considers the first two columns.  This file is optional.


OUTPUT

The following tables are are assumed to already exist in the sqlite3 db, created by the
hmmlib_to_sqlite3.py script.


table: hmm
----------
id = auto-int
accession => 'PF02977',
version   => 'PF02977.10',
name => 'S16'
hmm_com_name => 'ribosomal protein S16',
hmm_len => 81,
hmm_comment => "",
trusted_global_cutoff => ?,
trusted_domain_cutoff => ?,
noise_global_cutoff => ?,
noise_domain_cutoff => ?,
gathering_global_cutoff => 35.00,
gathering_domain_cutoff => 35.00,
ec_num => ?,
gene_symbol => ?,
isotype => ?


table: hmm_go
-------------
hmm_id = int
go_id = string (0003735)

table: hmm_ec
-------------
hmm_id = int
ec_id = 1.10.3.5

"""


import argparse
import os
import re
import sqlite3

def main():
    parser = argparse.ArgumentParser( description='Reads supplemental TIGRFAM info and updates an annotation source db file.')

    ## output file to be written
    parser.add_argument('-i', '--info_dir', type=str, required=True, help='Directory with TIGRFAM .INFO files' )
    parser.add_argument('-g', '--go_link', type=str, required=False, help='Path to the TIGRFAMS_GO_LINK file' )
    parser.add_argument('-sd', '--sqlite_db', type=str, required=True, help='Path to an output SQLite3 db to be updated' )
    args = parser.parse_args()

    if not os.path.isfile(args.sqlite_db):
        raise Exception("ERROR: SQLite3 database passed wasn't found.")

    # this creates it if it doesn't already exist
    conn = sqlite3.connect(args.sqlite_db)
    c = conn.cursor()

    ## first cache the ids for TIGRFAM entries already in the database
    #   like d{TIGR00009} = 14840
    tigrfam_id_cache = cache_current_tigrfam_ids(c)
    print("INFO: ID cache created for {0} TIGRFAM entries".format(len(tigrfam_id_cache)))

    for filename in os.listdir(args.info_dir):
        full_path = "{0}/{1}".format(args.info_dir, filename)
        if filename.endswith("INFO") and os.path.isfile(full_path):
            #print("INFO: Processing: {0}".format(full_path))
            process_info_file(full_path, c, tigrfam_id_cache)
            conn.commit()

    if args.go_link is not None:
        for line in open(args.go_link):
            line = line.rstrip()
            cols = line.split()

            if len(cols) == 3:
                if cols[0] in tigrfam_id_cache:
                    m = re.match("GO\:(\d+)", cols[1])
                    if m:
                        c.execute("""INSERT INTO hmm_go (hmm_id, go_id) VALUES (?, ?)""", (tigrfam_id_cache[cols[0]], m.group(1)) )

    conn.commit()
    c.close()


def process_info_file(filename, curs, id_cache):
    accession = None
    gene_symbol = None
    ec_numbers = list()
    
    for line in open(filename):
        try:
            line = line.rstrip()
        except UnicodeDecodeError:
            continue
        
        if line.startswith('AC  '):
            m = re.match('AC\s+(\S+)', line)
            if m:
                accession = m.group(1)
            else:
                raise Exception("ERROR: Unexpected format of AC line in file {0}".format(filename))
        elif line.startswith('EC  '):
            m = re.match('EC  (.+)', line)
            if m:
                for ec in m.group(1).split(' '):
                    ec_numbers.append(ec)
            else:
                raise Exception("ERROR: Unexpected format of EC line in file {0}".format(filename))
        elif line.startswith('GS  '):
            m = re.match('GS\s+(\S+)', line)
            if m:
                gene_symbol = m.group(1)
            else:
                raise Exception("ERROR: Unexpected format of GS line in file {0}".format(filename))

    #print("File: {0}\n\taccession: {1}\n\tgene_symbol: {2}\n\tec_numbers: {3}".format(filename, accession, gene_symbol, ec_numbers))
    if accession in id_cache:
        if gene_symbol is not None:
            curs.execute("""UPDATE hmm SET gene_symbol = ? WHERE accession = ?""", (gene_symbol, accession) )

        for ec_number in ec_numbers:
            #print("INFO: Inserting into hmm_ec: {0}, {1}".format(id_cache[accession], ec_number))
            curs.execute("""INSERT INTO hmm_ec (hmm_id, ec_id) VALUES (?, ?)""", (id_cache[accession], ec_number) )
    

def cache_current_tigrfam_ids(curs):
    cache = dict()
    for row in curs.execute("""SELECT hmm.id, hmm.accession from hmm where hmm.accession like 'TIGR%'""" ):
        cache[row[1]] = row[0]

    return cache




if __name__ == '__main__':
    main()







