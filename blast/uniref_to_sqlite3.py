#!/usr/bin/env python3

"""
Reads the uniref100 flat file and creates a SQLite3 database of commonly-accessed attributes
for each entry.

INPUT

When you BLAST against a database like UniRef100 or UniRef90 you'll get accessions like:

  UniRef100_A8PWR2      RepID=A8PWR2_MALGO
  UniRef100_M9MHS5      RepID=M9MHS5_PSEA3

The first is the accession and later in the description contains the link to the RepID.

  A8PWR2_MALGO

Entries for these are found in the uniprot_trembl.dat.gz file.

We're primarily looking for the following attributes:
  - id
  - accession
  - gene product names
  - gene symbol
  - GO terms
  - EC numbers

These look like this in the file:

ID   A8PWR1_MALGO            Unreviewed;       963 AA.
AC   Q9NRN7; B2R6D1; Q9C068; Q9P0Q3; Q9UG80; Q9Y389;

DE   SubName: Full=Putative uncharacterized protein;

DR   GO; GO:0000062; F:fatty-acyl-CoA binding; IEA:InterPro.
DR   GO; GO:0008289; F:lipid binding; IEA:UniProtKB-KW.
DR   GO; GO:0005737; C:cytoplasm; IEA:UniProtKB-SubCell.
DR   GO; GO:0005634; C:nucleus; IEA:Compara.

OS   Homo sapiens (Human).

GN   Name=amoA;

DE              EC=6.3.4.3;
DE              EC=1.5.1.5;
DE              EC=3.5.4.9;
DE              EC=6.3.4.3;

OUTPUT

The following tables are created in the SQLite3 db (which is created if it doesn't
already exist) (these are fake example data, and there are a lot of 1:many 
relationships here):

table: entry
----------
id = 001R_FRG3G
full_name = 11S globulin subunit beta
organism = Frog virus 3 (isolate Goorha)
symbol = FV3-001R

entry_acc
-----------------
id = 001R_FRG3G
accession = Q6GZX4
res_length = 121
is_characterized = 1

entry_go
----------------
id = 001R_FRG3G
go_id = string (0005634)

entry_ec
----------------
id = 001R_FRG3G
ec_num = 6.3.4.3

"""


import argparse
import gzip
import os
import re
import sqlite3

pfam2go = dict()

def main():
    parser = argparse.ArgumentParser( description='Reads a uniprot_trembl.dat file and creates a SQLite3 database of commonly-accessed attributes for each accession.')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the uniprot_trembl.dat file.  See INPUT section for details' )
    parser.add_argument('-o', '--output_db', type=str, required=True, help='Path to an output SQLite3 db to be created' )
    args = parser.parse_args()

    # this creates it if it doesn't already exist
    conn = sqlite3.connect(args.output_db)
    curs = conn.cursor()

    print("INFO: Creating tables ...")
    create_tables( curs )
    conn.commit()

    id = None
    accs = list()
    full_name = None
    organism = None
    symbol = None
    go_ids = list()
    ec_nums = list()
    record_count = 0
    res_length = 0
    is_characterized = 0
    is_compressed = False

    characterized_go_codes = ['EXP', 'IDA', 'IMP', 'IGI', 'IPI', 'IEP']

    if args.input.endswith('.gz'):
        ifh = gzip.open(args.input, 'rb')
        is_compressed = True
    else:
        ifh = open(args.input)
    
    print("INFO: Parsing DAT file ...")
    for line in ifh:
        if is_compressed:
            line = line.decode()
        
        line = line.rstrip()
        
        # is this the end of an entry?
        if re.match( "^//", line ):
            curs.execute("INSERT INTO entry (id, full_name, organism, symbol) values (?,?,?,?)", (id, full_name, organism, symbol))

            # nothing with hypothetical in the name can be characterized
            m = re.search('hypothetical', full_name, re.IGNORECASE)
            if m:
                is_characterized = 0
            
            for acc in accs:
                curs.execute("INSERT INTO entry_acc (id, accession, res_length, is_characterized) values (?,?,?,?)", (id, acc, res_length, is_characterized))

            for go_id in go_ids:
                curs.execute("INSERT INTO entry_go (id, go_id) values (?,?)", (id, go_id))

            for ec_num in ec_nums:
                curs.execute("INSERT INTO entry_ec (id, ec_num) values (?,?)", (id, ec_num))
            
            # reset
            id = None
            accs = list()
            full_name = None
            organism = None
            symbol = None
            go_ids = list()
            ec_nums = list()
            res_length = 0
            is_characterized = 0

            record_count += 1

            if record_count % 1000 == 0:
                print("{0} records processed ...".format(record_count))
                conn.commit()

        elif line.startswith("ID"):
            m = re.match("ID\s+(\S+).+?(\d+) AA\.", line)
            if m:
                id = m.group(1)
                res_length = m.group(2)
            else:
                raise Exception("ERROR: the following ID line didn't match the expected convention:\n {0}".format(line))
        elif line.startswith("AC"):
            ac_parts = line.split()
            for part in ac_parts:
                if part == 'AC':
                    continue
                else:
                    part = part.rstrip(';')
                    accs.append(part)
            
        elif line.startswith("DE"):
            m = re.match("DE\s+SubName: Full=(.+)", line.rstrip(';'))
            if m and full_name is None:
                full_name = m.group(1)
            else:
                m = re.match("DE\s+RecName: Full=(.+)", line.rstrip(';'))
                if m and full_name is None:
                    full_name = m.group(1)
                else:
                    m = re.search("EC=(\S+)", line.rstrip(';'))
                    if m:
                        ec_nums.append(m.group(1))

        elif line.startswith("DR"):
            m = re.search("GO:(\d+)\; .+?\; ([A-Z]+):\S+\.$", line)
            if m:
                go_ids.append(m.group(1))

                if m.group(2) in characterized_go_codes:
                    is_characterized = 1
                

        elif line.startswith("OS"):
            m = re.match("OS\s+(.+)\.$", line)
            if m:
                organism = m.group(1)
                
        elif line.startswith("GN"):
            m = re.search("^GN\s+Name=(.+?)\;", line)
            if m:
                symbol = m.group(1)
        
            
    conn.commit()

    create_indexes(curs)
    conn.commit()
    curs.close()
    print("INFO: Complete.")
    


def create_indexes( cursor ):
    print("INFO: Creating indexes ...")
    cursor.execute("CREATE INDEX idx_col_us_id  ON entry (id)")

    cursor.execute("CREATE INDEX idx_col_usa_id  ON entry_acc (id)")
    cursor.execute("CREATE INDEX idx_col_usa_acc ON entry_acc (accession)")
    
    cursor.execute("CREATE INDEX idx_col_usg_id ON entry_go (id)")
    cursor.execute("CREATE INDEX idx_col_usg_go ON entry_go (go_id)")

    cursor.execute("CREATE INDEX idx_col_use_id ON entry_ec (id)")
    cursor.execute("CREATE INDEX idx_col_use_ec ON entry_ec (ec_num)")


def create_tables( cursor ):
    cursor.execute("""
              CREATE TABLE entry (
                 id                text primary key,
                 full_name         text,
                 organism          text,
                 symbol            text
              )
    """)
    
    cursor.execute("""
              CREATE TABLE entry_acc (
                 id         text not NULL,
                 accession  text not NULL,
                 res_length integer,
                 is_characterized integer DEFAULT 0
              )
    """)
    
    cursor.execute("""
              CREATE TABLE entry_go (
                 id     text not NULL,
                 go_id  text not NULL
              )
    """)

    cursor.execute("""
              CREATE TABLE entry_ec (
                 id     text not NULL,
                 ec_num text not NULL
              )
    """)


if __name__ == '__main__':
    main()







