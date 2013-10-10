#!/usr/bin/env python3

"""
Reads the uniprot_sprot flat file and creates a SQLite3 database of commonly-accessed attributes
for each entry.

INPUT

This file (uncompressed):
ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz

We're primarily looking for the following attributes:
  - id
  - accessions for that id
  - gene product names
  - GO terms
  - EC numbers

These look like this in the file:

ID   ADPPT_HUMAN             Reviewed;         309 AA.
AC   Q9NRN7; B2R6D1; Q9C068; Q9P0Q3; Q9UG80; Q9Y389;

DE   RecName: Full=11S globulin subunit beta;
DE   RecName: Full=Uncharacterized protein 120L;
DE   RecName: Full=Putative RING finger protein 121R;

DR   GO; GO:0000062; F:fatty-acyl-CoA binding; IEA:InterPro.
DR   GO; GO:0008289; F:lipid binding; IEA:UniProtKB-KW.
DR   GO; GO:0005737; C:cytoplasm; IEA:UniProtKB-SubCell.
DR   GO; GO:0005634; C:nucleus; IEA:Compara.

DE            EC=6.3.4.3;
DE              EC=1.5.1.5;
DE              EC=3.5.4.9;
DE              EC=6.3.4.3;


OUTPUT

The following tables are created in the SQLite3 db (which is created if it doesn't
already exist) (these are fake example data, and there are a lot of 1:many 
relationships here):

table: uniprot_sprot
----------
id = 001R_FRG3G
full_name = 11S globulin subunit beta

uniprot_sprot_acc
-----------------
id = 001R_FRG3G
accession = Q6GZX4

uniprot_sprot_go
----------------
id = 001R_FRG3G
go_id = string (0005634)

uniprot_sprot_ec
----------------
id = 001R_FRG3G
ec_num = 6.3.4.3

"""


import argparse
import os
import re
import sqlite3

pfam2go = dict()

def main():
    parser = argparse.ArgumentParser( description='Reads a uniprot_sprot.dat file and creates a SQLite3 database of commonly-accessed attributes for each accession.')

    ## output file to be written
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the uniprot_sprot.dat file.  See INPUT section for details' )
    parser.add_argument('-o', '--output_db', type=str, required=True, help='Path to an output SQLite3 db to be created' )
    args = parser.parse_args()

    # this creates it if it doesn't already exist
    conn = sqlite3.connect(args.output_db)
    curs = conn.cursor()

    create_tables( curs )
    conn.commit()

    id = None
    accs = list()
    full_name = None
    go_ids = list()
    ec_nums = list()
    

    for line in open(args.input):
        line = line.rstrip()
        
        # is this the end of an entry?
        if re.match( "^//", line ):
            # save
            curs.execute("INSERT INTO uniprot_sprot (id, full_name) values (?,?)", (id, full_name))

            for acc in accs:
                curs.execute("INSERT INTO uniprot_sprot_acc (id, accession) values (?,?)", (id, acc))

            for go_id in go_ids:
                curs.execute("INSERT INTO uniprot_sprot_go (id, go_id) values (?,?)", (id, go_id))

            for ec_num in ec_nums:
                curs.execute("INSERT INTO uniprot_sprot_ec (id, ec_num) values (?,?)", (id, ec_num))
            
            # reset
            id = None
            accs = list()
            full_name = None
            go_ids = list()
            ec_nums = list()
            
        elif line.startswith("ID"):
            id = line.split()[1]
        elif line.startswith("AC"):
            ac_parts = line.split()
            for part in ac_parts:
                if part == 'AC':
                    continue
                else:
                    part = part.rstrip(';')
                    accs.append(part)
            
        elif line.startswith("DE"):
            m = re.match("DE\s+RecName: Full=(.+)", line.rstrip(';'))
            if m:
                full_name = m.group(1)
            else:
                m = re.search("EC=(\S+)", line.rstrip(';'))
                if m:
                    ec_nums.append(m.group(1))
        elif line.startswith("DR"):
            m = re.search("GO:(\d+)", line)
            if m:
                go_ids.append(m.group(1))
            

    conn.commit()
    curs.close()
    



def create_tables( cursor ):
    cursor.execute("""
              CREATE TABLE uniprot_sprot (
                 id                text primary key,
                 full_name         text
              )
    """)
    
    cursor.execute("""
              CREATE TABLE uniprot_sprot_acc (
                 id         text not NULL,
                 accession  text not NULL
              )
    """)
    
    cursor.execute("""
              CREATE TABLE uniprot_sprot_go (
                 id     text not NULL,
                 go_id  text not NULL
              )
    """)

    cursor.execute("""
              CREATE TABLE uniprot_sprot_ec (
                 id     text not NULL,
                 ec_num text not NULL
              )
    """)




if __name__ == '__main__':
    main()







