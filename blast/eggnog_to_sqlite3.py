#!/usr/bin/env python3

"""
Reads files from an EggNOG release and stores them as an indexed SQLite3 file for
rapid lookups by protein ID.  The use case for development here is our functional
annotation pipeline.

INPUT

COG.members.txt
---------------
COG1110 579137.Metvu_1264       1       865
COG1110 573064.Mefer_0049       6       869
COG1110 563041.HPG27_964        428     531
COG1110 550538.SG1400   318     571

COG.funccat.txt
---------------
COG1110 L
COG5155 DO
COG1119 P

fun.txt
-------
CELLULAR PROCESSES AND SIGNALING
 [D] Cell cycle control, cell division, chromosome partitioning 
 [Y] Nuclear structure 
 [V] Defense mechanisms

COG.description.txt
-------------------
COG1110 Reverse gyrase
COG0393 Uncharacterized conserved protein
COG4055 Methyl coenzyme M reductase, subunit D
COG1943 Transposase and inactivated derivatives


OUTPUT

The following tables are created in the SQLite3 db (which is created if it doesn't
already exist) 

cog
---
cog_id : COG1110
symbol: DO
desc : Reverse gyrase

funccat
------
symbol: D
desc  : Cell cycle control, cell division, chromosome partitioning

member
-------------
cog_id : COG1110
accession : 579137.Metvu_1264

Example command:
$ /home/jorvis/git/biocode/blast/eggnog_to_sqlite3.py -o eggnog.v3.db -f /local/db/repository/embl/eggNOG/20140314_151611/fun.txt -fc /local/db/repository/embl/eggNOG/20140314_151611/all.funccat/COG.funccat.txt -cd /local/db/repository/embl/eggNOG/20140314_151611/all.description/COG.description.txt -cm /local/db/repository/embl/eggNOG/20140314_151611/all.members/COG.members.txt

"""


import argparse
import os
import re
import sqlite3
import sys

pfam2go = dict()

def main():
    parser = argparse.ArgumentParser( description='Reads eggNOG release files and creates a SQLite3 database of commonly-accessed attributes for each accession.')

    ## output file to be written
    parser.add_argument('-f', '--fun', type=str, required=True, help='Path to the fun.txt file.  See INPUT section for details' )
    parser.add_argument('-fc', '--funccat', type=str, required=True, help='Path to the COG.funccat.txt file.  See INPUT section for details' )
    parser.add_argument('-cd', '--cog_description', type=str, required=True, help='Path to the COG.description.txt file.  See INPUT section for details' )
    parser.add_argument('-cm', '--cog_members', type=str, required=True, help='Path to the COG.members.txt file.  See INPUT section for details' )
    parser.add_argument('-o', '--output_db', type=str, required=True, help='Path to an output SQLite3 db to be created' )
    args = parser.parse_args()

    # key = symbol, value = description
    fun = dict()
    for line in open(args.fun):
        line = line.rstrip()
        m = re.match("\s*\[(.+?)\] (.+)", line)
        if m:
            fun[m.group(1)] = m.group(2)

    print("INFO: {0} functional categories loaded into memory".format(len(fun)))

    # key = COG_id, value = Full description
    cog_descriptions = dict()
    for line in open(args.cog_description):
        line = line.rstrip()
        m = re.match("(.+?)\s+(.+)", line)
        if m:
            cog_descriptions[m.group(1)] = m.group(2)

    print("INFO: {0} COG descriptions loaded into memory".format(len(cog_descriptions)))

    # this creates it if it doesn't already exist
    conn = sqlite3.connect(args.output_db)
    curs = conn.cursor()

    print("INFO: Creating tables ...")
    create_tables( curs )
    conn.commit()

    print("INFO: Loading tables ...")
    
    for sym in fun:
        curs.execute("INSERT INTO funccat (symbol, desc) values (?,?)", (sym, fun[sym]))

    funccat_lines = 0
    for line in open(args.funccat):
        line = line.rstrip()
        m = re.match("(.+?)\s+(.+)", line)
        if m:
            funccat_lines += 1
            cog_id = m.group(1)
            if cog_id in cog_descriptions:
                curs.execute("INSERT INTO cog (cog_id, symbol, desc) values (?,?,?)", (cog_id, m.group(2), cog_descriptions[cog_id]))
            else:
                raise Exception("ERROR: COG {0} found in funccat file but not in descriptions file".format(m.group(1)))
        else:
            print("This line didn't match: {0}".format(line))

    for line in open(args.cog_members):
        line = line.rstrip()
        cols = line.split("\t")
        if cols[0] in cog_descriptions:
            curs.execute("INSERT INTO member (cog_id, accession) values (?,?)", (cols[0], cols[1]))
        else:
            raise Exception("ERROR: COG {0} found in members file but not in descriptions file".format(cols[0]))
            

    # sanity check
    if funccat_lines != len(cog_descriptions):
        print("WARNING: The number of lines parsed from the funccat ({0}) and description file ({1}) weren't the same and were expected to be.".format(funccat_lines, len(cog_descriptions)))
            
    conn.commit()

    print("INFO: Creating indexes ...")
    create_indexes(curs)
    conn.commit()
    curs.close()
    print("INFO: Complete.")
    


def create_indexes( cursor ):
    cursor.execute("CREATE INDEX idx_member_cog_id ON member (cog_id)")
    cursor.execute("CREATE INDEX idx_member_accession ON member (accession)")
    cursor.execute("CREATE INDEX idx_cog_cog_id ON cog (cog_id)")



def create_tables( cursor ):
    cursor.execute("""
        CREATE TABLE member (
            cog_id      text not NULL,
            accession   text not NULL
        )
    """)
    
    cursor.execute("""
              CREATE TABLE funccat (
                 symbol   text not NULL,
                 desc     text not NULL
              )
    """)
    
    cursor.execute("""
              CREATE TABLE cog (
                 cog_id         text not NULL,
                 symbol         text not NULL,
                 desc           text not NULL
              )
    """)


if __name__ == '__main__':
    main()







