#!/usr/bin/env python3

"""
Reads an HMM file and creates a SQLite3 database of commonly-accessed attributes
for each HMM.

Several scripts perform lookups on HMM attributes such as annotation, trusted cutoff, etc and need
a rapid way to access this information.  This script reads an HMM library and generates an SQLite3
file of the most commonly accessed attributes.

INPUT

As described above, the input can be a library of HMMs or a single HMM file.  Records are separated
with a // symbol.  An example record:

    HMMER2.0
    NAME  S16
    ACC   TIGR00002
    DESC  S16: ribosomal protein S16
    LENG  81
    ALPH  Amino
    RF    no
    CS    no
    MAP   yes
    COM   hmmbuild --pam blosum62 --pamwgt 20 ogc_results//853/17288/17288_1_3.HMM ogc_results//853/17288/17288_1.msf
    COM   hmmcalibrate --cpu 1 --num 5000 --histfile ogc_results//853/17288/17288_1_3.hist ogc_results//853/17288/17288_1_3.HMM
    NSEQ  20
    DATE  Mon Feb 11 17:46:41 2002
    CKSUM 95
    GA    35.00 35.00
    TC    35.00 35.00
    NC    10.00 10.00
    XT      -8455     -4  -1000  -1000  -8455     -4  -8455     -4 
    NULT      -4  -8455
    NULE     595  -1558     85    338   -294    453  -1158    197    249    902  -1085   -142    -21   -313     45    531    201    384  -1998   -644 
    EVD   -40.960884   0.269394
    HMM        A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y    
             m->m   m->i   m->d   i->m   i->i   d->m   d->d   b->m   m->e
              -71      *  -4389
         1  -1123  -1482  -2599  -2123  -1302  -2623  -2482    844  -2004   1300   -331  -2482  -2123  -2004  -2331  -1982  -1123   2752  -2331  -1482     1
         -   -149   -500    233     43   -381    399    106   -626    210   -466   -720    275    394     45     96    359    117   -369   -294   -249 
         -     -3  -9476 -10518   -894  -1115   -701  -1378    -71      *

        ... [ and so on for all positions ] ...

             81   -439  -1642  -2419  -2205   2226  -2269  -1903   1024  -2082    446   1484  -2214  -2388  -1909  -2082  -1628    553   -751   3736   -339    91
         -      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      *      * 
         -      *      *      *      *      *      *      *      *      0 
//


OUTPUT

The following tables are created in the SQLite3 db (which is created if it doesn't already exist):

table: hmm
----------
id
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

Values for isotype should currently be limited to:

    equivalog
    subfamily
    hypoth_equivalog
    equivalog_domain
    subfamily_domain
    domain
    superfamily
    paralog
    hypoth_equivalog_domain
    paralog_domain
    exception
    repeat
    signature

The only isotype this script will set is 'pfam' if the accession starts with PF.  All
others are added by other scripts, such as tigrfam_info_to_mldbm.pl

Any other processes or scripts can add to this data structure to possibly improve annotation.  The
tigrfam_info_to_mldbm script, for example, adds ec_num, gene_symbol and isotype to the data
structure for TIGRFAM entries.

"""


import argparse
import os
import re
import sqlite3


def main():
    parser = argparse.ArgumentParser( description='Reads an HMM file and creates a SQLite3 database of commonly-accessed attributes for each HMM.')

    ## output file to be written
    parser.add_argument('-i', '--input_hmm', type=str, required=True, help='Path to an HMM collection.  See INPUT section for details' )
    parser.add_argument('-o', '--output_db', type=str, required=True, help='Path to an output SQLite3 db to be created' )
    args = parser.parse_args()

    # this creates it if it doesn't already exist
    conn = sqlite3.connect(args.output_db)
    c = conn.cursor()

    create_tables( c )
    conn.commit()

    # store attributes for each HMM as we encounter them.  These are inserted into the
    #  database at the end of each record
    atts = dict()

    for line in open(args.input_hmm):
        line = line.rstrip()
        
        # is this the end of an entry?
        if re.match( "^//", line ):
            add_hmm(atts, c)
            atts = dict()
        else:
            m = re.match("^([A-Z]+)\s+(.*)$", line)
            if m:
                key = m.group(1)
                atts[key] = m.group(2)

    conn.commit()
    c.close()
    


def add_hmm( atts, cur ):
    # the name attribute is a fall-back if acc isn't defined.
    if 'ACC' not in atts or len(atts['ACC']) == 0:
        atts['ACC'] = atts['NAME']

    print("INFO: adding hmm: {0}".format(atts['ACC']))
        
    trusted_global = None
    trusted_domain = None
    noise_global = None
    noise_domain = None
    gathering_global = None
    gathering_domain = None

    # set some defaults
    for key in ('CC', 'EC', 'GS', 'IT'):
        if key not in atts:
            atts[key] = None

    if 'TC' in atts:
        m = re.match("([0-9\.\-]*)\s+([0-9\.\-]*)", atts['TC'])
        if m:
            trusted_global = m.group(1)
            trusted_domain = m.group(2)
        
    if 'NC' in atts:
        m = re.match("([0-9\.\-]+)\s+([0-9\.\-]+)", atts['NC'])
        if m:
            noise_global = m.group(1)
            noise_domain = m.group(2)

    if 'GA' in atts:
        m = re.match("([0-9\.\-]+)\s+([0-9\.\-]+)", atts['GA'])
        if m:
            gathering_global = m.group(1)
            gathering_domain = m.group(2)

    # handle the product name.  if it contains NAME at the beginning, remove that
    product_name = atts['DESC']
    product_name = product_name.rstrip( "{0}: ".format(atts['NAME']) )

    # if this looks like a PFAM accession (PF03372.14), set the isotype to 'pfam'
    if re.match("^PF\d+\.\d+", atts['ACC']):
        atts['IT'] = 'pfam'

    columns = "name, hmm_com_name, hmm_len, hmm_comment, trusted_global_cutoff, trusted_domain_cutoff, " \
              "noise_global_cutoff, noise_domain_cutoff, gathering_global_cutoff, gathering_domain_cutoff, " \
              "ec_num, gene_symbol, isotype"
    cur.execute("INSERT INTO hmm ({0}) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)".format(columns), \
               (atts['NAME'], product_name, atts['LENG'], atts['CC'], trusted_global, trusted_domain, \
                noise_global, noise_domain, gathering_global, gathering_domain, \
                atts['EC'], atts['GS'], atts['IT'] ))


def create_tables( cursor ):
    cursor.execute("""
              CREATE TABLE hmm (
                 id                integer primary key autoincrement,
                 name              text,
                 hmm_com_name      text,
                 hmm_len           int,
                 hmm_comment       text,
                 trusted_global_cutoff    float,
                 trusted_domain_cutoff   float,
                 noise_global_cutoff      float,
                 noise_domain_cutoff     float,
                 gathering_global_cutoff  float,
                 gathering_domain_cutoff float,
                 ec_num            text,
                 gene_symbol       text,
                 isotype           text
              )
    """)

    cursor.execute("""
              CREATE TABLE hmm_go (
                 id     integer primary key autoincrement,
                 hmm_id int not NULL,
                 go_id  text
              )
    """)

    cursor.execute("""
              CREATE TABLE hmm_ec (
                 id     integer primary key autoincrement,
                 hmm_id int not NULL,
                 ec_id  text
              )
    """)




if __name__ == '__main__':
    main()







