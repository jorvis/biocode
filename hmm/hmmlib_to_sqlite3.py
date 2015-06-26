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

Another optional input is the pfam2go data file:
ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/external2go/pfam2go


OUTPUT

The following tables are created in the SQLite3 db (which is created if it doesn't
already exist) (these are fake example data):

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

pfam2go = dict()

def main():
    parser = argparse.ArgumentParser( description='Reads an HMM file and creates a SQLite3 database of commonly-accessed attributes for each HMM.')

    ## output file to be written
    parser.add_argument('-i', '--input_hmm', type=str, required=True, help='Path to an HMM collection.  See INPUT section for details' )
    parser.add_argument('-p', '--pfam2go', type=str, required=False, help='EBI pfam2go file' )
    parser.add_argument('-o', '--output_db', type=str, required=True, help='Path to an output SQLite3 db to be created' )
    args = parser.parse_args()

    if args.pfam2go is None:
        pfam2go = None
    else:
        pfam2go = load_pfam2go(args.pfam2go)
        print("DEBUG: loaded {0} pfam2go accessions".format(len(pfam2go)))
        
    if os.path.isfile(args.output_db):
        db_exists = True
    else:
        db_exists = False

    # this creates it if it doesn't already exist
    conn = sqlite3.connect(args.output_db)
    c = conn.cursor()

    if not db_exists:
        create_tables( c )
        conn.commit()

    # store attributes for each HMM as we encounter them.  These are inserted into the
    #  database at the end of each record
    atts = dict()

    for line in open(args.input_hmm):
        line = line.rstrip()
        
        # is this the end of an entry?
        if re.match( "^//", line ):
            add_hmm(atts, c, pfam2go)
            atts = dict()
        else:
            m = re.match("^([A-Z]+)\s+(.*)$", line)
            if m:
                key = m.group(1)
                atts[key] = m.group(2)

    conn.commit()

    create_indexes(c)
    
    c.close()
    

def load_pfam2go( path ):
    pf_acc = dict()
    
    for line in open(path):
        line = line.rstrip()
        m = re.match("Pfam:(PF\d+) .+ GO\:(\d+)", line)
        if m:
            acc = m.group(1)
            go_id = m.group(2)

            if acc not in pf_acc:
                pf_acc[acc] = list()

            pf_acc[acc].append(go_id)

    return pf_acc


def add_hmm( atts, cur, pf2g ):
    # the name attribute is a fall-back if acc isn't defined.
    if 'ACC' not in atts or len(atts['ACC']) == 0:
        atts['ACC'] = atts['NAME']

    accession = atts['ACC']
    version = None
    ec_number = None

    # if this is a PFAM accession, split it into accession and version and
    #  set the isotype to 'pfam'
    m = re.match("^(PF\d+)\.\d+", accession)
    if m:
        version = accession
        accession = m.group(1)
        atts['IT'] = 'pfam'

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
    m = re.match("{0}\: (.+)".format(atts['NAME']), product_name)
    if m:
        product_name = m.group(1)

    # if the name ends in ' (EC 4.6.1.3)' capture that, then take it off
    m = re.match("(.+) \(EC (\d+\.\d+\.\d+\.\d+)\)", product_name)
    if m:
        product_name = m.group(1)
        ec_number = m.group(2)

    # TIGRFAMS keep skipping the GS field and appending it to the beginning of the DESC
    #  example:  DESC  argG: argininosuccinate synthase
    # Let's make the rule that if the DESC starts with 7 or less consecutive non-whitespace
    # characters, then a ": ", we'll strip all that off as the gene symbol.  TIGRFAM accessions only.
    # If it's longer than 7 it appears they are appending what used to be the NAME, so strip it off.
    if accession.startswith('TIGR'):
        m = re.match("^(\S+)\: (.+)", product_name)
        if m:
            product_name = m.group(2)
            
            if atts['GS'] is None:
                if len(m.group(1)) <= 7:
                    atts['GS'] = m.group(1)
                
    print("INFO: adding hmm: {0} - {1} - {2}".format(accession, atts['NAME'], product_name))

    columns = "accession, version, name, hmm_com_name, hmm_len, hmm_comment, trusted_global_cutoff, trusted_domain_cutoff, " \
              "noise_global_cutoff, noise_domain_cutoff, gathering_global_cutoff, gathering_domain_cutoff, " \
              "ec_num, gene_symbol, isotype"
    cur.execute("INSERT INTO hmm ({0}) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)".format(columns), \
               (accession, version, atts['NAME'], product_name, atts['LENG'], atts['CC'], trusted_global, trusted_domain, \
                noise_global, noise_domain, gathering_global, gathering_domain, \
                atts['EC'], atts['GS'], atts['IT'] ))

    hmm_id = cur.lastrowid

    # if pfam2go was passed, does this accession have any associations?
    if pf2g is not None and accession in pf2g:
        for go_id in pf2g[ accession ]:
            cur.execute("INSERT INTO hmm_go (hmm_id, go_id) values (?,?)", (hmm_id, go_id))

    if ec_number is not None:
        cur.execute("INSERT INTO hmm_ec (hmm_id, ec_id) values (?,?)", (hmm_id, ec_number))


def create_tables( cursor ):
    cursor.execute("""
              CREATE TABLE hmm (
                 id                integer primary key autoincrement,
                 accession         text,
                 version           text,
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

def create_indexes( cursor ):
    cursor.execute("CREATE INDEX idx_hmm__accession ON hmm(accession)")
    cursor.execute("CREATE INDEX idx_hmm_ec__hmm_id ON hmm_ec(hmm_id)")
    cursor.execute("CREATE INDEX idx_hmm_go__hmm_id ON hmm_go(hmm_id)")


if __name__ == '__main__':
    main()







