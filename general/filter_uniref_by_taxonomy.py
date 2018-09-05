#!/usr/bin/env python3

"""
Filter the UniRef100 FASTA release by taxonomy.  Assumes you have already downloaded
this file (or other UniRef variant):

   ftp://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz

This utility relies on the taxadb python module:
  
   https://github.com/HadrienG/taxadb
   or
   pip3 install taxadb

The taxadb module requires you to first build a taxonomy database, which this script
uses.  If it's your first time, you'll need to install the taxadb Python module, then
build a database like this:

   $ taxadb download -o taxadb
   $ taxadb create -i taxadb --fast --dbname taxadb.sqlite

This script assumes FASTA headers like this:

   >UniRef100_Q6GZW6 Putative helicase 009L n=2 Tax=Frog virus 3 TaxID=10493 RepID=009L_FRG3G

Specifically, the TaxID=? part.  That ID is queried to get a taxonomy list like this:

   ['Viruses', 'dsDNA viruses, no RNA stage', 'Iridoviridae', 'Ranavirus', 'Frog virus 3']
   or
   ['cellular organisms', 'Bacteria', 'Terrabacteria group', 'Firmicutes', 'Clostridia', 'Clostridiales', 'Clostridiaceae', 'Clostridium', 'Clostridium sp. ATCC 29733']

You can pass a comma-separated string of clades from any of these levels, such as:

   -c Viruses,Bacteria
   
The entire path with be searched for this, but only for exact matches (case-sensitive).
Note - this currently won't work with any clade strings which contain a comma.

"""

import argparse
import os
import re
from taxadb.taxid import TaxID

def main():
    parser = argparse.ArgumentParser( description='Filter a Uniref FASTA file by taxonomy')

    ## output file to be written
    parser.add_argument('-i', '--input_fasta', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_fasta', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-c', '--clades', type=str, required=True, help='Comma-separated string of clades to be included.' )
    parser.add_argument('-db', '--taxadb', type=str, required=True, help='Path to the taxadb sqlite3 file' )
    args = parser.parse_args()

    if not os.path.exists(args.taxadb):
        raise Exception("ERROR:  Couldn't find taxadb specified")
    
    taxid = TaxID(dbtype='sqlite', dbname=args.taxadb)
    clades = args.clades.split(',')

    record_count = 0
    print_every = 1000

    clade_counter = dict()
    for clade in clades:
        clade_counter[clade] = 0

    # remembers for each ID if we're keeping it or not
    id_cache = dict()

    fout = open(args.output_fasta, 'wt')
    keep_entry = False

    for line in open(args.input_fasta):
        if line[0] == '>':
            record_count += 1
            if record_count % print_every == 0:
                print("{0} records processed ...".format(record_count), flush=True)
            
            m = re.search('TaxID=(\d+)', line)
            if m:
                tax_id = m.group(1)

                if tax_id in id_cache:
                    if id_cache[tax_id] == True:
                        keep_entry = True
                    else:
                        keep_entry = False
                else:
                    lineage = taxid.lineage_name(tax_id, reverse=True)
                    clade_found = False

                    if lineage is None:
                        keep_entry = False
                        continue

                    for clade in clades:
                        if clade in lineage:
                            clade_found = True
                            clade_counter[clade] += 1
                            break

                    if clade_found:
                        keep_entry = True
                        id_cache[tax_id] = True
                    else:
                        keep_entry = False
                        id_cache[tax_id] = False
                        
            else:
                keep_entry = False

        if keep_entry:
            fout.write(line)
        
    fout.close()

    print("Number of entries exported by clade:")

    for clade in clade_counter:
        print("\t{0}: {1}".format(clade, clade_counter[clade]))

        
if __name__ == '__main__':
    main()







