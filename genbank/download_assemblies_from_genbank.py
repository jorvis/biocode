#!/usr/bin/env python3

"""
This script is used to download the nucleotide assemblies (contigs/scaffolds) of one or more 
genome projects based on the parent Genbank ID of each.  It uses EUtils for this, such as:

  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=ACIN00000000&retmode=xml&rettype=gb

If you pass the --api_key option (required after May 1, 2018), the 1-second delay between accession 
list entries is skipped, since NCBI allows faster API calls if you have a key.  To get one, register:

  https://www.ncbi.nlm.nih.gov/books/NBK25497/#_chapter2_Usage_Guidelines_and_Requiremen_

How each ID is handled depends on the ID format.  Those like this are master records for assemblies:

   ACKQ00000000
   ACDJ00000000
   LSCU00000000
   ABWL00000000

While shorter ones are direct assembled contigs:

   CP002390
   CP007064
   KR131712

For each ID passed, this script first gets the master record XML file for it and looks for this section:

    <GBSeq_alt-seq>
      <GBAltSeqData>
        <GBAltSeqData_name>WGS</GBAltSeqData_name>
        <GBAltSeqData_items>
          <GBAltSeqItem>
          <GBAltSeqItem_first-accn>ACIN03000001</GBAltSeqItem_first-accn>
          <GBAltSeqItem_last-accn>ACIN03000020</GBAltSeqItem_last-accn>
          </GBAltSeqItem>
        </GBAltSeqData_items>
      </GBAltSeqData>
      <GBAltSeqData>
        <GBAltSeqData_name>WGS_SCAFLD</GBAltSeqData_name>
        <GBAltSeqData_items>
          <GBAltSeqItem>
          <GBAltSeqItem_first-accn>KI535340</GBAltSeqItem_first-accn>
          <GBAltSeqItem_last-accn>KI535343</GBAltSeqItem_last-accn>
          </GBAltSeqItem>
        </GBAltSeqData_items>
      </GBAltSeqData>
    </GBSeq_alt-seq>

Separate FASTA files are created for the WGS and WGS_SCAFLD ranges (if present).  For
each range, ESearch is used to get an id range:

  <eSearchResult>
    <Count>4</Count>
    <RetMax>4</RetMax>
    <RetStart>0</RetStart>
    <IdList>
      <Id>555661647</Id>
      <Id>555661645</Id>
      <Id>555661643</Id>
      <Id>555661641</Id>
    </IdList>
  ...

Then EFetch is used again to pull each of these IDs in FASTA format.  Whatever output
directory you pass will be used as the base, and a subdirectory will be created with the
accession name.  So, if I pass '-a ACIN00000000 -o ./' I'll get this output tree:

.
└── ACIN00000000
    ├── accession_ranges.xml
    ├── master.xml
    ├── wgs_accessions.fasta
    └── wgs_scaffold_accessions.fasta

"""

import argparse
import os
import requests
import sys
import time
import xml.etree.ElementTree as ET

def main():
    parser = argparse.ArgumentParser( description='Download one or more genomic assemblies from Genbank')
    parser.add_argument('-i', '--id_list', type=str, required=False, help='List of Genbank accessions in a file, one per line')
    parser.add_argument('-a', '--accession', type=str, required=False, help='Single Genbank accession')
    parser.add_argument('-o', '--output_directory', type=str, required=True, help='Base directory where output files will be written' )
    parser.add_argument('--skip_existing', dest='skip_existing', action='store_true', help="Skip any entries for which the output directory already exists.  Useful if restarting after a network error.")
    parser.add_argument('-ak', '--api_key', type=str, required=True, help='NCBI requires an API key for all E-utility calls as of May 1, 2018.  https://www.ncbi.nlm.nih.gov/books/NBK25500/#_chapter1_Announcement_')
    args = parser.parse_args()

    if args.id_list is None and args.accession is None:
        raise Exception("You must pass either --id_list or --accession")

    if not os.path.exists(args.output_directory):
        raise Exception("The passed --output_directory doesn't exist")
    
    accessions = get_accession_list_from_args(args)
    eutils_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={0}&retmode=xml&rettype=gb&tool=biocode'

    if args.api_key is not None:
        eutils_base += "&api_key={0}".format(args.api_key)

    for acc in accessions:
        url = eutils_base.format(acc)

        print("Processing accession: {0}".format(acc))
        r = requests.get(url)

        if r.status_code == 200:
            outdir = "{0}/{1}".format(args.output_directory, acc)

            if os.path.exists(outdir):
                if args.skip_existing:
                    print("\tSkipped because --skip_existing and directory exists")
                    continue
            else:
                os.mkdir(outdir)

            # ids like 'ACKQ00000000' are processed as master records
            if len(acc) == 12:
                ranges = process_master_xml(acc=acc, text=r.text, output_dir=args.output_directory)

                if ranges['wgs'] is not None:
                    ids = get_accessions_from_id_range(accs=ranges['wgs'], output_dir=outdir)
                    make_multi_fasta(ids=ids, output_dir=outdir, label='wgs_accessions')
                    
                if ranges['wgs_scaffold'] is not None:
                    ids = get_accessions_from_id_range(accs=ranges['wgs_scaffold'], output_dir=outdir)
                    make_multi_fasta(ids=ids, output_dir=outdir, label='wgs_scaffold_accessions')
            else:
                make_multi_fasta(ids=[acc], output_dir=outdir)
                
        else:
            print("\t{0}\t{1}\tERROR".format(acc, r.status_code))
            continue

        if args.api_key is None:
            time.sleep(1)


def get_accession_list_from_args(args):
    accessions = list()

    if args.id_list is None and args.accession is None:
        raise Exception("You must pass either --id_list or --accession")

    if args.accession is not None:
        accessions.append(args.accession)

    if args.id_list is not None:
        for line in open(args.id_list):
            line = line.rstrip()
            accessions.append(line)

    return accessions

def get_accessions_from_id_range(accs=None, output_dir=None):
    esearch_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore"

    # this gets me an XML ID range:
    #  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=KI535340:KI535343[accn]
    if accs[1] is None:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={0}[accn]&tool=biocode&retmax=100000".format(accs[0])
    else:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={0}:{1}[accn]&tool=biocode&retmax=100000".format(accs[0], accs[1])

    r = requests.get(url)
    acc_xml_path = "{0}/accession_ranges.xml".format(output_dir)

    if r.status_code == 200:
        ids_fh = open(acc_xml_path, 'wt')
        ids_fh.write(r.text)
        ids_fh.close()
    else:
        print("\t{0}\t{1}\tERROR".format(accs[0], accs[1]))

    tree = ET.parse(acc_xml_path)
    root = tree.getroot()

    # Make sure the record isn't too large.  This would require multiple fetches to get the entire thing, which isn't currently supported
    record_count = int(root.find('./Count').text)
    retmax = int(root.find('./RetMax').text)
    
    if record_count > retmax:
        raise Exception("ERROR: The number of records for this accession/range ({0}) is greater than what NCBI will return at once ({1}).".format(record_count, retmax))
        
    ids = list()

    for id_elm in root.find('./IdList'):
        ids.append(id_elm.text)

    return ids


def make_multi_fasta(ids=None, output_dir=None, label=None):
    url_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={0}&rettype=fasta&retmode=text&tool=biocode"
    ofh = open("{0}/{1}.fasta".format(output_dir, label), 'wt')

    for id in ids:
        url = url_base.format(id)
        #print("Getting URL: {0}".format(url))
        r = requests.get(url)
        
        if r.status_code == 200:
            ofh.write(r.text.rstrip())
            ofh.write("\n")
        else:
            print("\t{0}\t{1}\tERROR getting individual fasta".format(id, r.status_code))
            continue

    ofh.close()

def process_master_xml(acc=None, text=None, output_dir=None):
    acc_dir = "{0}/{1}".format(output_dir, acc)
    master_xml_path = "{0}/master.xml".format(acc_dir)

    if not os.path.exists(acc_dir):
        os.mkdir(acc_dir)

    master_fh = open(master_xml_path, 'wt')
    master_fh.write(text)
    master_fh.close()

    tree = ET.parse(master_xml_path)
    root = tree.getroot()

    wgs_range = None
    wgs_scaffold_range = None
    
    for item in root.findall('./GBSeq/GBSeq_alt-seq/GBAltSeqData'):
        name_elm = item.find('GBAltSeqData_name')
        name = name_elm.text

        GBAltSeqItem = item.find('GBAltSeqData_items')[0]
        first_acc = GBAltSeqItem.find('GBAltSeqItem_first-accn').text

        last_acc_elm = GBAltSeqItem.find('GBAltSeqItem_last-accn')

        if last_acc_elm is None:
            last_acc = None
        else:
            last_acc  = GBAltSeqItem.find('GBAltSeqItem_last-accn').text

        if name == 'WGS':
            wgs_range = [first_acc, last_acc]
        elif name == 'WGS_SCAFLD':
            wgs_scaffold_range = [first_acc, last_acc]

    return {'wgs': wgs_range, 'wgs_scaffold': wgs_scaffold_range}


if __name__ == '__main__':
    main()







