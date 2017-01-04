#!/usr/bin/env python3

"""
The input expected is a standards-compliant GFF3 file.  

  http://www.sequenceontology.org/gff3.shtml

Example of a gene with GO annotation:

  AAGK01000005	.	gene	161052	164192	.	+	.	ID=0ACE4B1023C0E198886936CD6BE3869B;locus_tag=TpMuguga_03g00084
  AAGK01000005	.	mRNA	161052	164192	.	+	.	ID=99A7DA3C31446D989D46A5B0C8EC5C0E;Parent=0ACE4B1023C0E198886936CD6BE3869B;locus_tag=TpMuguga_03g00084
  AAGK01000005	.	CDS	161337	162111	.	+	0	ID=470AEE91631AFCE9DBFD1CF9BA0E7365;Parent=99A7DA3C31446D989D46A5B0C8EC5C0E
  AAGK01000005	.	CDS	162179	163696	.	+	0	ID=470AEE91631AFCE9DBFD1CF9BA0E7365;Parent=99A7DA3C31446D989D46A5B0C8EC5C0E
  AAGK01000005	.	CDS	163868	164178	.	+	0	ID=470AEE91631AFCE9DBFD1CF9BA0E7365;Parent=99A7DA3C31446D989D46A5B0C8EC5C0E
  AAGK01000005	.	exon	163868	164192	.	+	.	ID=B1A06E24F7020FC4E92B7F2F1099D059;Parent=99A7DA3C31446D989D46A5B0C8EC5C0E
  AAGK01000005	.	exon	162179	163696	.	+	.	ID=BBBED213F1613EEBE471FA356BC818E4;Parent=99A7DA3C31446D989D46A5B0C8EC5C0E
  AAGK01000005	.	exon	161052	162111	.	+	.	ID=BD61FC870CA44840B1817F831519608A;Parent=99A7DA3C31446D989D46A5B0C8EC5C0E
  AAGK01000005	.	polypeptide	161052	164192	.	+	.	ID=99A7DA3C31446D989D46A5B0C8EC5C0E;Parent=99A7DA3C31446D989D46A5B0C8EC5C0E;gene_symbol=gcs-1;product_name=Glutamate-cysteine ligase;Ontology_term=GO:0004357,GO:0006750;Dbxref=EC:6.3.2.2

Most of the other options correspond to required fields in the GAF spec.  Here are
those arguments and the columns they correspond to:

  -db : 1
  
  
OUTPUT: 

GO Gene Association File 2.0 (GAF) format is specified here:

  http://geneontology.org/page/go-annotation-file-formats

This output is also appropriate for use in slimming tools such as map2slim:

  http://search.cpan.org/~cmungall/go-perl/scripts/map2slim

LIMITATIONS:

- I do not currently have support for representation of evidence_code / with_from in the GFF3, so
  it can't be supported here.  Instead, the default evidence code column (7) will just be IEA
  
"""

import argparse
import os
import re
import sys
import time

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Converts GFF3 files to GO Gene Association Format (GAF)')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-go', '--go_file', type=str, required=True, help='Gene Ontology (GO) file' )
    parser.add_argument('-db', '--database', type=str, required=True, help='Database issuing that IDs.  Example: UniProtKB' )
    parser.add_argument('-dbref', '--db_reference', type=str, required=True, help='DB reference, like PMID:2676709 (column 6)' )
    parser.add_argument('-ec', '--evidence_code', type=str, required=False, default='IEA', help='Like IEA (column 7)' )
    parser.add_argument('-t', '--taxon_id', type=int, required=True, help='NCBI taxon ID (column 13)' )
    parser.add_argument('-ad', '--annotation_date', type=str, required=False, help='Annotation date in YYYYMMDD format.  Default = GFF3 file datestamp' )
    parser.add_argument('-ab', '--assign_by', type=str, required=False, help='Assign by (column 15)  Defaults to --database argument value' )
    args = parser.parse_args()

    print("INFO: Parsing GFF3 objects", file=sys.stderr)
    (assemblies, features) = gff.get_gff3_features(args.input_file)

    print("INFO: Parsing GO file", file=sys.stderr)
    go_lookup = parse_go_file(args.go_file)

    annot_date = args.annotation_date
    if annot_date is None:
        annot_date = time.strftime('%Y%m%d', time.gmtime(os.path.getmtime(args.input_file)))

    assign_by = args.assign_by
    if assign_by is None:
        assign_by = args.database

    ofh = open(args.output_file, 'wt')
    
    ofh.write("!gaf-version: 2.0\n")
     
    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            for mRNA in gene.mRNAs():
                for polypeptide in mRNA.polypeptides():
                    for go_annot in polypeptide.annotation.go_annotations:
                        go_id = "GO:{0}".format(go_annot.go_id)
                        product = None
                        gene_sym = None
                        
                        if go_id not in go_lookup:
                            raise Exception("ERROR: GO ID {0} not found in provided go.obo file".format(go_id))

                        if polypeptide.annotation.product_name is not None: product = polypeptide.annotation.product_name
                        if polypeptide.annotation.gene_symbol is not None:  gene_sym = polypeptide.annotation.gene_symbol
                        
                        
                        # Aspect is F, P or C, depending on which component/ontology the term comes from
                        ofh.write("{0}\t{1}\t{1}\t\t{2}\t{3}\t{4}\t\t{5}\t{6}"
                                  "\t{7}\tprotein\ttaxon:{8}\t{9}\t{10}\t"
                                  "\t\n".format(args.database, polypeptide.id, go_id, args.db_reference,
                                                args.evidence_code, go_lookup[go_id], product, gene_sym,
                                                args.taxon_id, annot_date, assign_by))

    print("INFO: Conversion complete.", file=sys.stderr)

def parse_go_file(path):
    lookup = dict()
    primary_id = None
    alt_ids = list()
    namespace = None

    for line in open(path):
        line = line.rstrip()

        # if we found a new term, store the data from the last, then reset
        if line.startswith('[Term]'):
            if primary_id is not None:
                lookup[primary_id] = namespace
                
                for alt_id in alt_ids:
                    lookup[alt_id] = namespace

            primary_id = None
            alt_ids = list()
            namespace = None

        m = re.search('^id: (\S+)', line)
        if m:
            primary_id = m.group(1)
        else:
            m = re.search('^namespace: (\S+)', line)
            if m:
                if m.group(1) == 'biological_process':
                    namespace = 'P'
                elif m.group(1) == 'molecular_function':
                    namespace = 'F'
                elif m.group(1) == 'cellular_component':
                    namespace = 'C'

            else:
                m = re.search('^alt_id: (\S+)', line)
                if m:
                    alt_ids.append(m.group(1))

    print("INFO: \t{0} terms parsed".format(len(lookup)))
    return lookup
                    
if __name__ == '__main__':
    main()







