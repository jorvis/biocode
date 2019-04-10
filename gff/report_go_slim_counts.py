#!/usr/bin/env python3

"""
This script uses the following input:

- Annotated GFF3 file
- Full go.obo file
- Pickled go slim

And produces count output in plain text like this:

cellular_component	GO:0005575	cellular_component	650
cellular_component	GO:0016020	membrane	135
cellular_component	GO:0043226	organelle	51
cellular_component	GO:0032991	macromolecular complex	67
cellular_component	GO:0030054	cell junction	1
biological_process	GO:0006810	transport	195
biological_process	GO:0008150	biological_process	50
biological_process	GO:0008152	metabolic process	864
biological_process	GO:0050896	response to stimulus	112
biological_process	GO:0065007	biological regulation	118
biological_process	GO:0009987	cellular process	98
...


"""

import argparse
from biocode import utils, gff
import json
import os
import pickle
import re


def main():
    parser = argparse.ArgumentParser( description='Reporter of GO slim counts from annotated GFF3 file')
    parser.add_argument('-i', '--input_gff', type=str, required=True, help='Path to an input GFF3 file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-os', '--obo_slim_pickle', type=str, required=True, help='Pickled OBO slim file' )
    parser.add_argument('-fo', '--full_obo', type=str, required=True, help='Full OBO slim file' )
    args = parser.parse_args()

    source_go_terms = parse_go_terms_from_gff(args.input_gff)
    full_obo_index = parse_full_go(args.full_obo)
    slim_counts = map_to_slim(source_terms=source_go_terms, slim_map_file=args.obo_slim_pickle)

    out_fh = open(args.output_file, 'wt')
    
    for ns in slim_counts:
        if ns == 'unknown': continue
        for term in slim_counts[ns]:
            if term == 'unknown': continue
            print("{0}\t{1}\t{2}\t{3}".format(ns, term, full_obo_index[term]['name'], slim_counts[ns][term]), file=out_fh)

    out_fh.close()


# Goal is to create summary slim counts on main page like this:
#   http://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0130720
# Great viewer
#   http://visualdataweb.de/webvowl
def map_to_slim(source_terms=None, slim_map_file=None):
    # Slim terms to skip.  Usually the root Big Three terms
    skip_slim_terms = ['GO:0008150', 'GO:0003674', 'GO:0005575']
    slim_map = pickle.load( open(slim_map_file, 'rb') )

    counts = dict()

    for ns in slim_map:
        counts[ns] = {'unknown': 0}
    
    for gff_term in source_terms:
        gff_id = "GO:{0}".format(gff_term)
        slim_term = None
        
        for ns in slim_map:
            if gff_id in slim_map[ns]:
                slim_term = slim_map[ns][gff_id]
                if slim_term is None:
                    counts[ns]['unknown'] += source_terms[gff_term]
                else:
                    if slim_term not in counts[ns]:
                        counts[ns][slim_term] = 0

                    counts[ns][slim_term] += source_terms[gff_term]

                break

    return counts

def parse_full_go(file):
    id = None
    namespace = None
    name = None
    terms = dict()
    
    for line in open(file):
        line = line.rstrip()
        
        if line.startswith('[Term]'):
             if id is not None:
                 terms[id] = {'ns': namespace, 'name':name}

        elif line.startswith('id:'):
                id = line.split(' ')[1]

        elif line.startswith('namespace:'):
            namespace = line.split(' ')[1]

        elif line.startswith('name:'):
            m = re.match('name: (.+)', line)
            if m:
                name = m.group(1).rstrip()
            else:
                raise Exception("Failed to regex this line: {0}".format(line))

    return terms

def parse_go_terms_from_gff(file):
    terms = dict()
    assemblies, features = gff.get_gff3_features(file)
    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            for mRNA in gene.mRNAs():
                for polypeptide in mRNA.polypeptides():
                    annot = polypeptide.annotation
                    for go_annot in annot.go_annotations:
                        if go_annot.go_id in terms:
                            terms[go_annot.go_id] += 1
                        else:
                            terms[go_annot.go_id] = 1

    return terms

if __name__ == '__main__':
    main()







