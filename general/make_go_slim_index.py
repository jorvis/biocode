#!/usr/bin/env python3

"""
Use case example:
When we finish generating a genomic annotation we have a GFF3 file with functional
attributes, including GO terms.  We want to map these to a reduced set of GO terms
(GO slim) instead.  That DAG traversal can take a while, so this script is meant to
allow you to generate a key/value index GO slim.  

Input:
You provide a full go.obo file, and a file of selected terms.  Only the term id, name
and is_a relationship is really used.  The rest is ignored.

Output:
The output is a lookup with the highest-level keys being the different GO namespaces.  
Each of these values is a key-value lookup, where each key is one of the full set of GO terms
and the corresponding value is the closest slim term to which it maps.  Currently stored
as a python pickled dict() file or JSON text file, depending on option passed.

NOTES:
Terms with namespace: external are skipped.  If a match isn't found, None is the stored
value

"""

import argparse
import os
import pickle
import json
import re

import igraph


def main():
    parser = argparse.ArgumentParser( description='Creates a quick-lookup index of source terms to their slim counterparts')

    ## output file to be written
    parser.add_argument('-f', '--full_obo_file', type=str, required=True, help='Input file with ALL terms, such as go.obo' )
    parser.add_argument('-s', '--slim_obo_file', type=str, required=True, help='Input file with slim terms, such as go.slim.obo' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output file to be created.' )
    parser.add_argument('-of', '--output_format', type=str, required=False, default='pickle', help='Output format (pickle or json)' )
    args = parser.parse_args()

    # Parse the OBO file and store a term lookup as well as graphs for each namespace
    slim_terms, slim_graph = parse_obo_graph(args.slim_obo_file)
    full_terms, full_graph = parse_obo_graph(args.full_obo_file)

    # For each full term, walk up the tree until you hit one of the slim terms
    terms_processed = 0
    mapping = dict()
    for full_term_id in full_terms:
        full_term = full_terms[full_term_id]
        full_term_ns = full_term['ns']

        if full_term_ns not in mapping:
            mapping[full_term_ns] = dict()

        slim_id = check_parents(fg=full_graph, ft_ns=full_term_ns, term=full_term, depth=1,
                                 full_terms=full_terms, slim_terms=slim_terms)
        mapping[full_term_ns][full_term_id] = slim_id
        terms_processed += 1
        print("Processed {0} of {1} terms ... ".format(terms_processed, len(full_terms)), end='\r')

    if args.output_format == 'pickle':
        with open(args.output_file, 'wb') as f:
            pickle.dump(mapping, f, pickle.HIGHEST_PROTOCOL)
            
    elif args.output_format == 'json':
        with open(args.output_file, 'wt') as f:
            f.write(json.dumps(mapping))
    else:
        raise Exception("Error, --output_format must be either 'json' or 'pickle'")

    print("\nDone")

        
def check_parents(fg=None, ft_ns=None, term=None, depth=None, full_terms=None, slim_terms=None):
    """
    Does a recursive traversal up the GO tree until it finds a slim term
    """
    for v_idx in fg[ft_ns].neighbors(term['idx'], mode='out'):
        v_id = fg[ft_ns].vs[v_idx]['id']

        if v_id in slim_terms:
            return v_id
        else:
            return check_parents(fg=fg, ft_ns=ft_ns, term=full_terms[v_id], depth=depth + 1,
                                 full_terms=full_terms, slim_terms=slim_terms)

        return None
        

def parse_obo_graph(path):
    stored_pickles_found = False

    g = {'biological_process': igraph.Graph(directed=True), 
         'cellular_component': igraph.Graph(directed=True),
         'molecular_function': igraph.Graph(directed=True) }

    for ns in g:
        pickle_file_path = "{0}.{1}.graph".format(path, ns)
        if os.path.exists(pickle_file_path):
            print("Using stored ontology graph: {0}".format(pickle_file_path))
            g[ns] = igraph.Graph.Read_Pickle(fname=pickle_file_path)
            stored_pickles_found = True

    # key: GO:ID, value = {'ns': 'biological_process', 'idx': 25}
    terms = dict()

    if stored_pickles_found is True:
        print("Using stored terms data structure: {0}.terms".format(path))
        with open("{0}.terms".format(path), 'rb') as f:
            terms = pickle.load(f)

    # key: namespace, value=int
    next_idx = {'biological_process': 0, 
                'cellular_component': 0,
                'molecular_function': 0 }

    id = None
    namespace = None
    name = None

    # Pass through the file once just to get all the GO terms and their namespaces
    #  This makes the full pass far easier, since terms can be referenced which haven't
    #  been seen yet.

    if stored_pickles_found is False:
        for line in open(path):
            line = line.rstrip()
            if line.startswith('[Term]'):
                if id is not None:
                    # error checking
                    if namespace is None:
                        raise Exception("Didn't find a namespace for term {0}".format(id))

                    g[namespace].add_vertices(1)
                    idx = next_idx[namespace]
                    g[namespace].vs[idx]['id'] = id
                    g[namespace].vs[idx]['name'] = name
                    next_idx[namespace] += 1
                    terms[id] = {'ns': namespace, 'idx': idx}

                # reset for next term
                id = None
                namespace = None
                name = None

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
    
    id = None
    alt_ids = list()
    namespace = None
    name = None
    is_obsolete = False
    is_a = list()

    # Now actually parse the rest of the properties
    if stored_pickles_found is False:
        for line in open(path):
            line = line.rstrip()
            if line.startswith('[Term]'):
                if id is not None:
                    # make any edges in the graph
                    for is_a_id in is_a:
                        # these two terms should be in the same namespace
                        if terms[id]['ns'] != terms[is_a_id]['ns']:
                            raise Exception("is_a relationship found with terms in different namespaces")

                        #g[namespace].add_edges([(terms[id]['idx'], terms[is_a_id]['idx']), ])
                        # the line above is supposed to be able to instead be this, according to the 
                        # documentation, but it fails:
                        g[namespace].add_edge(terms[id]['idx'], terms[is_a_id]['idx'])

                # reset for this term
                id = None
                alt_ids = list()
                namespace = None
                is_obsolete = False
                is_a = list()

            elif line.startswith('id:'):
                id = line.split(' ')[1]
  
            elif line.startswith('namespace:'):
                namespace = line.split(' ')[1]

            elif line.startswith('is_a:'):
                is_a.append(line.split(' ')[1])

    if stored_pickles_found is False:
        for ns in g:
            pickle_file_path = "{0}.{1}.graph".format(path, ns)
            g[ns].write_pickle(fname=pickle_file_path)

        ## save the terms too so we don't have to redo that parse
        with open("{0}.terms".format(path), 'wb') as f:
            pickle.dump(terms, f, pickle.HIGHEST_PROTOCOL)

    return terms, g

                

if __name__ == '__main__':
    main()







