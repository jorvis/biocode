#!/usr/bin/env python3

"""
Creates an index of OBO terms in either Python Pickle format or JSON

Input:
You provide a full go.obo file.  Only the term id, name and is_a relationship is 
really used.  The rest is ignored.

Output:
The output is a lookup with the highest-level keys being the different GO namespaces.  
Each of these values is a key-value lookup, where each key is one of the full set of GO terms
and the corresponding value is the closest slim term to which it maps.  Currently stored
as a python pickled dict() file or JSON text file, depending on option passed.

NOTES:
Terms with namespace: external are skipped.  

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
    parser.add_argument('-i', '--obo_file', type=str, required=True, help='Input file with ALL terms, such as go.obo' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output file to be created.' )
    parser.add_argument('-f', '--output_format', type=str, required=False, default='pickle', help='Output format (pickle or json)' )
    args = parser.parse_args()

    # Parse the OBO file and store a term lookup as well as graphs for each namespace
    full_terms, full_graph = parse_obo_graph(args.obo_file)

    if args.output_format == 'pickle':
        out_data = { 'terms': full_terms, 'graph': full_graph }
        with open(args.output_file, 'wb') as f:
            pickle.dump(out_data, f, pickle.HIGHEST_PROTOCOL)
            
    elif args.output_format == 'json':
        out_data = { 'terms': full_terms }        
        with open(args.output_file, 'wt') as f:
            f.write(json.dumps(out_data))
    else:
        raise Exception("Error, --output_format must be either 'json' or 'pickle'")

    print("\nDone")

        
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
        alt_ids = list()

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
                    terms[id] = {'ns': namespace, 'idx': idx, 'name': name}

                    # duplicate this for any aliases
                    for alt_id in alt_ids:
                        terms[alt_id] = {'ns': namespace, 'idx': idx, 'name': name}

                # reset for next term
                id = None
                namespace = None
                name = None
                alt_ids = list()

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

            elif line.startswith('alt_id: '):
                alt_ids.append(line.split(' ')[1])
    
    id = None
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

                        g[namespace].add_edge(terms[id]['idx'], terms[is_a_id]['idx'])

                # reset for this term
                id = None
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







