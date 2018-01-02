#!/usr/bin/env python3

"""

Notes:  I wanted to use this C-based graph library, but it's apt source
URLs were missing files AND the site used to log issues was failing
Google App Auth.  Not a good sign?
  https://graph-tool.skewed.de

igraph notes:

NOTES:
Terms with namespace: external are skipped.

Developer:
Once parsed, the graphs can be saved with Graph.write_pickle() : http://igraph.org/python/doc/igraph.Graph-class.html#write_pickle

"""

import argparse
import os
import pickle
import re

import igraph
from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Creates a slim version of a set of ontology terms.')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Input file with one term per line. This is your file with ALL terms.' )
    parser.add_argument('-s', '--slim_terms', type=str, required=True, help='Plain text file with one term per line - your slim.' )
    parser.add_argument('-obo', '--ontology_file', type=str, required=True, help='Full obo file, providing the network of terms' )
    args = parser.parse_args()

    # Parse the OBO file and store a term lookup as well as graphs for each namespace
    terms, g = parse_obo_graph(args.ontology_file)

    # parse list of source GO terms
    source_go_terms = dict()
    assemblies, features = gff.get_gff3_features(args.input_file)
    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            for mRNA in gene.mRNAs():
                for polypeptide in mRNA.polypeptides():
                    annot = polypeptide.annotation
                    for go_annot in annot.go_annotations:
                        if go_annot.go_id in source_go_terms:
                            source_go_terms[go_annot.go_id] += 1
                        else:
                            source_go_terms[go_annot.go_id] = 1

    print("Slimming {0} unique source GO terms".format(len(source_go_terms)))

    # For testing, show all descendents of
    # GO:0008150 - biological_process
    # GO:0005575 - cellular_component
    # GO:0003674 - molecular_function
    biological_process_idx = terms['GO:0008150']['idx']
    slim_targets = {'unknown': 0}
    for edge_id in g['biological_process'].incident(biological_process_idx, mode='IN'):
        #print("Found edge id {0}".format(edge_id))
        # only looking at those getting more specific
        source_idx = g['biological_process'].es[edge_id].source
        target_idx = g['biological_process'].es[edge_id].target
        slim_targets[source_idx] = 0
  
    # how many do we find of each target?
    for source_id in source_go_terms:
        source_id = "GO:{0}".format(source_id)
        matching_targets = list()
        best_path_dist = 1000

        # check this annotated GO ID against the SLIM targets
        if source_id in terms:
            for target_idx in slim_targets:
                if target_idx != 'unknown':
                    paths = g['biological_process'].shortest_paths_dijkstra(source=terms[source_id]['idx'], target=target_idx)
                    path_dist = paths[0][0]
                    if path_dist != float('inf'):
                        if path_dist < best_path_dist:
                            best_path_dist = path_dist
                            matching_targets = [target_idx]
                        elif path_dist == best_path_dist:
                            matching_targets.append(target_idx)

        if len(matching_targets) > 0:
            for t_idx in matching_targets:
                slim_targets[t_idx] += 1
        else:
            slim_targets['unknown'] += 1

    print("Slim counts:")
    for id in slim_targets:
        print("\t{0} - {1}".format(id, slim_targets[id]))


def parse_obo_graph(path):
    stored_pickle_file_prefix = 'obo.graphs'
    stored_pickles_found = False

    g = {'biological_process': igraph.Graph(directed=True), 
         'cellular_component': igraph.Graph(directed=True),
         'molecular_function': igraph.Graph(directed=True) }

    for ns in g:
        pickle_file_path = "{0}.{1}".format(stored_pickle_file_prefix, ns)
        if os.path.exists(pickle_file_path):
            print("Using stored ontology graph: {0}".format(pickle_file_path))
            g[ns] = igraph.Graph.Read_Pickle(fname=pickle_file_path)
            stored_pickles_found = True

    # key: GO:ID, value = {'ns': 'biological_process', 'idx': 25}
    terms = dict()

    if stored_pickles_found is True:
        print("Using stored terms data structure: {0}".format(pickle_file_path))
        with open("{0}.terms".format(stored_pickle_file_prefix), 'rb') as f:
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
            pickle_file_path = "{0}.{1}".format(stored_pickle_file_prefix, ns)
            g[ns].write_pickle(fname=pickle_file_path)

        ## save the terms too so we don't have to redo that parse
        with open("{0}.terms".format(stored_pickle_file_prefix), 'wb') as f:
            pickle.dump(terms, f, pickle.HIGHEST_PROTOCOL)

    return terms, g

                

if __name__ == '__main__':
    main()







