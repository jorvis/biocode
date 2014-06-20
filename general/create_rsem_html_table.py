#!/usr/bin/env python3

"""
This script reads RSEM output and generates a single HTML page with a sortable,
searchable table of the output.

Author: Joshua Orvis (jorvis AT gmail)
"""

import argparse
import os

from jinja2 import Environment, FileSystemLoader

PATH = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_ENVIRONMENT = Environment(
    autoescape=False,
    loader=FileSystemLoader(os.path.join(PATH, '../data')),
    trim_blocks=False)

def main():
    parser = argparse.ArgumentParser( description='Creates an HTML page of RSEM output')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-t', '--title', type=str, required=False, default='', help='Puts a label at the top of the page')
    args = parser.parse_args()

    context = { 'transcripts':list(), 'title':args.title }

    for line in open(args.input_file):
        line = line.rstrip()
        if line.startswith("gene_id\ttranscript_id"):
            continue

        cols = line.split("\t")

        if len(cols) != 7:
            continue

        transcript = { 'gene_id':cols[0], 'id':cols[1], 'actual_length':cols[2],
                       'effective_length':cols[3], 'expected_count':cols[4],
                       'tpm':cols[5], 'fpkm':cols[6] }
        
        context['transcripts'].append(transcript)
        

    ofh = open(args.output_file, 'wt')
    tmpl = TEMPLATE_ENVIRONMENT.get_template('rsem_html_table.template').render(context)
    ofh.write(tmpl)
    

if __name__ == '__main__':
    main()







