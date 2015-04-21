#!/usr/bin/env python3

"""

This script is useful when there is a systematic update you wish to make to any of the
key/value pairs in column 9 of any feature.  The use case for its creation was an annotation
GFF3 file in which we wanted to update the gene product names for a selection of the
annotated polypeptides.

The --update file should be tab-delimited with only two columns.  The first is the ID
of the feature that needs updating, and the second is the value.  The attribute to update is
passed via the --attribute flag.  So if I had a file representing gene product names I wanted
to update like this:

0B1EFA9216B66647BA4D71E5C7357CB5    Fructose-bisphosphate aldolase 2
57EA1FA54D81C35437D76C7758679AE2    Ribosomal RNA small subunit methyltransferase H
3F038DEEA52E79D704E0EB9F1EDABF44    Conserved hypothetical protein

Then passed the --attribute='product_name'would determine which of the column 9 index keys
would have their values updated.

Usage example:

./update_selected_column9_values.py -i annotation.gff3 -u updates.tab -a 'product_name' -o annotation.v2.gff3

"""

import argparse
import biocodegff
import os


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='A GFF3 file' )
    parser.add_argument('-u', '--update_file', type=str, required=True, help='A two-column file (FeatureID, value)' )
    parser.add_argument('-a', '--attribute', type=str, required=True, help='The attribute value to update' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    outfh = open(args.output_file, 'wt')

    # first read in the values to be updated
    changes = dict()

    for line in open(args.update_file):
        cols = line.rstrip().split("\t")
        if len(cols) != 2:
            print("WARNING: Skipping the following update line because two columns were expected:\n{0}".format(line))
            continue

        changes[cols[0]] = cols[1]

    for line in open(args.input_file):
        cols = line.rstrip().split("\t")

        if len(cols) == 9:
            atts = biocodegff.column_9_dict(cols[8])
            
            if 'ID' in atts and atts['ID'] in changes:
                atts[args.attribute] = changes[atts['ID']]
                cols[8] = biocodegff.build_column_9_from_dict(atts)

        outfh.write("\t".join(cols) + "\n")


if __name__ == '__main__':
    main()







