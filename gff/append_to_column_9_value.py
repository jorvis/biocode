#!/usr/bin/env python3

"""

You might encounter GFF files with different attributes in the column 9 key=value
pairs, and want to systematically append a string onto each of them.  

Example use case:

Someone duplicated mRNA rows to create polypeptide rows in GFF, but failed to 
update the IDs so that they would be globally unique.  This script can be used to
append '-p' to each of the IDs so that they are like this:

./append_to_column_9_value.py -i input.gff3 -o output.gff3 -t polypeptide -k ID -a '_p'

This would change values like this:

gnl|WGS:AAGK|SeqID|gb|AAGK01000002  .  polypeptide  213984  215011  .  +  .  ID=5C1BDD5A6CB64E2F5259BC6F7505253A;

To this:

gnl|WGS:AAGK|SeqID|gb|AAGK01000002  .  polypeptide  213984  215011  .  +  .  ID=5C1BDD5A6CB64E2F5259BC6F7505253A_p;

NOTE:  The argparser fails if you attempt to pass a value which begins with a - symbol, 
even if in quotes.  This is why I did '_p' rather than '-p' in the example.

If you instead wanted to systematically update values in column 9 based on a series
of different inputs, you should use the biocode gff/update_selected_column9_values.py
script instead.

"""

import argparse

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Append a string to the end of existing key values in column 9 of GFF3 files')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-t', '--feature_type', type=str, required=False, help='Only update features of a given type (gff column 3)' )
    parser.add_argument('-k', '--key', type=str, required=True, help='Which column 9 key/attribute to update values for?' )
    parser.add_argument('-a', '--appended_text', type=str, required=True, help='Text to append' )
    args = parser.parse_args()

    # protect the user
    if args.input_file == args.output_file:
        raise Exception("ERROR:  Don't set --input_file and --output_file to be the same thing.  Bad things will happen.  Bad Things.")

    ofh = open(args.output_file, 'wt')
    replacement_count = 0

    for line in open(args.input_file):
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) == 9:
            if args.feature_type is None or args.feature_type == cols[2]:
                col9 = gff.column_9_dict(cols[8])
                if args.key in col9:
                    col9[args.key] = "{0}{1}".format(col9[args.key], args.appended_text)
                    cols[8] = gff.build_column_9_from_dict(col9)
                    replacement_count += 1

                ofh.write("{0}\t{1}\n".format("\t".join(cols[0:8]), cols[8]))
            else:
                ofh.write("{0}\n".format(line))
        else:
            ofh.write("{0}\n".format(line))

    print("INFO: Made {0} replacements in the file".format(replacement_count))

if __name__ == '__main__':
    main()







