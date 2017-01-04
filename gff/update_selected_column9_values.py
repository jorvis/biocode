#!/usr/bin/env python3

"""

This script is useful when there is a systematic update you wish to make to any of the
key/value pairs in column 9 of any feature.  The use case for its creation was an annotation
GFF3 file in which we wanted to update the gene product names for a selection of the
annotated polypeptides.

The --update file should be tab-delimited with only two columns.  The first is the ID (though
it can be something else) of the feature that needs updating, and the second is the value.  The
attribute to update is passed via the --attribute flag.  So if I had a file representing gene
product names I wanted to update like this:

  0B1EFA9216B66647BA4D71E5C7357CB5    Fructose-bisphosphate aldolase 2
  57EA1FA54D81C35437D76C7758679AE2    Ribosomal RNA small subunit methyltransferase H
  3F038DEEA52E79D704E0EB9F1EDABF44    Conserved hypothetical protein

Then passed the --attribute='product_name' would determine which of the column 9 index keys
would have their values updated.

Usage example:

  ./update_selected_column9_values.py -i annotation.gff3 -u updates.tab -a 'product_name' -o annotation.v2.gff3

--------------------------------------------------------------------
EXAMPLE 2:  Updates by attributes other than ID (and type filtering)

GFF3 rows like this:

  ChromosomeI  IGS     gene    408117  410018  .       -       .       ID=A784FBCA42CADAA2CEB6E5CE6526346A;locus_tag=BBM_01g01120
  ChromosomeI  IGS     mRNA    408117  410018  .       -       .       ID=E46A8E9563679B07C55E1994A417038A;Parent=A784FBCA42CADAA2CEB6E5CE6526346A;locus_tag=BBM_01g01120
  ChromosomeI  IGS     CDS     409547  410018  .       -       0       ID=E46A8E9563679B07C55E1994A417038A-CDS;Parent=E46A8E9563679B07C55E1994A417038A
  ChromosomeI  IGS     CDS     408117  409525  .       -       0       ID=E46A8E9563679B07C55E1994A417038A-CDS;Parent=E46A8E9563679B07C55E1994A417038A
  ChromosomeI  IGS     exon    408117  409525  .       -       .       ID=ECBEE0925D9BA9694B8E491CFB4DAB55;Parent=E46A8E9563679B07C55E1994A417038A
  ChromosomeI  IGS     exon    409547  410018  .       -       .       ID=8144E805D7E6673862A61298D92F3246;Parent=E46A8E9563679B07C55E1994A417038A
  ChromosomeI  IGS     polypeptide     408117  410018  .       -       .       ID=E46A8E9563679B07C55E1994A417038A;Parent=E46A8E9563679B07C55E1994A417038A;Dbxref=overlaps_old_locusTagID:BBM_I01120;product_name=hypothetical protein

Rather than ID, the user has an update file keyed on locus_tag values.  The problem here is that they
aren't unique to the feature graph (in the example above, locus_tag=BBM_01g01120 on the gene, mRNA
and polypeptide feature types.)  So the --type=polypeptide argument should be used.

Second, to tell the script to key the update on locus_tag rather than ID, we set --key=locus_tag
Here's the full example:

  ./update_selected_column9_values.py -i annotation.gff3 -u updates.tab -t polypeptide -k 'locus_tag' -a 'product_name' -o annotation.v2.gff3


"""

import argparse

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Updates 9th-column key/value pairs in GFF file using a batch-update file')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='A GFF3 file' )
    parser.add_argument('-u', '--update_file', type=str, required=True, help='A two-column file (FeatureID, value)' )
    parser.add_argument('-a', '--attribute', type=str, required=True, help='The attribute value to update' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-k', '--key', type=str, required=False, default='ID', help='Which key in the 9th column helps identify the row to be updated?' )
    parser.add_argument('-t', '--type', type=str, required=False, help='Filter rows updated based on the 3rd (type) column' )
    args = parser.parse_args()

    # protect the user
    if args.input_file == args.output_file:
        raise Exception("ERROR:  Don't set --input_file and --output_file to be the same thing.  Bad things will happen.  Bad Things.")

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
            if args.type is None or args.type == cols[2]:
                atts = gff.column_9_dict(cols[8])
            
                if args.key in atts and atts[args.key] in changes:
                    atts[args.attribute] = changes[atts[args.key]]
                    cols[8] = gff.build_column_9_from_dict(atts)

        outfh.write("\t".join(cols) + "\n")


if __name__ == '__main__':
    main()







