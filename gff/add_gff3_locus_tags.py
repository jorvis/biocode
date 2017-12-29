#!/usr/bin/env python3

"""
This script is used to take any GFF3 file and add locus tags.  Be sure you've sorted your features
first, as this script assigns them in the order they are already found.

NCBI's specification on locus_tag usage:

   http://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation#locus_tag
   The locus_tag is a systematic gene identifier that is assigned to each gene. The locus_tag must be unique
   for every gene of a genome. Each genome project (i.e. all chromosomes and plasmids) have the same unique
   locus_tag prefix to ensure that a locus_tag is specific for a particular genome project, which is why we
   require that the locus_tag prefix be registered. In addition, a gene may have a biological name, as assigned
   in the scientific literature, and described above. For example, OBB_0001 is a systematic gene identifier, 
   while abcD is the biological gene name. We recommend having the BioProject registration process auto-assign 
   a locus_tag prefix, as they are not meant to convey meaning. The locus_tag prefix must be 3-12 alphanumeric 
   characters and the first character may not be a digit. Additionally locus_tag prefixes are case-sensitive. 
   The locus_tag prefix is followed by an underscore and then an alphanumeric identification number that is 
   unique within the given genome. Other than the single underscore used to separate the prefix from the 
   identification number, no other special characters can be used in the locus_tag. Locus_tags must only be 
   used in combination with a gene feature.

   RNA features (rRNA, tRNA, ncRNA) must include a corresponding gene feature with a locus_tag qualifier.

You must give this script an identifier prefix (--prefix).  By default the number of genes will be scanned and
the numeric portion of the locus tag will be zero-padded to one order of magnitude higher than the current
gene count.  So if --prefix=Tparva and your file contains 9000 genes, the locus_tags will be formatted like:

   Tparva_00001 .. Tparva_09000

By default, the numeric portion is assigned serially.  If you pass the --interval argument, there will be gaps
within the IDs, which allows for spacing for gene additions, splits, etc.  For example, if --interval=5
you'd have:

   Tparva_00005
   Tparva_00010
   Tparva_00015
   ...

You can pass an ID file if you have a list of identifiers that you know you already want applied.  It looks
in the first column for feature IDs (column 9 ID attribute) and in the second column for the locus tag to use.

Another option, perhaps an esoteric one, is to pass a molecule mapping file which will embed a token in each
locus tag containing a marker for that molecule.  Take the example Tparva_00005 ID above, if a molecule
map were passed with these contents:

    tp.assembly.567468735.1	01
    tp.assembly.567478335.1	02
    tp.assembly.567492885.1	03

The script will read the molecule name in column 1 of the GFF, match it to this table, and instead set
a locus tag of Tparva_01g00005
"""

import argparse
import re
import sys

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Adds locus tag identifiers to GFF3 features')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='TA file of source molecules' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Optional output file path (else STDOUT)' )
    parser.add_argument('-p', '--prefix', type=str, required=True, help='The prefix portion of IDs to be generated')
    parser.add_argument('-a', '--padding', type=int, required=True, help='Specify the minimum with to reserve for the numeric portion of the IDs.  Smaller numbers will be zero-padded.' )
    parser.add_argument('-n', '--interval', type=int, required=False, default=1, help='Interval between generated identifiers' )
    parser.add_argument('-s', '--starting_id', type=int, required=False, default=0, help='Initial numeric portion of IDs to be generated (do not zero-pad)' )
    parser.add_argument('-d', '--id_file', type=str, required=False, help='Pass a 2-column file of IDs to retain (in case you have mapped genes, for example)')
    parser.add_argument('-m', '--molecule_map', type=str, required=False, help='Pass a 2-column file of molecule->token identifiers (see documentation)')
    parser.add_argument('-c', '--custom', type=str, required=False, help='For custom parsing steps.  Most should ignore this.')

    args = parser.parse_args()
    check_arguments(args)

    # used to store locus_tags associated with each gene (so children can inherit)
    gene_loci = dict()
    next_id = args.starting_id
    last_molecule = None

    id_mapping  = parse_mapping_file( args.id_file )
    mol_mapping = parse_mapping_file( args.molecule_map )
    loci_assigned = list()

    ## if using Joana's custom options, check assumptions
    if args.custom == 'joana':
        if args.molecule_map is None or args.id_file is None:
            raise Exception("ERROR: Expected --molecule_map and --id_file options when using --custom=joana")
        else:
            ## need to process the ID map to reformat IDs
            for id in id_mapping:
                # TP05_0002 -> TpMuguga_05g00002
                m = re.match('TP(\d\d)_(\d+)', id_mapping[id])
                if m:
                    id_mapping[id] = "{0}_{1}g0{2}".format(args.prefix, m.group(1), m.group(2) )
                    
    elif args.custom == 'bmicroti':
        microti_map = { 'I':'01', 'II':'02', 'III':'03', 'IV':'04' }
        
        if  args.molecule_map is None or args.id_file is None:
            raise Exception("ERROR: Expected --molecule_map and --id_file options when using --custom=bmicroti")
        else:
            for id in id_mapping:
                m = re.match('BBM_(\D+)(\d+)', id_mapping[id])
                if m:
                    print("Changing id from {0} to ".format(id))
                    id_mapping[id] = "{0}_{1}g{2}".format(args.prefix, microti_map[m.group(1)], m.group(2) )
                    print(id_mapping[id])
                else:
                    raise Exception("ERROR: id ({0}) didn't match expected convention.".format(id_mapping[id]))
                    
        
    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output_file is not None:
        fout = open(args.output_file, 'wt')

    last_number_portion_assigned = 0

    for line in open(args.input_file):
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) != 9:
            fout.write(line + "\n")
            continue

        if last_molecule is None or (args.molecule_map is not None and mol_mapping[cols[0]] != mol_mapping[last_molecule]):
            print("Found molecule {0}, resetting id counter from {1}".format(cols[0], next_id) )
            next_id = args.starting_id
            last_molecule = cols[0]

        # grab the ID column if any
        id = gff.column_9_value(cols[8], 'ID')
        parent = gff.column_9_value(cols[8], 'Parent')
        type = cols[2]

        if type == 'gene':
            while True:
                if id in id_mapping:
                    locus_id = id_mapping[id]
                else:
                    if args.molecule_map is None:
                        locus_id = "{0}_{1}".format(args.prefix, str(next_id).zfill(args.padding))
                    else:
                        if cols[0] in mol_mapping:
                            if args.custom == 'bmicroti':
                                locus_id = "{0}_{2}g{1}".format(args.prefix, str(int(last_number_portion_assigned) + 1).zfill(args.padding), mol_mapping[cols[0]])
                            else:
                                locus_id = "{0}_{2}g{1}".format(args.prefix, str(next_id).zfill(args.padding), mol_mapping[cols[0]])
                        else:
                            raise Exception("ERROR: --molecule_map passed but {0} wasn't found in it.".format(cols[0]) )

                    next_id += args.interval

                cols[8] = gff.set_column_9_value(cols[8], 'locus_tag', locus_id)

                ## make sure this wasn't generated already (possibly conflict between --id_file and an
                #   auto-generated ID?
                if locus_id not in loci_assigned:
                    break
                else:
                    print("DEBUG: Duplicate ID assigned ({0}), trying again.".format(locus_id) )

            loci_assigned.append(locus_id)
            gene_loci[id] = locus_id

            m = re.search(r"(\d+)$", locus_id)
            if m:
                last_number_portion_assigned = m.group(1)
            
        elif type.endswith('RNA'):
            if parent in gene_loci:
                cols[8] = gff.set_column_9_value(cols[8], 'locus_tag', gene_loci[parent])
            else:
                raise Exception("ERROR: found RNA {0} whose parent {1} wasn't found yet".format(id, parent))
        
        fout.write("\t".join(cols) + "\n")



def parse_mapping_file( file ):
    id_map = dict()

    if file is None:
        return id_map

    for line in open(file):
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) == 2:
            id_map[cols[0]] = cols[1]
    
    return id_map



def check_arguments( args ):
    # Check the acceptable values for format
    pass




if __name__ == '__main__':
    main()







