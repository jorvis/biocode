#!/usr/bin/env python3

"""
The annotation pipeline is pretty liberal about the input molecule identifiers.  Once
we've reached the EVM step we assign IGS-convention identifiers.  This script will
export the new GFF3 file as well as a plain text file mapping the old and new IDs.


INPUT
-----
jcf7180000271712        EVM     gene    407     1013    .       +       .       ID=evm.TU.jcf7180000271712.1;Name=EVM%20prediction%20jcf7180000271712.1
jcf7180000271712        EVM     mRNA    407     1013    .       +       .       ID=evm.model.jcf7180000271712.1;Parent=evm.TU.jcf7180000271712.1
jcf7180000271712        EVM     exon    407     466     .       +       .       ID=evm.model.jcf7180000271712.1.exon1;Parent=evm.model.jcf7180000271712.1
jcf7180000271712        EVM     CDS     407     466     .       +       0       ID=cds.evm.model.jcf7180000271712.1;Parent=evm.model.jcf7180000271712.1
jcf7180000271712        EVM     exon    531     729     .       +       .       ID=evm.model.jcf7180000271712.1.exon2;Parent=evm.model.jcf7180000271712.1
jcf7180000271712        EVM     CDS     531     729     .       +       0       ID=cds.evm.model.jcf7180000271712.1;Parent=evm.model.jcf7180000271712.1
jcf7180000271712        EVM     exon    799     1013    .       +       .       ID=evm.model.jcf7180000271712.1.exon3;Parent=evm.model.jcf7180000271712.1
jcf7180000271712        EVM     CDS     799     1013    .       +       2       ID=cds.evm.model.jcf7180000271712.1;Parent=evm.model.jcf7180000271712.1

ID naming convention expected from EVM (where {asm} is the input assembly ID):

gene: evm.TU.{asm}.1
   mRNA: evm.model.{asm}.1
      exon: evm.model.{asm}.1.exon1
      exon: evm.model.{asm}.1.exon2
      CDS : cds.evm.model.{asm}.1  (CDS properly share IDs across fragments)

      
"""

import argparse

from biocode import gff

next_id = { 'gene':1, 'mRNA':1, 'exon':1, 'CDS':1 }

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_gff', type=str, required=True, help='Path to an output GFF file to be created with new IDs' )
    parser.add_argument('-p', '--id_prefix', type=str, required=True, help='Will be used as the base for all IDs generated' )
    parser.add_argument('-m', '--output_map', type=str, required=False, help='This will create a tab-delimited mapping of old and new IDs' )
    args = parser.parse_args()

    ofh = open(args.output_gff, 'w')

    if args.output_map is None:
        map_ofh = None
    else:
        map_ofh = open(args.output_map, 'w')

    idmap = dict()

    for line in open(args.input_file):
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) != 9:
            ofh.write(line + "\n")
            continue

        feat_id   = gff.column_9_value(cols[8], 'ID')
        parent_id = gff.column_9_value(cols[8], 'Parent')
        
        if feat_id in idmap:
            new_feat_id = idmap[feat_id]
        else:
            new_feat_id = get_new_id(args.id_prefix, cols[2], feat_id, map_ofh)
            idmap[feat_id] = new_feat_id

        if parent_id is None:
            cols[8] = "ID={0}".format(new_feat_id)
        else:
            if parent_id in idmap:
                new_parent_id = idmap[parent_id]
            else:
                new_parent_id = get_new_id(args.id_prefix, cols[2], parent_id, map_ofh)
                idmap[parent_id] = new_parent_id

            cols[8] = "ID={0};Parent={1}".format(new_feat_id, new_parent_id)

        ofh.write( "\t".join(cols) + "\n" )


def get_new_id(prefix, feat_type, old_id, ofh):
    new_id = "{0}.{1}.{2}.1".format(prefix, feat_type, next_id[feat_type])
    next_id[feat_type] += 1

    if ofh is not None:
        ofh.write("{0}\t{1}\n".format(old_id, new_id))
        
    return new_id

            
if __name__ == '__main__':
    main()







