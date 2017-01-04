#!/usr/bin/env python3

"""
This script is used to take any GFF3 file and re-generate feature identifiers within
it to match the convention used at IGS.  This is:

   $prefix.$type.$id.$version

The mode here defines what the identifier.  For example, if using --mode=sequential for 
an organism (--prefix) of b_microti, the ID for the 100th gene might be:

The mode here defines what the identifier.  Examples:

    --mode=sequential  (the default)
         b_microti.gene.100.1

    --mode=uuid    (UUID, as specified in RFC 4122, using uuid4())
         b_microti.gene.8d8f9231-262e-48e7-b066-a84b6a939746.1

    --mode=hex8
         b_microti.gene.c08ca446.1

    --mode=hex12
         b_microti.gene.191ccac20a56.1

The only values that are replaced are in the ID and Parent attributes in the 9th column.

"""

import argparse
import os
import sys
from binascii import hexlify
from collections import defaultdict
from uuid import uuid4

from biocode import gff

## constants
next_ids_sequential = defaultdict(lambda: 1)

def main():
    parser = argparse.ArgumentParser( description='Generates new identifiers in GFF3 files following the IGS identifier convention.')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='TA file of source molecules' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Optional output file path (else STDOUT)' )
    parser.add_argument('-p', '--prefix', type=str, required=True, help='The prefix portion of IDs to be generated')
    parser.add_argument('-m', '--mode', type=str, required=False, default='sequential', help='ID modes (see embedded documentation): sequential, uuid, hex8, hex12')

    args = parser.parse_args()
    check_arguments(args)

    id_map = dict()
    
    ## output will either be a file or STDOUT
    fout = sys.stdout
    if args.output_file is not None:
        fout = open(args.output_file, 'wt')

    for line in open(args.input_file):
        line = line.rstrip()

        cols = line.split("\t")

        if len(cols) != 9:
            fout.write(line + "\n")
            continue

        # grab the ID column if any
        id = gff.column_9_value(cols[8], 'ID')
        parent = gff.column_9_value(cols[8], 'Parent')
        new_id = None
        new_parent = None
        type = cols[2]

        if id is not None:
            if id in id_map:
                new_id = id_map[id]
            else:
                new_id = get_new_id(args.prefix, type, args.mode)
                id_map[id] = new_id

            cols[8] = cols[8].replace("ID={0}".format(id), "ID={0}".format(new_id))

        if parent is not None:
            if parent in id_map:
                new_parent = id_map[parent]
            else:
                raise Exception("ERROR: parent ({0}) referenced before it was used as an ID".format(parent))

            cols[8] = cols[8].replace("Parent={0}".format(parent), "Parent={0}".format(new_parent))

        #print("DEBUG: old_id:{0} - old_parent:{1}, new_id:{2} - new_parent:{3}".format(id, parent, new_id, new_parent))
        fout.write("\t".join(cols) + "\n")

    #>>> binascii.hexlify(os.urandom(4))
    #b'c08ca446'

    #>>> uuid.uuid4()
    #UUID('37cd0fbf-bdc3-49bc-8351-a7ebc5a93ea5')


def get_new_id(prefix, type, mode):
    new_id = "{0}.{1}.".format(prefix, type)
    
    if mode == 'sequential':
        new_id += str(next_ids_sequential[type])
        next_ids_sequential[type] += 1
    elif mode == 'uuid':
        new_id += str(uuid4())
    elif mode == 'hex8':
        new_id += hexlify(os.urandom(4)).decode('ascii')
    elif mode == 'hex12':
        new_id += hexlify(os.urandom(6)).decode('ascii')

    new_id += '.1'

    return new_id

def check_arguments( args ):
    # Check the acceptable values for format
    mode_options = ('sequential', 'uuid', 'hex8', 'hex12')
    if args.mode not in mode_options:
        raise Exception("ERROR: The --mode provided ({0}) isn't supported.  Please check the documentation again.".format(args.mode))




if __name__ == '__main__':
    main()







