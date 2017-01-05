#!/usr/bin/env python3

"""

OUTPUT:
Using the base name passed, several files will be created temporarily, but the final output
path will be:

    OUTPUT_BASE.bam

ASSUMPTIONS:
This assumes very simple headers only as input like this:

    @HD     VN:1.0  SO:coordinate
    @SQ     SN:gi|551236890|ref|NZ_ATVB01000022.1|  LN:17321
    @SQ     SN:gi|657682046|ref|NZ_JNIV01000093.1|  LN:18254
    @SQ     SN:gi|655339533|ref|NZ_ATYY01000023.1|  LN:86451
    @SQ     SN:gi|654719175|ref|NZ_AXAT01000007.1|  LN:74141

The first line will be checked for the version number, the sort order will be changed to 'unknown'.

It's probably best to redo a sort operation after running this, such as:

    samtools sort test.bam test_sorted

"""

import argparse
import datetime
import re
import subprocess

from biocode import utils


def main():
    parser = argparse.ArgumentParser( description='Merge multiple BAM files, retaining headers')

    ## output file to be written
    parser.add_argument('-i', '--input_list', type=str, required=True, help='Path to an input list file to be read' )
    parser.add_argument('-o', '--output_base', type=str, required=True, help='Base name/path of output file to be created' )
    args = parser.parse_args()

    bam_files = utils.read_list_file(args.input_list)

    ## take the header from the first file, change the sorting definition
    first_file = bam_files[0]
    run_command("samtools view -H {0} | head -n 1 > {0}.firstline".format(first_file))
    
    ofh = open("{0}.sam".format(args.output_base), 'wt')
    for line in open("{0}.firstline".format(first_file)):
        ## we expect: @HD     VN:1.0  SO:coordinate
        m = re.match("^@HD\s+VN\:(\S+)\s+SO:", line)
        if m:
            ofh.write("@HD\tVN:{0}\tSO:unknown\n".format(m.group(1)))
    ofh.close()

    ## merge headers
    for bam_file in bam_files:
        run_command("samtools view -H {0} | grep -v '@HD' >> {1}.sam".format(bam_file, args.output_base))

    ## merge sequence files
    for bam_file in bam_files:
        run_command("samtools view {0} >> {1}.sam".format(bam_file, args.output_base))

    ## convert SAM to BAM
    run_command("samtools view -bS {0}.sam > {0}.bam".format(args.output_base))

    ## delete the SAM
    run_command("rm {0}.sam".format(args.output_base))

def run_command(cmd):
    print("INFO: Running command: {0}",format(cmd), flush=True)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
       raise Exception("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))

if __name__ == '__main__':
    main()







