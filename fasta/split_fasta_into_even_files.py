#!/usr/bin/env python3

"""

If you have a large FASTA file with a great number of entries, this script is used to
split the file into multiple new files with an even number of records (or as close
as possible.)

Use case example: I have a reference genome collection with 2000 genomes and want to
search this using bowtie.  But the bowtie index memory footprint is too large for
a large search, so I needed to split the reference genome file into 4 even parts
and index them individually.

Usage example:

./split_fasta_into_even_files.py -i reference_genomes.fna -fc 4

The script uses the same file name as before but appends '.partN' where 'N' is replaced
by an increasing digit to indicate the file number.

"""

import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser( description='Split a large FASTA file into new evenly sized files')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input FASTA file to be read' )
    parser.add_argument('-fc', '--file_count', type=int, required=True, help='Number of files to split the records into' )
    parser.add_argument('-o', '--output_directory', type=str, required=False, help='Directory where output will be written.  Default: input file location' )
    args = parser.parse_args()

    # First, we need to know how many entries are in the file
    total_record_count = 0
    for line in open(args.input_file):
        if line.startswith(">"):
            total_record_count += 1

    # It doesn't make sense for there to be fewer records than file count passed
    if total_record_count < args.file_count:
        raise Exception("Error: You asked for {0} split files to be created but there were only {1} " \
                        "entries in the input file.".format(args.file_count, total_record_count))

    if args.output_directory:
        output_dir = args.output_directory
    else:
        output_dir = os.path.dirname(args.input_file)

    if output_dir is None or len(output_dir) == 0:
        output_dir = "./"

    print("INFO: Output directory is: ({0})".format(output_dir))

    basename = os.path.basename(args.input_file)
                        
    print("INFO: There were {0} records found in the input file.".format(total_record_count))
    min_records_per_file = int(total_record_count / args.file_count)
    file_count_to_create = int(total_record_count / min_records_per_file)

    # words matter
    if file_count_to_create == 1:
        print("INFO: {0} file will be created.".format(file_count_to_create))
    else:
        print("INFO: {0} files will be created.".format(file_count_to_create))

    print("INFO: Most files will have {0} records in each.".format(min_records_per_file))

    file_part_num = 1
    current_fragment_record_count = 0
    current_fh = open("{0}/{1}.part{2}".format(output_dir, basename, file_part_num), 'w')

    for line in open(args.input_file):
        if line.startswith(">"):
            current_fragment_record_count += 1
        
        if current_fragment_record_count == min_records_per_file and file_part_num < file_count_to_create:
            if file_part_num <= args.file_count:
                if current_fh is not None:
                    current_fh.close()

                file_part_num += 1
                current_fh = open("{0}/{1}.part{2}".format(output_dir, basename, file_part_num), 'w')
                current_fragment_record_count = 0

        current_fh.write(line)

    current_fh.close()


if __name__ == '__main__':
    main()







