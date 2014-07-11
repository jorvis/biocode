#!/usr/bin/env python3

"""
OVERVIEW

Some analysis tasks produce an interleaved FASTA or FASTQ file (such as digital normalization).
Use this script to split that interleaved file back out into R1, R2 and singleton files.

The read mates will be kept in pairwise order within the R1 and R2 files.

The -o option defines the base name of the output files to which the directionality
and extensions will be added.  So given the options:

  -o ATCC_30222
  -t fastq

The following files will be created:

  ATCC_30222.R1.fastq
  ATCC_30222.R1.fastq
  ATCC_30222.single.fastq

Note: No format conversion will happen here.  If your input is FASTA your output
will also be FASTA.

HEADERS

This script handles read headers in either of the two following formats:

   @SN7001163:74:C0YGBACXX:1:1101:1099:2196 1:N:0:ATCACGA
   @SN7001163:74:C0YGBACXX:1:1101:1099:2196/1

The first of these is what comes off our Illumina machines, with the direction being the first
number after the whitespace.  The latter format with /1 or /2 is what most software expects.

READ ORDER

This script expects that any reads which have pairs will have the /1 direction first immediately
followed by the /2 mate.  If one of these is absent the singlet will be written to the 'single'
file.  If the /1 and /2 mates are not consecutive, they will both be written as singlets.

Original author: Priti Kumari
Heavy edits: Joshua Orvis
"""

import argparse
import os
import re


def interface():
    ## read input fastq file
    parser = argparse.ArgumentParser( description='Script to split output of Diginorm.')
    parser.add_argument('-i', '--input_file',type=str, required=True, help='The interleaved fastq/fasta file to split.')
    parser.add_argument('-t', '--input_type',type=str, required=True, help='The type of input file(fasta/fastq).')
    ## output file to be written
    parser.add_argument('-o', '--output', type=str, required=True, help='Base path/name for the output files to be created' )
    parser = parser.parse_args()
    return parser

def process_fasta(args):
    read1_file  = args.output + ".R1.fasta"
    read2_file  = args.output + ".R2.fasta"
    single_file = args.output + ".single.fasta"
    
    if args.input_file.endswith('.gz'):
        inp = gzip.open(args.input_file,'rb')
        fout1= gzip.open(read1_file,'wb')
        fout2= gzip.open(read2_file,'wb')
        fout3= gzip.open(single_file,'wb')
    else:
        inp =  open(args.input_file,'rU')
        fout1= open(read1_file,'w')
        fout2= open(read2_file,'w')
        fout3= open(single_file,'w')

    flag = 0
    while 1 :
        if flag == 0:
            Read1=''
        Read2=''
        c=2
        d=2
        if flag == 0 :
            while (c > 0):
                 inp_line = inp.readline()
                 Read1 += inp_line
                 c = c-1
            if Read1 == '': break
            id1 =  Read1.split('\n')[0].split('/')[0]
            if not id1.startswith('>'):
                print ("Error!! Unknown header: not '>'")
                break
            tag1 = Read1.split('\n')[0].split('/')[1]
        while (d > 0):
            inp_line = inp.readline()           
            Read2 += inp_line
            d = d-1
        if not Read1 == '' and  Read2 == '':
            fout3.write(Read1)
            break
        
        id2  = Read2.split('\n')[0].split('/')[0]
        if not id1.startswith('>'):
                print ("Error!! Unknown header: not '>'")
                break
            
        tag2 = Read2.split('\n')[0].split('/')[1]
        if tag1 == '1' and tag2 == '2' :
            if id1 == id2:
                fout1.write(Read1)
                fout2.write(Read2)
                flag=0
            else:
                fout3.write(Read1)
                fout3.write(Read2)
                flag=0
                
        elif tag1 == '1' and tag2 == '1':
            fout3.write(Read1)
            flag=1
            Read1=Read2
            tag1=tag2
            id1=id2
            
        elif tag1 == '2' and tag2 == '2':
            fout3.write(Read1)
            fout3.write(Read2)
            flag=0
            
        elif tag1 == '2' and tag2 == '1':
            fout3.write(Read1)
            Read1=Read2
            tag1=tag2
            id1=id2
            flag=1
        

def process_fastq(args):
    read1_file  = args.output + ".R1.fastq"
    read2_file  = args.output + ".R2.fastq"
    single_file = args.output + ".single.fastq"

    input_is_compressed = False
    
    if args.input_file.endswith('.gz'):
        inp = gzip.open(args.input_file,'rb')
        read1_out_fh = gzip.open(read1_file,'wb')
        read2_out_fh = gzip.open(read2_file,'wb')
        singles_out_fh = gzip.open(single_file,'wb')
        input_is_compressed = True
    else:
        inp =  open(args.input_file,'rU')
        read1_out_fh = open(read1_file,'w')
        read2_out_fh = open(read2_file,'w')
        singles_out_fh = open(single_file,'w')

    line_count = 0
    last_base_name = None
    last_dir = None
    last_base_buffer = list()
    current_fh = None
        
    for line in inp:
        if input_is_compressed:
            line = line.decode()

        line_count += 1

        if re.match(r'^\s*$', line): continue

        if line_count % 4 == 1:
            # this should be a header line
            spaced_match = re.search(r'^\@(\S+) ([12])', line)
            
            if spaced_match:
                base = spaced_match.group(1)
                dir = spaced_match.group(2)
            else:
                traditional_match = re.search(r'^\@(\S+)\/([12])', line)
                if traditional_match:
                    base = traditional_match.group(1)
                    dir = traditional_match.group(2)
                else:
                    raise Exception("ERROR: Expected this to be a header line, but format wasn't recognized: {0}".format(line))

            ## if this is /2, was /1 found?
            if int(dir) == 1:
                # if the last direction was 1, purge it as a singleton
                if last_dir == 1:
                    for buffer_line in last_base_buffer:
                        singles_out_fh.write(buffer_line)
                
                last_base_name = base
                last_dir = 1
                last_base_buffer = [line]
                current_fh = None
                
            elif int(dir) == 2:
                if last_base_name == base:
                    if last_dir == 1:
                        # purge the /1
                        for buffer_line in last_base_buffer:
                            read1_out_fh.write(buffer_line)

                        # then write the /2
                        read2_out_fh.write(line)
                        current_fh = read2_out_fh
                    else:
                        raise Exception("ERROR: Were there two {0}/2 reads in a row?".format(last_base_name))
                else:
                    # this must be a /2 where the /1 is missing
                    singles_out_fh.write(line)
                    current_fh = singles_out_fh

                last_base_name = base
                last_dir = 2
        else:
            if current_fh is None:
                last_base_buffer.append(line)
            else:
                current_fh.write(line)


if __name__ == '__main__':
    args = interface()
    if args.input_type == 'fastq':
        process_fastq(args)
    elif args.input_type == 'fasta':
        process_fasta(args)
    else:
        print ("Error! Input file format incorrect.  Expected to be 'fastq' or 'fasta'");

