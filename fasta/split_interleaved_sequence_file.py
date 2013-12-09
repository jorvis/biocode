#!/usr/bin/env python3

"""
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

Warning:  This script requires that the read IDs end with /1 or /2 to discern their direction.
  
Original author: Priti Kumari
"""



import argparse
import os


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
        c=4
        d=4
        if flag == 0 :
            while (c > 0):
                 inp_line = inp.readline()
                 Read1 += inp_line
                 c = c-1
            if Read1 == '': break
            id1 =  Read1.split('\n')[0].split(' ')[0]
            if not id1.startswith('@'):
                print ("Error!! Unknown header: not '@'")
                break
            
            tag1 = Read1.split('\n')[0].split(' ')[1].split(':')[0]
        while (d > 0):
            inp_line = inp.readline()
            Read2 += inp_line
            d = d-1
        if not Read1 == '' and  Read2 == '':
            fout3.write(Read1)
            break
            
        id2  = Read2.split('\n')[0].split(' ')[0]
        if not id2.startswith('@'):
                print ("Error!! Unknown header: not '@'")
                break
        tag2 = Read2.split('\n')[0].split(' ')[1].split(':')[0]
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
            

if __name__ == '__main__':
    args = interface()
    if args.input_type == 'fastq':
        process_fastq(args)
    elif args.input_type == 'fasta':
        process_fasta(args)
    else:
        print ("Error! Input file format incorrect.");

