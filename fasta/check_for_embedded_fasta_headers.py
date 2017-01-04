#!/usr/bin/env python3

"""
There are occasionally issues with large FASTA files in which one entry gets truncated
and the next entry starts immediately on the next character of the line, causing
merged entries.  Example:

    >gi1006569 DnaA [Vibrio harveyi]
    MSSSLWLQCLQQLQEELPATEFSMWVRPLQAELNDNTLTLFAPNRFVLDWVRDKYLNSINRLLQEYCGND
    VPSLRFEVGSRRVAAPKPAPTRTPADVAAESSAPAQLQARKPVHKTWDDDAQVIADINHRSNVNPKHKFN
    NFVEGKSNQLGLAAARQVSDNPGAAYNPLFLYGGTGLGKTHLLHAVGNAIVDNNPNAKVVYMHSERFVQD
    MVKALQNNAIEEFKRYYRSVDALLIDDIQFFANKERSQEEFFHTFNALLEGNQQIILTSDRYPKEISGVE
    DRLKSRFGWGLTVAIEPPELETRVAILMKKAEDHQIHLADEVAFFIAKRLRSNVRELEGALNRVIANANF
    TGRPITIDFVREALRDLLALQEKLVTIDNIQKTVAEYYKIKVADLLSKRRSRSVARPRQLAMALAKELTN
    HSLPEIGDAFGGRDHTTVLHAC>gi409247 DNA repair protein
    MVSLTFKNFKKEKVPLDLEPSNTILETKTKLAQSISCEESQIKLIYSGKVLQDSKTVSECGLKDGDQVVF
    MVSQKKSTKTKVTEPPIAPESATTPGRENSTEASPSTDASAAPAATAPEGSQPQEEQTATTERTESASTP
    GFVVGTERNETIERIMEMGYQREEVERALRAAFNNPDRAVEYLLMGIPENLRQPEPQQQTAAAAEQPSTA
    ATTAEQPAEDDLFAQAAQGGNASSGALGTTGGATDAAQGGPPGSIGLTVEDLLSLRQVVSGNPEALAPLL
    ENISARYPQLREHIMANPEVFVSMLLEAVGDNMQDVMEGADDMVEGEDIEVTGEAAAAGLGQGEGEGSFQ
    VDYTPEDDQAISRLCELGFERDLVIQVYFACDKNEEAAANILFSDHAD
    >gi4587293 RecN [Deinococcus radiodurans]
    MVEPPPPDAAPTGPRLSRLEIRNLATITQLELELGGGFCAFTGETGAGKSIIVDALGLLLGGRANHDLIR
    SGEKELLVTGFWGDGDESEADSASRRLSSAGRGAARLSGEVVSVRELQEWAQGRLTIHWQHSAVSLLSPA
    NQRGLLDRRVTKEAQAYAAAHAAWREAVSRLERLQASQRERARQIDLLAFQVQEISEVSPDPGEEEGLNT
    ELSRLSNLHTIAQAAAGGVELLSDGDLNAAGLIGEAVRALNAGAKYDETVMQLQNELRAALESVQAIAGE
    LRDVAEGSAADPEALDRVEARLSALSKLKNKYGPTLEDVVEFGAQAAEELAGLEEDERDAGSLQADVDAL
    HAELLKVGQALDAAREREAEPLVDSLLAVIRELGMPHARMEFALSALAEPAAYGLSDVLLRFSANPGEEL
    GPLSDVASGGELSRVMLAVSTVLGADTPSVVFDEVDAGIGGAAAIAVAEQLSRLADTRQVLVVTHLAQIA
    ARAHHHYKVEKQVEDGRTVSHVRLLTGDERLEEIARMLSGNTSEAALEHARELLAG


Running on this file gives:

$ ./check_for_embedded_fasta_headers.py -f /tmp/test.faa 
INFO: processing file: /tmp/test.faa
/tmp/test.faa	1

	HSLPEIGDAFGGRDHTTVLHAC>gi409247 DNA repair protein
"""

import argparse
import sys

from biocode.utils import read_list_file


def main():
    parser = argparse.ArgumentParser( description='Looks for FASTA entries which seem to have been embedded within another')
    parser.add_argument('-l', '--input_list', type=str, required=False, help='A list file with the paths of all files' )
    parser.add_argument('-f', '--fasta_file', type=str, required=False, help='The single FASTA file to be processed' )
    args = parser.parse_args()

    if args.input_list is None and args.fasta_file is None:
        raise Exception("ERROR: You must specify either --input_list or --fasta_file to define input");

    files_to_process = list()

    if args.input_list is not None:
        for file in read_list_file( args.input_list ):
            files_to_process.append(file)

    if args.fasta_file is not None:
        files_to_process.append(args.fasta_file)

    for file in files_to_process:
        fail_lines = process_fasta_file( file )

        if ( len(fail_lines) > 0 ):
            print( "{0}\t{1}\n".format(file, len(fail_lines)) )

            for line in fail_lines:
                print("\t{0}\n".format(line))

        


def process_fasta_file( file ):
    embedded_header_lines = list()

    print("INFO: processing file: {0}".format(file), file=sys.stderr )

    ifh = open(file)
    
    for line in ifh:
        ## if the line contains a '>' but it's not at the beginning, flag the entry:
        if line.find('>') != -1 and not line.startswith('>'):
            embedded_header_lines.append(line)
    
    
    return embedded_header_lines


if __name__ == '__main__':
    main()







