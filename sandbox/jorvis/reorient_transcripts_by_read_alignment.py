#!/usr/bin/env python3

"""

Put some general, high-level documentation here

http://onetipperday.blogspot.com/2012/04/understand-flag-code-of-sam-format.html

Better detailed version:
http://blog.nextgenetics.net/?e=18

"""

import argparse
import os


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    #parser.add_argument('-fi', '--fasta_in', type=str, required=False, help='Path to a FASTA file representing sequences that were aligned against.  If this is passed, you should also pass the -fo argument' )
    #parser.add_argument('-fo', '--fasta_out', type=str, required=False, help='If passed along with -fi, the orientation-corrected sequences will be written here.' )
    args = parser.parse_args()

    total_read_mappings = 0
    last_transcript_id = None
    counts = { '1':{'T':0,'F':0}, '2':{'T':0,'F':0} }
    transcript_count = 0
    correct_orientation_count = 0
    incorrect_orientation_count = 0

    for line in open(args.input_file):
        if line.startswith('@'): continue
        
        cols = line.split("\t")
        if len(cols) < 5: continue

        read_dir = cols[0][-1]
        transcript_id = cols[2]
        total_read_mappings += 1

        flag = cols[1]
        if int(flag) & 16:
            seq_revcomped = 'T'
        else:
            seq_revcomped = 'F'

        #print("DEBUG: match:{2}, SEQ_revcomped={0}, read_dir={1}".format(seq_revcomped, read_dir, transcript_id))

        if transcript_id == last_transcript_id:
            counts[read_dir][seq_revcomped] += 1
        else:
            transcript_count += 1
            
            if last_transcript_id is not None:
                ## determine transcript orientation
                ## Given an RF library, the 1:T count should outnumber the 1:F one
                if counts['1']['T'] > counts['1']['F']:
                    correct_orientation_count += 1
                else:
                    incorrect_orientation_count += 1
                
                ## report counts
                #print("{0}\t1-T:{1}\t1-F:{2}\t2-T:{3}\t2-F:{4}".format(last_transcript_id, counts['1']['T'], counts['1']['F'], counts['2']['T'], counts['2']['F']))

            ## reset
            last_transcript_id = transcript_id
            counts = { '1':{'T':0,'F':0}, '2':{'T':0,'F':0} }


    print("Total transcripts: {0}".format(transcript_count))
    print("Total reads mapped: {0}".format(total_read_mappings))
    print("Transcripts in correct orientation: {0}".format(correct_orientation_count))
    print("Transcripts in reverse orientation: {0}".format(incorrect_orientation_count))
#    print("R1 reads mapped where SEQrevcomp=FALSE: {0}".format(dir1_reads_mapped))
#    print("R2 reads mapped where SEQrevcomp=TRUE: {0}".format(dir2_reads_mapped))
    
        
    



if __name__ == '__main__':
    main()







