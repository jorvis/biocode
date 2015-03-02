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
import biocodeutils
import datetime
import os
from queue import Queue
import re
import subprocess
import threading
import time

# lock to serialize console output
lock = threading.Lock()
q = Queue()

def main():
    parser = argparse.ArgumentParser( description='Merge multiple BAM files, retaining headers')

    ## output file to be written
    parser.add_argument('-i', '--input_list', type=str, required=True, help='Path to an input list file to be read' )
    parser.add_argument('-o', '--output_base', type=str, required=True, help='Base name/path of output file to be created' )
    parser.add_argument('-t', '--threads', type=int, required=False, default=1, help='Count of threads to use' )
    args = parser.parse_args()

    ## duh.
    if args.threads < 1:
        raise Exception("ERROR: the thread count must be >= 1")

    bam_files = biocodeutils.read_list_file(args.input_list)

    # Create the queue and thread pool.
    for i in range(args.threads):
         t = threading.Thread(target=worker)
         t.daemon = True  # thread dies when main thread (only non-daemon thread) exits.
         t.start()

    # stuff work items on the queue (in this case, just a number).
    start = time.perf_counter()
    for file_path in bam_files:
        q.put(file_path)

    q.join()       # block until all tasks are done

    ## take the header from the first file, change the sorting definition
    first_file = bam_files[0]
    ofh = open("{0}.sam".format(args.output_base), 'wt')
    for line in open(first_file):
        ## we expect: @HD     VN:1.0  SO:coordinate
        m = re.match("^@HD\s+VN\:(\S+)\s+SO:", line)
        if m:
            ofh.write("@HD\tVN:{0}\tSO:unknown\n".format(m.group(1)))
        break
    ofh.close()

    ## merge headers
    for bam_file in bam_files:
        run_command("grep -v '@HD' {0}.header >> {1}.sam".format(bam_file, args.output_base) )

    ## merge sequence file
    for bam_file in bam_files:
        run_command("cat {0}.alignments >> {1}.sam".format(bam_file, args.output_base) )

    ## convert SAM to BAM
    run_command("samtools view -bS {0}.sam > {0}.bam".format(args.output_base))

    # "Work" took .1 seconds per task.
    # 20 tasks serially would be 2 seconds.
    # With 4 threads should be about .5 seconds (contrived because non-CPU intensive "work")
    print('time:',time.perf_counter() - start)

def do_work(bam_path):
    # Make sure the whole print completes or threads can mix up output in one line.
    with lock:
        print("{0}\tProcessing: {1}".format(threading.current_thread().name, bam_path), flush=True)

    # generate the header file
    run_command("samtools view -H {0} > {0}.header".format(bam_path))

    # generate the alignments file
    run_command("samtools view {0} > {0}.alignments".format(bam_path))

    with lock:
        print("{0}\tCompleted: {1}".format(threading.current_thread().name, bam_path), flush=True)

# The worker thread pulls an item from the queue and processes it
def worker():
    while True:
        bam_path = q.get()
        do_work(bam_path)
        q.task_done()

def run_command(cmd):
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
       raise Exception("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))

if __name__ == '__main__':
    main()







