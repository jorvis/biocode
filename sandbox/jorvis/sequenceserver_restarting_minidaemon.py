#!/usr/bin/env python3

"""

The SequenceServer software is great but dies all the time. Whether from unexpected sequence
composition or whatever, this just logs and restarts it.

"""

import argparse
import os
import subprocess


def main():
    parser = argparse.ArgumentParser( description='SequenceServer restarter')
    parser.add_argument('-m', '--mode', type=str, required=True, help='Mode is either legacy, genomic or transcriptomic' )
    args = parser.parse_args()

    while True:

        if args.mode == 'legacy':
            cmd = "sequenceserver -c ~/.sequenceserver.old.conf"
        elif args.mode == 'genomic':
            cmd = "/opt/sequence-server2/bin/sequenceserver -c ~/.sequenceserver.genomic.conf"
        elif args.mode == 'transcriptomic':
            cmd = "sequenceserver"
        else:
            raise Exception("Unrecognized --mode value. Check the usage")

        print("Executing command {0}".format(cmd))
        return_code = subprocess.call(cmd, shell=True)
        print("Restarting ...")

if __name__ == '__main__':
    main()







