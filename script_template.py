#!/usr/bin/env python3

"""

Put some general, high-level documentation here

"""

import argparse
import os


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## Add arguments for your script here.  Already added are -i and -o for input/output.
    #  More information on argument parsing here: https://docs.python.org/3/howto/argparse.html
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    # Your code goes here


if __name__ == '__main__':
    main()







