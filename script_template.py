#!/usr/bin/env python3.2

import argparse
import os


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()




if __name__ == '__main__':
    main()







