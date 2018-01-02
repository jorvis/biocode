#!/usr/bin/env python3

"""

This is a utility for users who have installed Biocode via PyPi to 
see the list of scripts which are available to them.  Prints categories first
and then the list of scripts under each category.

"""

import argparse
import json
from operator import itemgetter
import os


def main():
    parser = argparse.ArgumentParser( description='Prints a list of the scripts which are part of the Biocode install')
    args = parser.parse_args()

    index_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/biocode_script_index.json')
    print("Opening this file: {0}".format(index_path))

    if not os.path.exists(index_path):
        error_msg = "ERROR: Failed to find expected index path ({0}). ".format(index_path)
        error_msg += "This is built dynamically when installed via PyPi.  If you got biocode from GitHub, this utility won't work."
        raise Exception(error_msg)

    with open(index_path) as f:
        json_data = json.loads(f.read())

    for cat in sorted(json_data):
        print("{0}:".format(cat))

        for item in sorted(json_data[cat], key=itemgetter('name')):

            print("\t{0}\n\t\t{1}".format(item['name'], item['desc']))
    

if __name__ == '__main__':
    main()







