#!/usr/bin/env python3

"""

Repeatedly asks for user input and then reports whether that could be written
in the alphabet of single-letter amino acid codes.

Selenocysteine is included.

Assumes we're upper-casing all input

"""

import os
import sys


def main():
    codes = 'ARNDCQEGHILKMFPSTWYVU'

    while True:
        user_input = input("Enter text to test (q to quit): ")
        user_input = user_input.rstrip()

        if user_input == 'q':
            sys.exit(0)

        user_input = user_input.upper()

        bad_letters = 0

        for letter in user_input:
            if letter == ' ':
                continue
            elif letter not in codes:
                print(f"Not an amino acid letter: {letter}")
                bad_letters += 1

        if bad_letters == 0:
            print("Legal amino acid string!\n")
        else:
            print("\n")
        

if __name__ == '__main__':
    main()







