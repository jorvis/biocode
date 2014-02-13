#!/usr/bin/env python3

'''
Loads features from a test GFF3 file and uses them to test the coordinate comparison
functions within the biothings API.

Author:  Joshua Orvis
'''

import argparse
import os
import biocodegff
import biothings

def main():
    bin_dir = os.path.abspath(os.path.dirname(__file__))
    test_gff_file = bin_dir + '/biothings_coordinate_comparisons.data'
    
    (assemblies, features) = biocodegff.get_gff3_features( test_gff_file )


###########################################################################################

    if features['TP03_0010'] < features['TP03_0012.t01_polypeptide']:
        print("INFO: < positive check successful")
    else:
        print("ERROR: < check unsuccessful")

    if features['TP03_0012'] < features['TP03_0012.t01_polypeptide']:
        print("ERROR: < check unsuccessful")
    else:
        print("INFO: < negative check successful")

###########################################################################################
        
    if features['TP03_0012'] > features['TP03_0010']:
        print("INFO: > positive check successful")
    else:
        print("ERROR: > check unsuccessful")

    if features['TP03_0010'] > features['TP03_0012.t01_polypeptide']:
        print("ERROR: > check unsuccessful")
    else:
        print("INFO: > negative check successful")
        
###########################################################################################

    if features['TP03_0012.t01_exon-auto15079'] <= features['TP03_0012.t01_polypeptide']:
        print("INFO: <= positive check successful")
    else:
        print("ERROR: <= check unsuccessful")

    if features['TP03_0010'] <= features['TP03_0012']:
        print("ERROR: <= check unsuccessful")
    else:
        print("INFO: <= negative check successful")

###########################################################################################

    if features['TP03_0012.t01_exon-auto15085'] >= features['TP03_0012.t01_polypeptide']:
        print("INFO: >= positive check successful")
    else:
        print("ERROR: >= check unsuccessful")

    if features['TP03_0010'] >= features['TP03_0012']:
        print("ERROR: >= check unsuccessful")
    else:
        print("INFO: >= negative check successful")

###########################################################################################

    if features['TP03_0012.t01_exon-auto15079'].overlaps_with(features['TP03_0012.t01_polypeptide']):
        print("INFO: overlaps_with() positive check successful")
    else:
        print("ERROR: overlaps_with() positive check unsuccessful")

    if features['TP03_0002'].overlaps_with(features['TP03_0010']):
        print("ERROR: overlaps_with() negative check unsuccessful")
    else:
        print("INFO: overlaps_with() negative check successful")


if __name__ == '__main__':
    main()






