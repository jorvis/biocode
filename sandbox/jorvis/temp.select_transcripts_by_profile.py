#!/usr/bin/env python3

"""
This script is currently quite experimental, and requires editing to use for each organism.  If
you can't see what to edit in the code you probably shouldn't be trying it.

Your input is a GFF3 file containing gene models.

"""

import argparse
import math
import random

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Generates a set of transcripts based on a user-defined exon-complexity profile')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output ID list file to create' )
    parser.add_argument('-ni', '--not_included', type=str, required=False, help='Writes the ID list of genes not included' )
    parser.add_argument('-c', '--count', type=int, required=True, help='Count of transcripts to pull' )
    parser.add_argument('-e', '--exclude', type=str, required=False, help='List of IDs to exclude' )
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_file)
    exclude = list()

    if args.exclude is not None:
        for line in open(args.exclude):
            line = line.rstrip()
            exclude.append(line)

    ofh = open(args.output_file, 'wt')

    if args.not_included is None:
        ni_ofh = None
    else:
        ni_ofh = open(args.not_included, 'wt')



    # profile for 99-892
    #profile = { 1:12.2, 2:13.2, 3:15.4, 4:14.1, 5:11.5, 6:8.77, 7:6.99, 8:4.74, 9:3.81, 10:2.45 }

    # profile for 99-880
    #profile = { 1:19.7, 2:17.9, 3:20.6, 4:15.6, 5:10.1, 6:6.62, 7:4.14, 8:1.58, 9:2.07 }

    # c. hominis TU502
    #profile = { 1:89.7, 2:6.8, 3:2.6, 4:0.78 }

    # c. baileyi TAMU 10GZ1
    #profile = { 1:85.6 , 2:9.85 , 3:4.06 , 4:0.25 , 5:0.25 }
    profile = { 1:85.6 , 2:9.85 , 3:4.06 }
    
    
    mRNAs = dict()
    unselected_mRNAs = dict()
    target = dict()
    selected = dict()

    for CDS_count in profile:
        target[CDS_count] = math.trunc( args.count * (profile[CDS_count]/100) )
        selected[CDS_count] = list()
        mRNAs[CDS_count] = list()
        unselected_mRNAs[CDS_count] = list()

    # fill the bins of each target size, then fill a reservoir to select any additional ones from
    reservoir = list()
    total_mRNAs_selected = 0
    total_mRNA_count = 0
        
    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            for mRNA in gene.mRNAs():
                total_mRNA_count += 1
                
                if mRNA.id in exclude:
                    continue
                
                CDS_count = mRNA.CDS_count()

                if CDS_count not in profile:
                    reservoir.append(mRNA)
                    continue

                mRNAs[CDS_count].append(mRNA)

    for CDS_count in profile:
        # make sure this many were found
        if target[CDS_count] <= len(mRNAs[CDS_count]):
            selected[CDS_count] = random.sample( mRNAs[CDS_count], target[CDS_count] )
        else:
            print("WARN: Not enough mRNAs of length {0} to meet profile request".format(CDS_count))
            selected[CDS_count] = mRNAs[CDS_count]
        
        total_mRNAs_selected += len(selected[CDS_count])
        unselected_mRNAs[CDS_count] = list(set(mRNAs[CDS_count]) & set(set(mRNAs[CDS_count]) ^ set(selected[CDS_count])))

        for mRNA in unselected_mRNAs[CDS_count]:
            reservoir.append(mRNA)

    print("INFO: selected CDS profile:")
    for ccount in sorted(selected):
        print("CDS_count:{0}, target:{3}, gathered:{1}, unselected:{4}, target_perc:{2}".format( \
                ccount, len(selected[ccount]), (len(selected[ccount])/target[ccount]), target[ccount], \
                len(unselected_mRNAs[ccount]) \
              ) )

        for mRNA in selected[ccount]:
            ofh.write("{0}\n".format(mRNA.parent.id) )

    print("Total selected according to profile: {0}".format(total_mRNAs_selected) ) 
    # now, from the rounding portions fill the rest randomly
    sample_from_reservoir = random.sample(reservoir, args.count - total_mRNAs_selected)
    sample_ids_from_reservoir = list()

    for mRNA in sample_from_reservoir:
        sample_ids_from_reservoir.append(mRNA.id)
        reservoir.remove(mRNA)
    
    ofh.write( "\n".join( sample_ids_from_reservoir ) )
    ofh.write("\n")

    total_mRNAs_selected += len(sample_from_reservoir)

    print("Total selected randomly afterwards: {0}".format(len(sample_from_reservoir)) )

    if ni_ofh is not None:
        for mRNA in reservoir:
            ni_ofh.write("{0}\n".format(mRNA.parent.id) )
        
        

    
                
        
if __name__ == '__main__':
    main()







