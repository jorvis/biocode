#!/usr/bin/env python3

"""



"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    tpms = list()
    
    with open(args.input_file) as fh:
        # read off the header line
        fh.readline()
        
        for line in fh:
            line = line.rstrip()
            cols = line.split("\t")
            tpms.append(float(cols[4]))

    # stores the count of transcripts at each TPM value or lower
    tpm_c = dict()
    last_tpm = None
    tpm_count = 0

    for tpm in sorted(tpms, reverse=True):
        tpm_count += 1

        if tpm == last_tpm:
            tpm_c[tpm] += 1
        else:
            tpm_c[tpm] = tpm_count

        last_tpm = tpm

    tpm_for_df = dict({'tpm':[], 'count':[]})
    
    for tpm in tpm_c:
        tpm_for_df['tpm'].append(tpm)
        tpm_for_df['count'].append(tpm_c[tpm])
        print("{0}\t{1}".format(tpm, tpm_c[tpm]))

    df = pd.DataFrame(tpm_for_df)

    plt.figure(figsize=(20,9)) 
    lineplot = sns.lineplot(data=df, x='tpm', y='count')
    lineplot.set(xscale="log")
    fig = lineplot.get_figure()
    fig.savefig(args.output_file)
    
        

if __name__ == '__main__':
    main()







