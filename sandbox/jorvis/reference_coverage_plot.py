#!/usr/bin/env python3

"""

Initially written to plot results from Priti's pipeline, where the goal was to determine
how well an assembled transcriptome (Trinity) produced contigs which covered a set of
reference transcripts.  Steps in her pipeline:

## First create symlinks:

reference.fasta - The smaller set of reference transcripts for which you want to calculate coverage
assembly.fasta - The transcriptome assem

formatdb -i assembly.fasta -p F
blastall -p blastn -d assembly.fasta -i reference.fasta -m 9 -e 1 -o assembly.btab -a 8
grep -v "#" assembly.btab | awk '{print $1"\t"$7"\t"$8"\t"$2"\t0\t+"}' > assembly.bed
bedtools sort -i assembly.bed > assembly.sorted.bed
bedtools merge -i assembly.sorted.bed > assembly.merged.bed
~/git/biocode/fasta/fasta_size_report.pl -f reference.fasta | cut -f 1,2 > transcript.lengths
#/local/projects-t3/aplysia/priti_analysis/Scripts/comp_transcript.pl --i assembly.merged.bed --i1 transcript.lengths > assembly.all.per_cov.txt
/usr/local/projects/aplysia/priti_analysis/Scripts/comp_transcript.pl --i assembly.merged.bed --i1 transcript.lengths > assembly.all.per_cov.txt

#/local/projects-t3/aplysia/priti_analysis/Scripts/longest_transcript.pl --i assembly.btab > assembly.longest.btab
/usr/local/projects/aplysia/priti_analysis/Scripts/longest_transcript.pl --i assembly.btab > assembly.longest.btab
grep -v "#" assembly.longest.btab | awk '{print $1"\t"$7"\t"$8"\t"$2"\t0\t+"}' > assembly.longest.bed
bedtools sort -i assembly.longest.bed > assembly.longest.sorted.bed
#/local/projects-t3/aplysia/priti_analysis/Scripts/comp_transcript.pl --i assembly.longest.sorted.bed --i1 transcript.lengths > assembly.longest.per_cov.txt
/usr/local/projects/aplysia/priti_analysis/Scripts/comp_transcript.pl --i assembly.longest.sorted.bed --i1 transcript.lengths > assembly.longest.per_cov.txt

# then to actually plot these:
~/git/biocode/sandbox/jorvis/reference_coverage_plot.py -i assembly.all.per_cov.txt,assembly.longest.per_cov.txt -l All,Longest -t "Repeat filtered Oases + Trinity, post-TGICL"

"""

import argparse
import os
import plotly.plotly as py
import plotly.graph_objs as go

def main():
    parser = argparse.ArgumentParser( description='Generates a graphic showing how well reference transcripts are covered by a transcript assembly')

    ## output file to be written
    parser.add_argument('-i', '--input_files', type=str, required=True, help='Comma-separated list of cov files to be plotted' )
    parser.add_argument('-l', '--labels', type=str, required=True, help='Labels for each cov file passed' )
    parser.add_argument('-t', '--title', type=str, required=False, default='Transcript coverage', help='Title for the plot' )
    args = parser.parse_args()

    cov_files = args.input_files.split(",")
    labels = args.labels.split(",")
    colors = ['rgb(49,130,189)',    #blue
              'rgb(204,204,204)',   #light grey
              'rgb(50, 171, 96)',   #green
              'rgb(222,45,38)',     #red
              'rgb(142, 124, 195)', #purple
              'rgb(100,100,100)' ]

    #print("Got {0} coverage files".format(len(cov_files)))
    #print("Got {0} labels".format(len(labels)))

    if len(labels) > len(colors):
        raise Exception("Sorry, this many datasets is not yet supported (only because not enough colors were defined in code.)")

    # This stores the positions of the labels
    label_position = dict()
    
    traces = []
    file_idx = 0
    for file in cov_files:
        xvals = []
        yvals = []
        
        for line in open(file):
            cols = line.rstrip().split("\t")
            xvals.append(cols[0])
            yvals.append(float(cols[1]))
        
        trace = go.Bar(
            x=xvals,
            y=yvals,
            name=labels[file_idx],
            marker=dict(
                color=colors[file_idx]
            )
        )

        traces.append(trace)
        file_idx += 1

    layout = go.Layout(
        title=args.title,
        xaxis=dict(
            # set x-axis' labels direction at 45 degree angle
            tickangle=-65
        ),
        yaxis=dict(
            title='Percent coverage',
            titlefont=dict(
                size=16,
                color='rgb(107, 107, 107)'
            ),
            tickfont=dict(
                size=16,
                color='rgb(107, 107, 107)'
            )
        ),
        legend=dict(
            #x=0,
            #y=1.2,
            bgcolor='rgba(255, 255, 255, 0)',
            bordercolor='rgba(255, 255, 255, 0)',
            font=dict(
                size=20,
                color='#000'
            )
        ),
        barmode='group',
        bargap=0.15,
        bargroupgap=0.1
    )
    fig = go.Figure(data=traces, layout=layout)
    plot_url = py.plot(fig, filename='angled-text-bar')


if __name__ == '__main__':
    main()







