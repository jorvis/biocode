#!/usr/bin/env python3

"""

Initially written to plot results from Priti's pipeline, where the goal was to determine
how well an assembled transcriptome (Trinity) produced contigs which covered a set of
reference transcripts.  Steps in her pipeline:

export BASE=A8_muscle.trinity.standard
export REF=toi_20200330.fasta
export TITLE=PASA

formatdb -p F -i $BASE.fasta
blastall -p blastn -i $REF -m 9 -e 1 -d $BASE.fasta -o $BASE.blast.m9
/home/jorvis/git/biocode/blast/calculate_query_coverage_by_blast.py -f $REF -b $BASE.blast.m9 -o $BASE
sort $BASE.cov.all.perc.txt > $BASE.cov.all.perc.sorted.txt
sort $BASE.cov.longest.perc.txt > $BASE.cov.longest.perc.sorted.txt

# then to actually plot these:
/usr/bin/python3 ~/git/biocode/sandbox/jorvis/reference_coverage_plot.py -i $BASE.cov.longest.perc.sorted.txt,$BASE.cov.all.perc.sorted.txt -l "Longest,All" -t "${TITLE} transcript coverage" -rf $REF -qf $BASE.fasta -o $BASE.both.png

/usr/bin/python3 ~/git/biocode/sandbox/jorvis/reference_coverage_plot.py -i $BASE.cov.longest.perc.sorted.txt -l "Longest" -t "${TITLE} - Longest transcript coverage" -rf $REF -qf $BASE.fasta -o $BASE.longest.png

/usr/bin/python3 ~/git/biocode/sandbox/jorvis/reference_coverage_plot.py -i $BASE.cov.all.perc.sorted.txt -l "All" -t "${TITLE} - All transcript coverage" -rf $REF -qf $BASE.fasta -o $BASE.all.png

--stacked option:
Rather than just show coverage with the max Y value at 100, use of this option creates a
stacked bar chart which gives information about the contig matching the reference and their
relative lengths.  This allows us to see if a contig covers a reference transcript completely
BUT is also far longer than the reference.

If you're having errors in the image generation, you have to follow the complete cluster that is setting
up local static image export in plotly:

https://plotly.com/python/static-image-export/

For me this involved:

sudo /usr/bin/pip3 install plotly==4.5.4
sudo /usr/bin/pip3 install psutil requestsipywidgets

And renaming this to 'orca' and putting it in PATH:

wget https://github.com/plotly/orca/releases/download/v1.3.1/orca-1.3.1.AppImage

"""

import argparse
import biocode.utils
import os
import plotly
import plotly.offline as py
import plotly.graph_objs as go

plotly.io.orca.config.executable = '/opt/bin/orca'

def main():
    parser = argparse.ArgumentParser( description='Generates a graphic showing how well reference transcripts are covered by a transcript assembly')

    ## output file to be written
    parser.add_argument('-i', '--input_files', type=str, required=True, help='Comma-separated list of cov files to be plotted' )
    parser.add_argument('-l', '--labels', type=str, required=True, help='Labels for each cov file passed' )
    parser.add_argument('-t', '--title', type=str, required=False, default='Transcript coverage', help='Title for the plot' )
    parser.add_argument('-s', '--stacked', dest='stacked', action='store_true')
    parser.set_defaults(stacked=False)
    parser.add_argument('-rf', '--ref_fasta', required=False, help='Only needed if passing --stacked')
    parser.add_argument('-qf', '--qry_fasta', required=False, help='Only needed if passing --stacked')
    parser.add_argument('-mb', '--margin_bottom', type=int, required=False, default=120, help='Size of the bottom margin, in case X labels are being cut off')
    parser.add_argument('-o', '--output_image', type=str, required=False, help='Name for PNG file to be created. If not passed, will post to plotly site' )
    args = parser.parse_args()

    cov_files = args.input_files.split(",")
    labels = args.labels.split(",")
    colors = ['rgb(49,130,189)',    #blue
              'rgb(204,204,204)',   #light grey
              'rgb(50, 171, 96)',   #green
              'rgb(222,45,38)',     #red
              'rgb(142, 124, 195)', #purple
              'rgb(100,100,100)',   #darker grey
              'rgb(255,255,61)',    #yellow
              'rgb(255,169,58)'     #orange
             ]

    #print("Got {0} coverage files".format(len(cov_files)))
    #print("Got {0} labels".format(len(labels)))

    if len(labels) > len(colors):
        raise Exception("Sorry, this many datasets is not yet supported (only because not enough colors were defined in code.)")

    # This stores the positions of the labels
    label_position = dict()

    if len(cov_files) > 1 and args.stacked == True:
        raise Exception("Use of the --stacked option requires a single input file")

    # Only used if doing a single-file stacked bar chart
    if args.stacked == True:
        ref_sizes = biocode.utils.fasta_sizes_from_file(args.ref_fasta)
        qry_sizes = biocode.utils.fasta_sizes_from_file(args.qry_fasta)
    
    traces = []
    file_idx = 0
    for file in cov_files:
        xvals = []
        yvals = []
        stacked_yvals = []
        
        for line in open(file):
            cols = line.rstrip().split("\t")
            ref_id = cols[0]
            xvals.append(ref_id)
            yvals.append(float(cols[1]))

            if args.stacked == True:
                qry_id = cols[2]
                if qry_sizes[qry_id] > ref_sizes[ref_id]:
                    # what percentage larger than the reference is the query?
                    rel_perc = (qry_sizes[qry_id] / ref_sizes[ref_id]) * 100

                    # for the stacked bars we have to subtract the current yval cov, since this adds to it
                    rel_perc_adj = rel_perc - float(cols[1])
                    stacked_yvals.append(rel_perc_adj)
                else:
                    stacked_yvals.append(0)
        
        trace = go.Bar(
            x=xvals,
            y=yvals,
            name=labels[file_idx],
            marker=dict(
                color=colors[file_idx]
            )
        )

        traces.append(trace)

        if args.stacked == True:
            trace2 = go.Bar(
                x=xvals,
                y=stacked_yvals,
                name=labels[file_idx],
                marker=dict(
                    color='rgb(200,200,200)'
                )
            )

            traces.append(trace2)

        file_idx += 1

    if args.stacked == True:
        barmode = 'stack'
    else:
        barmode = 'group'

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
        barmode=barmode,
        bargap=0.15,
        bargroupgap=0.1,
        width=1500,
        height=800,
        margin = go.Margin(b = args.margin_bottom, pad=5)
    )
    fig = go.FigureWidget(data=traces, layout=layout)

    if args.output_image is None:
        plot_url = py.plot(fig, filename='angled-text-bar')
        print("Navigate to {0} for your image".format(plot_url))
    else:
        #py.image.save_as(fig, filename=args.output_image)
        fig.write_image(args.output_image)
        print("Output written to file: {0}".format(args.output_image))


if __name__ == '__main__':
    main()







