#!/usr/bin/env python3

import argparse

from biocode import utils, gff


def main():
    '''
    This script reports statistics on the areas of a genome where features aren't - introns and
    intergenic space.  Pass a valid GFF3 file (along with FASTA data) and get a report like this:

    Molecule count: 9

    Gene count: 4171
    Intergenic space count: 4061
    Average intergenic space distance: 361.7 bp
    Median intergenic space distance: 245 bp
    Minimum intergenic space distance: 0 bp
    Maximum intergenic space distance: 6272 bp

    Intron count: 10533
    Intron space count: 989024
    Average intron size: 93.9 bp
    Median intron size: 63 bp
    Minimum intron size: 2 bp
    Maximum intron size: 1676 bp


    Optionally, you can pass the path to a PNG file to be created using the --histogram parameter,
    which will generate a size distribution histogram with two overlaying plots - one representing
    the distribution of intergenic region sizes and the other the intron lengths.  Because these
    can often have long tails, you can limit both the Y- and X-axes values with the --ylimit and
    --xlimit options, respectively.

    FASTA:
    If your FASTA isn't embedded at the end of your GFF3 file after a ##FASTA directive you'll need
    to specify the --fasta option in this script and pass it as a separate file.

    Definitions:
    Intergenic space was a little ambiguous to me as I started writing this.  Does one count the space from
    the beginning of the contig until the first gene, or only between them?  What about short contigs which
    have no annotated genes at all?  From the Sequence Ontology:

    SO:0000605: A region containing or overlapping no genes that is bounded on either side by a gene, or
    bounded by a gene and the end of the chromosome.

    To my reading, this includes contig ends but not gene-less contigs.  To that end, I include the
    former in intergenic space reporting but include the latter as a separate statistic.

    Author: Joshua Orvis (jorvis AT gmail)
    '''
    parser = argparse.ArgumentParser( description='Reports statistics of reference gene coverage and extension by aligned RNA-seq transcript data.')

    ## output file to be written
    parser.add_argument('-i', '--input_gff3', type=str, required=True, help='GFF3 file of a reference annotation' )
    parser.add_argument('-g', '--histogram', type=str, required=False, help='Optional path to a histogram of intron/intergenic space size distribution to be created (PNG)' )
    parser.add_argument('-x', '--xlimit', type=int, required=False, help='Use this if you want to limit the X-axis of the histogram (feature length)' )
    parser.add_argument('-y', '--ylimit', type=int, required=False, help='Use this if you want to limit the Y-axis of the histogram (feature count)' )
    parser.add_argument('-f', '--fasta', type=str, required=False, help='Required if you don\'t have GFF3 with embedded FASTA')
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_gff3)

    if args.fasta is not None:
        seqs = utils.fasta_dict_from_file(args.fasta)
        for seq_id in seqs:
            if seq_id in assemblies:
                assemblies[seq_id].residues = seqs[seq_id]['s']
                assemblies[seq_id].length = len(assemblies[seq_id].residues)

    ## things to keep stats on and report
    total_molecule_count = len(assemblies)
    total_gene_count = 0
    
    ## this number is NOT just the total genes N - 1, since there can be multiple molecules
    #   genes can overlap, etc.
    total_intergenic_space_count = 0
    
    total_intergenic_space_residues = 0
    intergenic_distances = list()

    total_contig_residues = 0
    empty_contig_residues = 0

    total_intron_count = 0
    total_intron_residues = 0
    intron_sizes = list()

    ############################
    ## Calculation section
    ############################

    for asm_id in assemblies:
        #print("DEBUG: processing assembly: {0}".format(asm_id))
        assembly = assemblies[asm_id]
        genes = sorted(assembly.genes())
        total_gene_count += len(genes)
        previous_gene_loc = None

        # we should have a length here
        if assembly.length is None or assembly.length == 0:
            raise Exception("ERROR: Detected assembly with undefined or 0 length: {0}".format(assembly.id))

        if total_gene_count == 0:
            empty_contig_residues += assembly.length
            continue

        total_contig_residues += assembly.length
        first_gene_loc = None
        last_gene_loc = None

        for gene in genes:
            gene_loc = gene.location_on(assembly)

            # if this is the first gene, track the number of bases from the start of the molecule here
            if first_gene_loc is None:
                total_intergenic_space_count += 1
                intergenic_distance = gene_loc.fmin
                total_intergenic_space_residues += intergenic_distance
                intergenic_distances.append(intergenic_distance)
                first_gene_loc = gene_loc

            if previous_gene_loc is not None:
                ## skip this gene if it overlaps the previous
                if gene_loc.fmin < previous_gene_loc.fmax:
                    if gene_loc.fmax > previous_gene_loc.fmax:
                        previous_gene_loc = gene_loc

                else:
                    total_intergenic_space_count += 1
                    intergenic_distance = gene_loc.fmin - previous_gene_loc.fmax
                    total_intergenic_space_residues += intergenic_distance
                    intergenic_distances.append(intergenic_distance)
                    
            for mRNA in gene.mRNAs():
                introns = mRNA.introns( on=assembly )

                for intron in sorted(introns):
                    total_intron_count += 1
                    intron_loc = intron.location_on(assembly)
                    intron_size = intron_loc.fmax - intron_loc.fmin

                    #if intron_size > 0:
                        #print("\tDEBUG: found mRNA:{0} intron {1}-{2} ({3} bp)".format(mRNA.id, intron_loc.fmin, intron_loc.fmax, intron_size))

                    if intron_size < 0:
                        print("\tWARN: Intron size ({1}) < 0 reported in gene {0}".format(gene.id, intron_size))
                    
                    intron_sizes.append(intron_size)
                    total_intron_residues += intron_size
                
            previous_gene_loc = gene_loc
            last_gene_loc = previous_gene_loc
        
        if last_gene_loc is not None:
            total_intergenic_space_count += 1
            intergenic_distance = assembly.length - last_gene_loc.fmax
            total_intergenic_space_residues += intergenic_distance
            intergenic_distances.append(intergenic_distance)

    if total_intergenic_space_count == 0:
        avg_intergenic_space_dist = None
        intergenic_distances = None
        median_int_space_dist = None
    else:
        avg_intergenic_space_dist = total_intergenic_space_residues / total_intergenic_space_count
        intergenic_distances = sorted(intergenic_distances)
        median_int_space_dist = intergenic_distances[ int(len(intergenic_distances)/2) ]

    avg_intron_size = total_intron_residues / total_intron_count
    intron_sizes = sorted(intron_sizes)
    median_intron_size = intron_sizes[int(len(intron_sizes)/2)]
            
    ############################
    ## Reporting section
    ############################

    print("\nMolecule count: {0}".format(total_molecule_count))
    print("Gene count: {0}".format(total_gene_count) )

    print("\nTotal molecule bases: {0} bp".format(total_contig_residues) )
    print("Empty molecule bases: {0} bp".format(empty_contig_residues) )

    if total_intergenic_space_count > 0:
        print("Intergenic space count: {0}".format(total_intergenic_space_count) )
        print("Average intergenic space distance: {0:.1f} bp".format(avg_intergenic_space_dist) )
        print("Median intergenic space distance: {0} bp".format(median_int_space_dist) )
        print("Minimum intergenic space distance: {0} bp".format(intergenic_distances[0]) )
        print("Maximum intergenic space distance: {0} bp\n".format(intergenic_distances[-1]) )
    else:
        print("There were no intergenic spaces found.  This might mean there were no molecules with at least 2 genes.")
 
    print("Intron count: {0}".format(total_intron_count) )
    print("Intron space count: {0}".format(total_intron_residues) )

    print("Average intron size: {0:.1f} bp".format(avg_intron_size) )
    print("Median intron size: {0} bp".format(median_intron_size) )
    print("Minimum intron size: {0} bp".format(intron_sizes[0]) )
    print("Maximum intron size: {0} bp\n".format(intron_sizes[-1]) )
    
    ############################
    ## Graphics section (optional)
    ############################
    if args.histogram is not None:
        import matplotlib.pyplot as plt

        plt.xlabel('length (bp)')
        plt.ylabel('count')
        plt.title('Distribution of intron size and intergenic distances')
        plt.hist(intergenic_distances, bins=50, histtype='stepfilled', color='b', label='Intergenic distances' )
        plt.hist(intron_sizes, bins=50, histtype='stepfilled', color='r', alpha=0.5, label='Intron sizes' )

        if args.xlimit is not None:
            plt.xlim([0, args.xlimit])
        
        if args.ylimit is not None:
            plt.ylim([0, args.ylimit])

        plt.legend(loc='best')
        plt.savefig(args.histogram)

if __name__ == '__main__':
    main()








