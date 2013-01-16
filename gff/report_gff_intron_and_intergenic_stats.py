#!/usr/bin/env python3.2

import argparse
import biothings
import biocodegff

def main():
    '''
    This script reports statistics on the areas of a genome where features aren't - introns and
    intergenic space.  Pass a valid GFF3 file and get a report like this:

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
    '''
    parser = argparse.ArgumentParser( description='Reports statistics of reference gene coverage and extension by aligned RNA-seq transcript data.')

    ## output file to be written
    parser.add_argument('-i', '--input_gff3', type=str, required=True, help='GFF3 file of a reference annotation' )
    args = parser.parse_args()

    (assemblies, features) = biocodegff.get_gff3_features( args.input_gff3 )

    ## things to keep stats on and report
    total_molecule_count = len(assemblies)
    total_gene_count = 0
    
    ## this number is NOT just the total genes N - 1, since there can be multiple molecules
    #   genes can overlap, etc.
    total_intergenic_space_count = 0
    
    total_intergenic_space_residues = 0
    intergenic_distances = list()

    total_intron_count = 0
    total_intron_residues = 0
    intron_sizes = list()

    ############################
    ## Calculation section
    ############################
    
    for asm_id in assemblies:
        assembly = assemblies[asm_id]
        genes = assembly.genes()
        total_gene_count += len(genes)
        last_gene_loc = None

        for gene in sorted(genes):
            gene_loc = gene.location_on(assembly)

            if last_gene_loc is not None:
                ## skip this gene if it overlaps the previous
                if gene_loc.fmin < last_gene_loc.fmax:
                    if gene_loc.fmax > last_gene_loc.fmax:
                        last_gene_loc = gene_loc

                else:
                    total_intergenic_space_count += 1
                    intergenic_distance = gene_loc.fmin - last_gene_loc.fmax
                    total_intergenic_space_residues += intergenic_distance
                    intergenic_distances.append(intergenic_distance)
                    
            for mRNA in gene.mRNAs():
                introns = mRNA.introns( on=assembly )

                for intron in sorted(introns):
                    total_intron_count += 1
                    intron_loc = intron.location_on(assembly)
                    intron_size = intron_loc.fmax - intron_loc.fmin
                    intron_sizes.append(intron_size)
                    total_intron_residues += intron_size
                
            last_gene_loc = gene_loc


    avg_intergenic_space_dist = total_intergenic_space_residues / total_intergenic_space_count
    intergenic_distances = sorted(intergenic_distances)
    median_int_space_dist = intergenic_distances[ int(len(intergenic_distances)/2) ]

    avg_intron_size = total_intron_residues / total_intron_count
    intron_sizes = sorted(intron_sizes)
    median_intron_size = intron_sizes[int(len(intron_sizes)/2)]
            
    ############################
    ## Reporting section
    ############################

    print("\nMolecule count: {0}\n".format(total_molecule_count))
    print("Gene count: {0}".format(total_gene_count) )
    print("Intergenic space count: {0}".format(total_intergenic_space_count) )

    print("Average intergenic space distance: {0:.1f} bp".format(avg_intergenic_space_dist) )
    print("Median intergenic space distance: {0} bp".format(median_int_space_dist) )
    print("Minimum intergenic space distance: {0} bp".format(intergenic_distances[0]) )
    print("Maximum intergenic space distance: {0} bp\n".format(intergenic_distances[-1]) )
 
    print("Intron count: {0}".format(total_intron_count) )
    print("Intron space count: {0}".format(total_intron_residues) )

    print("Average intron size: {0:.1f} bp".format(avg_intron_size) )
    print("Median intron size: {0} bp".format(median_intron_size) )
    print("Minimum intron size: {0} bp".format(intron_sizes[0]) )
    print("Maximum intron size: {0} bp\n".format(intron_sizes[-1]) )
    
    

if __name__ == '__main__':
    main()








