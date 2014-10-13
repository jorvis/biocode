#!/usr/bin/env python3

'''
Author:  Joshua Orvis

If the GFF has gene.locus_tag annotation that will be used for the FASTA headers, else it will
be the ID attributes.

'''

import argparse
import os
import biocodegff
import biocodeutils
import sys

def main():
    parser = argparse.ArgumentParser( description='')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    parser.add_argument('-f', '--fasta', type=str, required=False, help='Required if you don\'t have GFF3 with embedded FASTA')
    args = parser.parse_args()

    (assemblies, features) = biocodegff.get_gff3_features( args.input_file )

    if args.fasta is not None:
        seqs = biocodeutils.fasta_dict_from_file( args.fasta )
        for seq_id in seqs:
            if seq_id in assemblies:
                assemblies[seq_id].residues = seqs[seq_id]['s']
                assemblies[seq_id].length = len(assemblies[seq_id].residues)

    ## output will either be a file or STDOUT
    ofh = sys.stdout
    if args.output_file is not None:
        ofh = open(args.output_file, 'wt')

    for assembly_id in assemblies:
        assembly = assemblies[assembly_id]
        
        for gene in assembly.genes():
            # for debugging only
            #if gene.id != '668964A84D72393D98D75CDE21350D8E':
                #continue

            if gene.locus_tag is None:
                gene_label = gene.id
            else:
                gene_label = gene.locus_tag
            
            gene_seq = gene.get_residues().upper()
            gene_loc = gene.location_on(assembly)

            ## we have to do this here because of the coordinates
            if gene_loc.strand == -1:
                gene_seq = "".join(reversed(gene_seq))
            
            #print("INFO: Processing gene with length {0} at {1}-{2}".format(len(gene_seq), gene_loc.fmin, gene_loc.fmax))

            if len(gene.mRNAs()) > 1:
                raise Exception("ERROR: script doesn't currently support multi-isoform genes, but found one: {0}".format(gene.id))
            
            for mRNA in gene.mRNAs():
                introns = mRNA.introns( on=assembly )
                
                for intron in introns:
                    intron_loc = intron.location_on(assembly)

                    # this helps us get where the intron is on the gene
                    offset = gene_loc.fmin
                    
                    #print("INFO:\tfound intron at {0}-{1}".format(intron_loc.fmin, intron_loc.fmax))
                    lower_mid = gene_seq[intron_loc.fmin - offset:intron_loc.fmax - offset].lower()
                    #print("INFO:\tlower-casing offset adjusted coordinates: {0}-{1}".format(intron_loc.fmin - offset, intron_loc.fmax - offset))
                    #print("INFO:\tgenerating lower case seq of length: {0}\n".format(len(lower_mid)) )
                    gene_seq = gene_seq[0:intron_loc.fmin - offset] + lower_mid + gene_seq[intron_loc.fmax - offset:]


            ## make sure to switch it back
            if gene_loc.strand == -1:
                gene_seq = "".join(reversed(gene_seq))
                    
            #print("INFO: Got gene with length {0} after modification".format(len(gene_seq)))
            ofh.write(">{0}\n{1}\n".format(gene_label, biocodeutils.wrapped_fasta(gene_seq)))



if __name__ == '__main__':
    main()






