#!/usr/bin/env python3

'''
Author:  Joshua Orvis

If the GFF has gene.locus_tag annotation that will be used for the FASTA headers, else it will
be the ID attributes.

By default this writes the sequences of the mRNA features with introns included.  If you
pass --type=CDS it will instead do the CDS range.

'''

import argparse
import sys

from biocode import utils, gff


def main():
    parser = argparse.ArgumentParser( description='')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Path to an output file to be created' )
    parser.add_argument('-f', '--fasta', type=str, required=False, help='Required if you don\'t have GFF3 with embedded FASTA')
    parser.add_argument('-t', '--type', type=str, required=False, default='mRNA', choices=['mRNA', 'CDS'], help='Feature type to export (mRNA or CDS)')
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_file)

    # set this to None if you don't want the debug print statements
    #debugging_gene = 'D9AE6116893A0D5711D56C0F1E6CF58C'
    debugging_gene = None

    if args.fasta is not None:
        seqs = utils.fasta_dict_from_file(args.fasta)
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

            if debugging_gene is not None:
                debug_mode = True
                if gene.id != debugging_gene: continue
            else:
                debug_mode = False

            if gene.locus_tag is None:
                gene_label = gene.id
            else:
                gene_label = gene.locus_tag
            
            gene_seq = gene.get_residues().upper()
            gene_loc = gene.location_on(assembly)

            ## we have to do this here because of the coordinates
            if gene_loc.strand == -1:
                gene_seq = "".join(reversed(gene_seq))

            if debug_mode:
                print("INFO: Processing gene with length {0} at {1}-{2}".format(len(gene_seq), gene_loc.fmin, gene_loc.fmax))

            if len(gene.mRNAs()) > 1:
                #raise Exception("ERROR: script doesn't currently support multi-isoform genes, but found one: {0}".format(gene.id))
                print("ERROR: skipping gene {0} because it appears to have multiple isoforms (not currently supported)".format(gene.id))
                continue

            
            for mRNA in gene.mRNAs():
                introns = mRNA.introns( on=assembly )

                # this helps us get where the intron is on the gene
                offset = gene_loc.fmin
                
                for intron in introns:
                    intron_loc = intron.location_on(assembly)
                    lower_mid = gene_seq[intron_loc.fmin - offset:intron_loc.fmax - offset].lower()
                    gene_seq = gene_seq[0:intron_loc.fmin - offset] + lower_mid + gene_seq[intron_loc.fmax - offset:]

                    if debug_mode:
                        print("INFO:\tfound intron at {0}-{1}".format(intron_loc.fmin, intron_loc.fmax))
                        print("INFO:\tlower-casing offset adjusted coordinates: {0}-{1}".format(intron_loc.fmin - offset, intron_loc.fmax - offset))
                        print("INFO:\tgenerating lower case seq of length: {0}\n".format(len(lower_mid)) )

                if debug_mode:
                    print("INFO: seq length before CDS processing is: {0}".format(len(gene_seq)))

                ## do we need to trim down to the CDS range?
                if args.type == 'CDS':
                    CDSs = sorted(mRNA.CDSs())
                    CDS_min = CDSs[0].location_on(assembly).fmin
                    CDS_max = CDSs[-1].location_on(assembly).fmax

                    if debug_mode:
                        print("INFO: Calculated CDS range, with introns, should be: {0}-{1}={2}".format(CDS_max, CDS_min, CDS_max - CDS_min))

                    if gene_loc.fmin != CDS_min or gene_loc.fmax != CDS_max:
                        fmin_chomp = CDS_min - offset
                        fmax_chomp = gene_loc.fmax - CDS_max

                        if debug_mode:
                            print("gene:{0} coords:{1}-{2} ({3}), CDS coords: {4}-{5}".format(gene.id, gene_loc.fmin, \
                                                                                      gene_loc.fmax, gene_loc.strand, \
                                                                                      CDS_min, CDS_max \
                                                                                     ))

                            print("\tfmin_chomp:{0}, fmax_chomp:{1}".format(fmin_chomp, fmax_chomp))
                            print("\tpulling range: gene_seq[{0} : {1}]".format(fmin_chomp, len(gene_seq) - fmax_chomp))
                            
                        gene_seq = gene_seq[fmin_chomp : len(gene_seq) - fmax_chomp]

                        if debug_mode:
                            print("\tGene {0} CDS seq: {1}".format(gene.id, gene_seq))

            ## make sure to switch it back
            if gene_loc.strand == -1:
                gene_seq = "".join(reversed(gene_seq))
                    
            #print("INFO: Got gene with length {0} after modification".format(len(gene_seq)))
            ofh.write(">{0}\n{1}\n".format(gene_label, utils.wrapped_fasta(gene_seq)))



if __name__ == '__main__':
    main()






