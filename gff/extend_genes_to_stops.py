#!/usr/bin/env python3

'''
This script reads a GFF3 file and FASTA file (or FASTA embedded in the GFF) and
attempts to extend the termini of genes which do not end with a stop codon.  It
is splice-aware, so this works for prokaryotic or eukaryotic models.

Example input gene model:
gnl|WGS:AAGK|SeqID|gb|AAGK01000001      .       gene    128156  128682  .       -       .       ID=210D90E7BF533A9A08EF22D9F1F1CED9_1_2;locus_tag=TpMuguga_01g00060
gnl|WGS:AAGK|SeqID|gb|AAGK01000001      .       mRNA    128156  128682  .       -       .       ID=0D5133EC9A534184443BC6C9D05E565B;Parent=210D90E7BF533A9A08EF22D9F1F1CED9_1_2;locus_tag=TpMuguga_01g00060
gnl|WGS:AAGK|SeqID|gb|AAGK01000001      .       CDS     128349  128616  .       -       0       ID=5BBE8CD7940726C5808483FC5EC65D28;Parent=0D5133EC9A534184443BC6C9D05E565B
gnl|WGS:AAGK|SeqID|gb|AAGK01000001      .       CDS     128156  128278  .       -       0       ID=5BBE8CD7940726C5808483FC5EC65D28;Parent=0D5133EC9A534184443BC6C9D05E565B
gnl|WGS:AAGK|SeqID|gb|AAGK01000001      .       exon    128156  128278  .       -       .       ID=3E36768E591D34F98C62893C736848FB;Parent=0D5133EC9A534184443BC6C9D05E565B
gnl|WGS:AAGK|SeqID|gb|AAGK01000001      .       exon    128349  128682  .       -       .       ID=3A9F94FDF76013E7216D47419026FDF1;Parent=0D5133EC9A534184443BC6C9D05E565B
gnl|WGS:AAGK|SeqID|gb|AAGK01000001      .       polypeptide     128156  128682  .       -       .       ID=0D5133EC9A534184443BC6C9D05E565B_p;Parent=0D5133EC9A534184443BC6C9D05E565B;product_name=hypothetical protein

It provides output of this process like this:

   gene:210D90E7BF533A9A08EF22D9F1F1CED9_1_2, mRNA: 0D5133EC9A534184443BC6C9D05E565B is missing a stop
           mRNA:128155-128682, CDS end: 128156
           Extending.GGT(128153).GAG(128150).TAA(128147) Found a stop

Then a summary at the end like:

   Total mRNAs found:4079
   mRNAs initially with terminal stops: 3953
   mRNAs successfully extended: 104

You can limit how far extensions will go using the --extension_limit argument

Follow the GFF3 specification!

Author:  Joshua Orvis
'''

import argparse
from biocode import utils, gff

def main():
    parser = argparse.ArgumentParser( description='Extends GFF gene models to the first in-frame stop')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-g', '--genome_fasta', type=str, required=False, help='Optional.  You must specify this unless the FASTA sequences for the molecules are embedded in the GFF')
    parser.add_argument('-o', '--output_gff', type=str, required=False, help='Optional.  Writes an output GFF3 file with CDS (and containing features) extended to nearest stop')
    parser.add_argument('-el', '--extension_limit', type=int, required=False, default=100, help='Optional.  Limits how far an extension will happen looking for an in-frame stop codon')
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_file)

    # deal with the FASTA file if the user passed one
    if args.genome_fasta is not None:
        utils.add_assembly_fasta(assemblies, args.genome_fasta)

    total_mRNAs = 0
    mRNAs_with_terminal_stops = 0
    stop_codons = ['TAG', 'TAA', 'TGA']
    mRNAs_corrected = 0

    for assembly_id in sorted(assemblies):
        print("Assembly {0} has length {1}".format(assembly_id, assemblies[assembly_id].length))
        for gene in sorted(assemblies[assembly_id].genes()):
            for mRNA in gene.mRNAs():
                coding_seq = mRNA.get_CDS_residues()
                total_mRNAs += 1
                translation = utils.translate(coding_seq)

                if translation.endswith('*'):
                    mRNAs_with_terminal_stops += 1
                else:
                    print("gene:{1}, mRNA: {0} is missing a stop".format(mRNA.id, gene.id))
                    print("\tCDS: {0}".format(coding_seq))
                    print("\tcoding sequence ends with {0}, last three a.a.: {1}".format(coding_seq[-3:], translation[-3:]))
                    mRNA_loc = mRNA.location_on(assemblies[assembly_id])
                    
                    CDSs = sorted(mRNA.CDSs())
                    CDS_frame_overhang = len(coding_seq) % 3
                    print("\tCDS frame overhang: {0}".format(CDS_frame_overhang))
                    codon_step_size = 3

                    if mRNA_loc.strand == 1:
                        # get the in-frame end coordinate of the last CDS position
                        CDS_pos = CDSs[-1].location_on(assemblies[assembly_id]).fmax - CDS_frame_overhang
                        mRNA_limit = mRNA_loc.fmax + args.extension_limit
                    else:
                        # get the in-frame end coordinate of the last CDS position
                        CDS_pos = CDSs[0].location_on(assemblies[assembly_id]).fmin + CDS_frame_overhang
                        mRNA_limit = mRNA_loc.fmin - args.extension_limit
                        codon_step_size = -3

                    print("\tmRNA:{0}-{1} ({3}), CDS end: {2}.  Extending ... \n\t".format(mRNA_loc.fmin, mRNA_loc.fmax, CDS_pos, mRNA_loc.strand), end='')

                    new_stop_found = False

                    # We have to step backwards to start if on the reverse strand
                    CDS_pos += codon_step_size

                    while True:
                        if (mRNA_loc.strand == 1 and CDS_pos > mRNA_limit) or (mRNA_loc.strand == -1 and CDS_pos < mRNA_limit):
                            print("  Reached the mRNA limit")
                            break
                        elif CDS_pos < 1:
                            print("  Reached beginning of the molecule")
                            break
                        else:
                            next_codon = assemblies[assembly_id].residues[CDS_pos:CDS_pos + 3]
                            
                            if mRNA_loc.strand == -1:
                                next_codon = utils.reverse_complement(next_codon)
                                print(".{0}({1}-{2})".format(next_codon, CDS_pos, CDS_pos - 3), end='')
                            else:
                                print(".{0}({1}-{2})".format(next_codon, CDS_pos - 3, CDS_pos), end='')
                        
                            if next_codon in stop_codons:
                                if mRNA_loc.strand == 1:
                                    mRNA.extend_stop(on=assemblies[assembly_id], to=(CDS_pos + 3))
                                    print(" Found a stop, extending to: {0} ({1})".format(CDS_pos + 3, mRNA_loc.strand))
                                else:
                                    mRNA.extend_stop(on=assemblies[assembly_id], to=CDS_pos)
                                    print(" Found a stop, extending to: {0} ({1})".format(CDS_pos, mRNA_loc.strand))

                                new_stop_found = True
                                break

                        CDS_pos += codon_step_size

                    if new_stop_found == True:
                        print("\tCDS_pos: UPDATE: {0}".format(CDS_pos))
                        mRNAs_corrected += 1
                    else:
                        print("\tCDS_pos:   SAME: {0}".format(CDS_pos))


    print("\nTotal mRNAs found:{0}".format(total_mRNAs))
    print("mRNAs initially with terminal stops: {0}".format(mRNAs_with_terminal_stops))
    print("mRNAs successfully extended: {0}".format(mRNAs_corrected))

    ofh = open(args.output_gff, 'wt')
    gff.print_gff3_from_assemblies(assemblies=assemblies, ofh=ofh)


if __name__ == '__main__':
    main()






