#!/usr/bin/env python3

'''
The phase column of the GFF3 specification is often ignored.  It states:

    For features of type "CDS", the phase indicates where the feature begins with
    reference to the reading frame. The phase is one of the integers 0, 1, or 2,
    indicating the number of bases that should be removed from the beginning of
    this feature to reach the first base of the next codon. In other words, a phase
    of "0" indicates that the next codon begins at the first base of the region
    described by the current line, a phase of "1" indicates that the next codon
    begins at the second base of this region, and a phase of "2" indicates that
    the codon begins at the third base of this region. This is NOT to be confused
    with the frame, which is simply start modulo 3.

    For forward strand features, phase is counted from the start field. For reverse
    strand features, phase is counted from the end field.

    The phase is REQUIRED for all CDS features.

This script reads a GFF3 file and FASTA file (or FASTA embedded in the GFF) and
checks the coding phase of each CDS to report in-frame stops.  If one is found, the
other possible values of phase are tested.  If all contain stops, a warning issued.
If another phase value besides the current one has no stops, and the current one does,
it is corrected.

Warning: Depending on what you have in your GFF3 this script could be pretty
destructive.  It is considering only gene models and their children, so these features
are kept: gene, mRNA, CDS, exon.  Any other feature types and comments will be omitted
from the output file in the current implementation.

Follow the GFF3 specification!

Author:  Joshua Orvis
'''

import argparse

from biocode import utils, gff


def main():
    parser = argparse.ArgumentParser( description='Checks the CDS features against a genome sequence to report/correct phase columns.')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input GFF3' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-g', '--genome_fasta', type=str, required=False, help='Optional.  You must specify this unless the FASTA sequences for the molecules are embedded in the GFF')
    parser.add_argument('-s', '--source', type=str, required=False, default='.', help='Optional.  Sets the value for column 2 in all rows.  Default = .' )
    args = parser.parse_args()

    (assemblies, features) = gff.get_gff3_features(args.input_file)

    fout = open(args.output_file, mode='wt', encoding='utf-8')
    fout.write("##gff-version 3\n")

    # deal with the FASTA file if the user passed one
    if args.genome_fasta is not None:
        process_assembly_fasta(assemblies, args.genome_fasta)

    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            for mRNA in gene.mRNAs():
                for CDS in mRNA.CDSs():
                    check_and_update_phase(CDS)

            gene.print_as(fh=fout, source=args.source, format='gff3')

    fasta_header_written = False

    for assembly_id in assemblies:
        if assemblies[assembly_id].length > 0:
            if fasta_header_written is False:
                fout.write("##FASTA\n")
                fasta_header_written = True

            fout.write(">{0}\n".format(assemblies[assembly_id].id) )
            fout.write("{0}\n".format(utils.wrapped_fasta(assemblies[assembly_id].residues)))

def check_and_update_phase(CDS):
    loc = CDS.location()
    CDS.get_residues()

    best_phase = None
    orig_phase_stop_count = None
    best_phase_stop_count = None

    for phase in [ 0, 1, 2 ]:
        protein_seq = utils.translate(CDS.residues[phase:]).rstrip('*')
        stop_count = protein_seq.count('*')

        if phase == loc.phase:
            orig_phase_stop_count = stop_count

        if best_phase is None or stop_count < best_phase_stop_count:
            best_phase = phase
            best_phase_stop_count = stop_count
            continue
        
    if best_phase != loc.phase:
        print("INFO: CDS {0} at coordinate:{1}, phase:{2} had {3} stops.  Updating to phase:{4} which had {5}".format( \
               CDS.id, loc.fmin, loc.phase, orig_phase_stop_count, best_phase, best_phase_stop_count) )
        loc.phase = best_phase

    
    


def process_assembly_fasta(mols, fasta_file):
    fasta_seqs = utils.fasta_dict_from_file(fasta_file)

    for mol_id in mols:
        # check if the FASTA file provides sequence for this
        if mol_id in fasta_seqs:
            mol = mols[mol_id]
            mol.residues = fasta_seqs[mol_id]['s']
            mol.length   = len(mol.residues)
        

if __name__ == '__main__':
    main()






