#!/usr/bin/env python3

import argparse
import os
import re
import biothings
import biocodegff
import biocodeutils

"""
C. hominis
The goal is to annotate the best assembly (TU502) and then map annotation to  UK01 strain.

IV. Train Augustus using Cufflinks transcripts as models
   A. First, filter Cufflinks transcripts
      1. Identify all cufflinks transcripts that fully contain one and only one published C. hominis gene model:
         a. Must contain start & stop codon coordinates of gene in published annotation;
         b. Must not overlap any other C. hominis gene model;
         c. Must not overlap any other transcript;
         d. Must have AAT hit, where reciprocal AAT alignment coverage of annotated C. hominis
            gene model with identity is 98% 

   B. Augustus training & optimization
      1. Randomly select ~300 models for training & 150 for evaluation, selecting so that the exon count
         profile matches that in IV.1.d above.
      2. Train Augustus to generate a parameter file;
      3. Do optimization of parameter file to increase Sn/Sp.

Test command:

Just run it.  For documentation, everything is embedded.

  $ cd /usr/local/projects/cryptosporidium
  $ ~/git/biocode/sandbox/jorvis/custom.c_hominis_training_set.py


NOTES:
We have new assemblies, so the reference ones had to be mapped using gmap/gsnap:

"""

def main():
    parser = argparse.ArgumentParser( description='Custom script to select Cryptosporidium hominis training genes')

    #parser.add_argument('-o', '--output_id_list', type=str, required=False, help='List of IDs from organism1 that passed' )
    args = parser.parse_args()

    # reference CDS, mapped to our genome using GMAP with -f gff3_gene
    #ref_gff_file = '/usr/local/projects/cryptosporidium/structural_prediction/gmap/ref_vs_our_c_hominis_TU502/gmap.gff'
    ref_gff_file = '/usr/local/projects/cryptosporidium/structural_prediction/cegma/c_hominis_TU502/output.cegma.gff3'
    gmap_pct_id_cutoff = 97
    gmap_cov_cutoff = 100

    # output of cufflinks assembly, translated to GFF3 with my script: biocode/gff/convert_cufflinks_gtf_to_gff3.py
    cufflinks_model_gff_file = '/usr/local/projects/cryptosporidium/rna_seq/c_hominis_TU502/tophat_cufflinks/cufflinks_transcript_assembly.models.gff3';

    print("INFO: Parsing C. hominis TU502 reference annotation (mapped with gmap):")
    (ref_assemblies, ref_features) = biocodegff.get_gff3_features( ref_gff_file )
    print("INFO: Parsed {0} features on {1} assemblies".format(len(ref_features), len(ref_assemblies)) )

    print("INFO: Parsing Cufflinks assemblies")
    (cuff_assemblies, cuff_features) = biocodegff.get_gff3_features( cufflinks_model_gff_file )
    print("INFO: Parsed {0} features on {1} assemblies".format(len(cuff_features), len(cuff_assemblies)) )

    cuffs_matching_refs = list()
    
    ## A.1.a: Must contain start/stop coordinates of gene in published annotation
    for asm_id in cuff_assemblies:
        cuff_assembly = cuff_assemblies[asm_id]

        # there may not have been any mapped ref genes on this assembly
        if asm_id not in ref_assemblies:
            continue

        ref_assembly = ref_assemblies[asm_id]

        for cuff_gene in cuff_assembly.genes():
            cuff_gene_loc = cuff_gene.location_on(cuff_assembly)
            overlaps_found = 0

            for ref_gene in ref_assembly.genes():
                ref_gene_loc = ref_gene.location_on(ref_assembly)

                if cuff_gene_loc.fmin == ref_gene_loc.fmin and cuff_gene_loc.fmax == ref_gene_loc.fmax:
                    print("DEBUG: comparing cuff at {0}-{1} with ref at {2}-{3}".format( \
                        cuff_gene_loc.fmin, cuff_gene_loc.fmax, ref_gene_loc.fmin, ref_gene_loc.fmax) )
                
                #if cuff_gene.has_same_coordinates_as( thing=ref_gene ):
                #    print("DEBUG: found an overlap")
                #    overlaps_found += 1

            if overlaps_found == 1:
                cuffs_matching_refs.append(cuff_gene)

    print("INFO: There were {0} cufflinks transcripts with only one matching ref gene (exact min/max coordinates)\n".format(len(cuffs_matching_refs)) )

    


if __name__ == '__main__':
    main()







