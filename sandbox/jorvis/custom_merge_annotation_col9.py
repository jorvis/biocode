#!/usr/bin/env python3

"""

The original R1 annotation mapped to the 4 chrm assembly is here:
/usr/local/projects/b_microti/annotation/R1_4chrm_Genoscope_re-annotation/Bmicroti_new_karyotype_annotation.phasecorrected.gff3

The gff3 of the new annotation, which we want to modify, is here:
/usr/local/projects/b_microti/annotation/Final_annotation_20141120/Bmi_sorted_with_functional_locustags_fixed.gff3

In summary, the goal of the script is to add in the last (free for all) field of the new annotation the information about the genes in the original annotation overlapped by each gene in the new (final) annotation.

The qualifier is “overlaps_old_locusTagID= <old-locus-tag-id>”.
We do not differentiate between strand on which the original genes are located.
But we may also include strand, maybe as follows
overlaps_old_locusTagID= <old-locus-tag-id> in strand=<strand>; <same info for another gene here if two are overlapped>
"""

import argparse

from biocode import gff


def main():
    parser = argparse.ArgumentParser( description='Script for reporting of possible polycistronic genes transcripts based on a reference annotation and RNA-seq transcript assemblies')

    ## output file to be written
    parser.add_argument('-r', '--reference_file', type=str, required=True, help='GFF3 file of a reference annotation' )
    parser.add_argument('-q', '--query_file', type=str, required=True, help='GFF3 file with alternative annotation (such as an RNA-seq assemby)' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    (ref_assemblies, ref_feats) = gff.get_gff3_features(args.reference_file)
    (qry_assemblies, qry_genes) = gff.get_gff3_features(args.query_file)

    for assembly_id in ref_assemblies:
        # we expect to find this assembly ID in the qry set too
        if assembly_id not in qry_assemblies:
            print("WARN: expected to find assembly_id {0} in both reference and query sets".format(assembly_id))
            continue
        
        for ref_gene in ref_assemblies[assembly_id].genes():
            overlaps = list()
            polypeptides = ref_gene.polypeptides()

            if len(polypeptides) == 0:
                print("WARN: skipped gene {0} because it has no polypeptides".format(ref_gene.id))
                continue
                
            ref_annot = ref_gene.polypeptides()[0].annotation
            
            for qry_gene in qry_assemblies[assembly_id].genes():
                overlap = ref_gene.overlaps_with(qry_gene)
                
                if overlap:
                    #print("DEBUG: {0} and {1} appear to overlap".format(ref_gene.id, qry_gene.id) )
                    overlaps.append(overlap)
                    # add a dbxref to the gene
                    ref_annot.add_dbxref("overlaps_old_locusTagID:{0}".format(qry_gene.id))

            if len(overlaps) > 0:
                print("INFO: ref_gene {0} had {1} overlaps".format(ref_gene.id, len(overlaps)))
    
    gff.print_gff3_from_assemblies(assemblies=ref_assemblies, ofh=open(args.output_file, 'w'))

    

if __name__ == '__main__':
    main()




