#!/usr/bin/env python3

'''
In the original GFF annotation, the CDS always extend to the boundaries of the gene/mRNA
rather than the actual boundary of the coding sequence.  This script takes as input the
flawed GFF annotation as well as the ILRI GFF sequence (whose polypeptide sequences are
correct) and adjusts the CDS accordingly.

TEST CASES:

tp.assembly.567497685.1 GenBank gene    251791  253104  .       -       .       ID=TP03_0744
tp.assembly.567497685.1 GenBank mRNA    251791  253104  .       -       .       ID=TP03_0744.t01;Parent=TP03_0744
tp.assembly.567497685.1 GenBank exon    251791  251870  .       -       .       ID=TP03_0744.t01_exon-auto19995;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank exon    251998  252228  .       -       .       ID=TP03_0744.t01_exon-auto19996;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank exon    252288  252771  .       -       .       ID=TP03_0744.t01_exon-auto19997;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank exon    252839  252891  .       -       .       ID=TP03_0744.t01_exon-auto19998;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank exon    253075  253104  .       -       .       ID=TP03_0744.t01_exon-auto19999;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank CDS     251791  251870  .       -       .       ID=TP03_0744.t01_CDS-auto19995;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank CDS     251998  252228  .       -       .       ID=TP03_0744.t01_CDS-auto19996;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank CDS     252288  252771  .       -       .       ID=TP03_0744.t01_CDS-auto19997;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank CDS     252839  252891  .       -       .       ID=TP03_0744.t01_CDS-auto19998;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank CDS     253075  253104  .       -       .       ID=TP03_0744.t01_CDS-auto19999;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank polypeptide     251815  252844  .       -       .       ID=TP03_0744.p01;Parent=TP03_0744.t01;D

tp.assembly.567497685.1 GenBank gene    251791  253104  .       -       .       ID=TP03_0744
tp.assembly.567497685.1 GenBank mRNA    251791  253104  .       -       .       ID=TP03_0744.t01;Parent=TP03_0744
tp.assembly.567497685.1 GenBank CDS     251791  251870  .       -       0       ID=TP03_0744.t01_CDS-auto19995;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank CDS     251998  252228  .       -       0       ID=TP03_0744.t01_CDS-auto19996;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank CDS     252288  252771  .       -       0       ID=TP03_0744.t01_CDS-auto19997;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank CDS     252839  252891  .       -       0       ID=TP03_0744.t01_CDS-auto19998;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank CDS     253075  253104  .       -       0       ID=TP03_0744.t01_CDS-auto19999;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank exon    251791  251870  .       -       .       ID=TP03_0744.t01_exon-auto19995;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank exon    251998  252228  .       -       .       ID=TP03_0744.t01_exon-auto19996;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank exon    252288  252771  .       -       .       ID=TP03_0744.t01_exon-auto19997;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank exon    252839  252891  .       -       .       ID=TP03_0744.t01_exon-auto19998;Parent=TP03_0744.t01
tp.assembly.567497685.1 GenBank exon    253075  253104  .       -       .       ID=TP03_0744.t01_exon-auto19999;Parent=TP03_0744.t01

Author:  Joshua Orvis
'''

from biocode import gff, things


def main():
    flawed_gff_file = 'canonical.flawed.gff3'
    ilri_gff = 'Theileria-all-Theileria1_ourids.gff'
    source = 'GenBank'
    out_gff = 'canonical.corrected.gff3'
    
    fout = open(out_gff, mode='wt', encoding='utf-8')
    fout.write("##gff-version 3\n")

    (assemblies, features) = gff.get_gff3_features(flawed_gff_file)

    print("INFO: loaded {0} assemblies and {1} features".format(len(assemblies), len(features)))

    polypeptides = dict()

    for line in open(ilri_gff):
        cols = line.split("\t")

        if len(cols) != 9 or cols[2] != 'polypeptide':
            continue

        id = gff.column_9_value(cols[8], 'ID')
        parent = gff.column_9_value(cols[8], 'Parent')
        polypeptides[parent] = things.Polypeptide(id=id, parent=parent)
        polypeptides[parent].locate_on(target=assemblies[cols[0]], fmin=int(cols[3]) - 1, fmax=int(cols[4]), strand=cols[6])

    print("DEBUG: loaded {0} polypeptides from ILRI file".format(len(polypeptides)) )

    for assembly_id in assemblies:
        for gene in assemblies[assembly_id].genes():
            for mRNA in gene.mRNAs():
                if mRNA.id not in polypeptides:
                    print("DEBUG: {0} not found as a parent to any polypeptide".format(mRNA.id))
                else:
                    polypeptide = polypeptides[mRNA.id]

                # pull this outside of the iteration since iterating might delete some
                CDSs = mRNA.CDSs()
                    
                for CDS in CDSs:
                    keep = True
                    
                    if CDS < polypeptide:
                        mRNA.delete_CDS(CDS)
                    elif CDS <= polypeptide:
                        CDS.location().fmin = polypeptide.location().fmin
                    if CDS > polypeptide:
                        mRNA.delete_CDS(CDS)
                    elif CDS >= polypeptide:
                        CDS.location().fmax = polypeptide.location().fmax
                        #print("WARN: found a CDS {0}:{1}-{2} outside the range of the polypeptide {3}:{4}-{5}".format( \
                        #        CDS.id, CDS.locations[0].fmin, CDS.locations[0].fmax, \
                        #        polypeptide.id, polypeptide.locations[0].fmin, polypeptide.locations[0].fmax))                    

            gene.print_as(fh=fout, source=source, format='gff3')
                    


if __name__ == '__main__':
    main()






