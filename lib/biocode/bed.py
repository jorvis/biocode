import sys

import biocode.utils
import biocode.things

def print_bed_from_assemblies(assemblies=None, ofh=None):
    """
    Utility function to write a BED file from a list() of biothings.Assembly objects.

    No headers are used so the output is compatible with other tools like bedToBigBed

    References:

      BED format descriptions:
      - http://genome.ucsc.edu/FAQ/FAQformat#format1
      - http://bedtools.readthedocs.io/en/latest/content/general-usage.html

      bedToBigBed:
      - http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed

      Representing gene models in BED:
      - http://transvar.org/6111/gene_models.pdf

    """
    if type(ofh).__name__ != 'TextIOWrapper':
        ofh = sys.stdout #TextIOWrapper

    for assembly_id in assemblies:
        current_assembly = assemblies[assembly_id]

        for gene in sorted(assemblies[assembly_id].genes()):
            rnas_found = 0
            mRNAs = gene.mRNAs()
            
            for mRNA in mRNAs:
                mRNA_loc = mRNA.location_on(current_assembly)
                rnas_found += 1

                if rnas_found > 1:
                    gene.remove_mRNA(mRNA)
                    
                    print("INFO: splitting mRNA off gene {0}".format(gene.id))
                    new_gene = biocode.things.Gene(id="{0}_{1}".format(gene.id, rnas_found))
                    
                    new_gene.locate_on(target=current_assembly, fmin=mRNA_loc.fmin, fmax=mRNA_loc.fmax, strand=mRNA_loc.strand)
                    new_gene.add_RNA(mRNA)
                    print_biogene(gene=new_gene, fh=ofh, obo_dict=go_idx, lab_name=lab_name)

            if len(mRNAs) > 1:
                gene_loc = gene.location_on(current_assembly)
                mRNA_loc = mRNAs[0].location_on(current_assembly)
                gene_loc.fmin = mRNA_loc.fmin
                gene_loc.fmax = mRNA_loc.fmax
                gene_loc.strand = mRNA_loc.strand

            print_biogene(gene=gene, fh=ofh)


def print_biogene( gene=None, fh=None, on=None, obo_dict=None, lab_name=None ):
    '''
    This method accepts a Gene object located on an Assembly object (from things.py) and prints
    the gene in BED format.
    '''
    if gene is None:
        raise Exception( "ERROR: The print_biogene() function requires a biogene to be passed via the 'gene' argument" );

    ## we can auto-detect the molecule if the user didn't pass one
    #   and if there's only one.
    if on is None:
        on = gene.location().on

    gene_loc = gene.location_on( on )

    if gene.locus_tag is None:
        locus_tag = '.'
    else:
        locus_tag = gene.locus_tag

    

    fh.write("{0}\t{1}\t{2}\t{3}\t.\t{4}\t".format(
        on.id, gene_coords[0], gene_coords[1], locus_tag, gene_loc.
    ))

    for RNA in gene.RNAs():
        RNA_loc = RNA.location_on(on)

        if RNA_loc is None:
            raise Exception("ERROR: Expected RNA {0} to be located on {1} but it wasn't".format(RNA.id, on.id))

        # In this format, the RNA is supposed to be provided as segmented coordinates
        exons_printed = 0
        exons = sorted(RNA.exons())
        if RNA_loc.strand == -1:
            exons.reverse()
            
        for exon in exons:
            exon_loc = exon.location_on(on)
            
            if exon_loc is None:
                raise Exception("ERROR: Expected exon {0} to be located on {1} but it wasn't".format(exon.id, on.id) )

            exon_coords = biocode.utils.interbase0_to_humancoords(exon_loc.fmin, exon_loc.fmax, exon_loc.strand)

            if exons_printed == 0:
                fh.write("{0}\t{1}\t{2}\n".format(exon_coords[0], exon_coords[1], RNA.__class__.__name__))
            else:
                fh.write("{0}\t{1}\n".format(exon_coords[0], exon_coords[1]))

            exons_printed += 1

        if annot is not None and exons_printed > 0:
            if RNA.__class__.__name__ == 'mRNA':
                fh.write("\t\t\tprotein_id\t{0}\n".format(official_protein_id))
                
            fh.write("\t\t\ttranscript_id\t{0}\n".format(official_transcript_id))
            fh.write("\t\t\tproduct\t{0}\n".format(annot.product_name))
            if 'Note' in annot.other_attributes:
                # Sometimes this is a string, sometimes a list.  Only strings have the lower function:
                if hasattr(annot.other_attributes['Note'], 'lower'):
                    fh.write("\t\t\tnote\t{0}\n".format(annot.other_attributes['Note']))
                else:
                    fh.write("\t\t\tnote\t{0}\n".format(",".join(annot.other_attributes['Note'])))

        CDS_printed = 0
        CDSs = sorted(RNA.CDSs())
        if RNA_loc.strand == -1:
            CDSs.reverse()
        for CDS in CDSs:
            CDS_loc = CDS.location_on(on)

            if CDS_loc is None:
                raise Exception("ERROR: Expected CDS {0} to be located on {1} but it wasn't".format(CDS.id, on.id) )

            CDS_coords = biocode.utils.interbase0_to_humancoords(CDS_loc.fmin, CDS_loc.fmax, CDS_loc.strand)

            if CDS_printed == 0:
                fh.write("{0}\t{1}\tCDS\n".format(CDS_coords[0], CDS_coords[1]))
            else:
                fh.write("{0}\t{1}\n".format(CDS_coords[0], CDS_coords[1]))

            CDS_printed += 1

        if annot is not None and CDS_printed > 0:
            fh.write("\t\t\tprotein_id\t{0}\n".format(official_protein_id))
            fh.write("\t\t\ttranscript_id\t{0}\n".format(official_transcript_id))
            fh.write("\t\t\tproduct\t{0}\n".format(annot.product_name))

            if len(annot.ec_numbers) > 0:
                for ec_number in annot.ec_numbers:
                    fh.write("\t\t\tEC_number\t{0}\n".format(ec_number.number))
        
            if obo_dict is not None and len(annot.go_annotations) > 0:
                for go_term in annot.go_annotations:
                    go_id = "GO:{0}".format(go_term.go_id)
                    
                    """
                    There should actually be a distinction here between ISA (from blast/rapsearch matches)
                    and ISM (from HMM searches) but this is currently not encoded into the GFF, which is
                    our only source so far.  TODO: devise a way to pass this through.  For now, safer
                    to just use the all-encompassing 'ISS'

                    Reference:  http://geneontology.org/page/iss-inferred-sequence-or-structural-similarity
                    """
                    if go_id in obo_dict:
                        if obo_dict[go_id]['n'] == 'molecular_function':
                            fh.write("\t\t\tgo_function\t{1}|{0}||ISS\n".format(go_term.go_id, obo_dict[go_id]['p']))
                        elif obo_dict[go_id]['n'] == 'biological_process':
                            fh.write("\t\t\tgo_process\t{1}|{0}||ISS\n".format(go_term.go_id, obo_dict[go_id]['p']))
                        elif obo_dict[go_id]['n'] == 'cellular_component':
                            fh.write("\t\t\tgo_component\t{1}|{0}||ISS\n".format(go_term.go_id, obo_dict[go_id]['p']))
                    else:
                        raise Exception("ERROR: Didn't find GO id ({0}) in passed OBO file".format(go_id))

            if len(annot.dbxrefs) > 0:
                for dbxref in annot.dbxrefs:
                    fh.write("\t\t\tdb_xref\t{0}:{1}\n".format(dbxref.db, dbxref.identifier))
        

    











