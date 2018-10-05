from biocode import gff

def print_biogene( gene=None, fh=None, source=None, on=None ):
    '''
    This method accepts a Gene object located on an Assembly object (from things.py) and prints
    the feature graph for that gene in NCBI-spec GFF3 format, including the gene, mRNA, CDS and 
    exon features.
    '''
    ## handle defaults
    if source is None:
        source = '.'

    if gene is None:
        raise Exception( "ERROR: The print_biogene() function requires a biogene to be passed via the 'gene' argument" );

    ## we can auto-detect the molecule if the user didn't pass one
    #   and if there's only one.
    if on is None:
        on = gene.location().on

    gene_loc = gene.location_on( on )
    strand = '0'

    if gene_loc.strand == 1:
        strand = '+'
    elif gene_loc.strand == -1:
        strand = '-'

    gene_annot_atts = dict()

    if gene.locus_tag is None:
        raise Exception("""
           In order to create NCBI-spec GFF3 gene features must have locus_tag attributes.  If you
           need help, the biocode/gff/add_gff3_locus_tags.py script does this.
        """)

    gene_partiality_string = gff._partiality_string(gene_loc)
    if gene_partiality_string is not None:
        gene_annot_atts['Partial'] = gene_partiality_string

    # we need to pull the gene symbol from the RNA (for non-coding features) or the polypeptide
    rnas = sorted(gene.RNAs())
    if rnas[0].annotation and rnas[0].annotation.gene_symbol:
        gene_annot_atts['gene'] = rnas[0].annotation.gene_symbol
    else:
        polypeptides = sorted(rnas[0].polypeptides())
        if len(polypeptides) and polypeptides[0].annotation:
            gene_annot_atts['gene'] = polypeptides[0].annotation.gene_symbol

    columns = ['.']*9
    columns[0:3] = [gene_loc.on.id, source, 'gene']
    columns[3:7] = [str(gene_loc.fmin + 1), str(gene_loc.fmax), '.', strand]
    columns[8] = gff.build_column_9( id=gene.id, parent=None, other=gene_annot_atts )

    ## print the gene line
    fh.write( "\t".join(columns) + "\n" )

    for RNA in rnas:
        RNA_loc = RNA.location_on( on )

        if RNA_loc is None:
            raise Exception("ERROR: Expected RNA {0} to be located on {1} but it wasn't".format(RNA.id, on.id))

        rna_annot_atts = dict()
        if RNA.locus_tag is not None:
            rna_annot_atts['locus_tag'] = RNA.locus_tag

        # NCBI specification requires RNA features to have an transcript_id attribute
        rna_annot_atts['transcript_id'] = RNA.id

        rna_partiality_string = gff._partiality_string(RNA_loc)
        if rna_partiality_string is not None:
            rna_annot_atts['Partial'] = rna_partiality_string

        if not RNA.annotation:
            polypeptides = RNA.polypeptides()
            if len(polypeptides) and polypeptides[0].annotation:
                RNA.annotation = polypeptides[0].annotation

        # ncRNAs might have annotation elements at the RNA level
        if RNA.annotation:
            # tRNAs should have a product and anticodon
            if RNA.annotation.product_name:
                rna_annot_atts['product'] = RNA.annotation.product_name

            if hasattr(RNA, 'anticodon') and RNA.anticodon:
                rna_annot_atts['anticodon'] = RNA.anticodon

            if RNA.annotation.gene_symbol:
                rna_annot_atts['gene'] = RNA.annotation.gene_symbol

            dbxref_strs = list()
            for dbxref in RNA.annotation.dbxrefs:
                dbxref_strs.append("{0}:{1}".format(dbxref.db, dbxref.identifier))

            if len(dbxref_strs) > 0:
                rna_annot_atts['Dbxref'] = ",".join(dbxref_strs)

        columns[2] = RNA.__class__.__name__
        columns[3:5] = [str(RNA_loc.fmin + 1), str(RNA_loc.fmax)]
        columns[8] = gff.build_column_9( id=RNA.id, parent=RNA.parent.id, other=rna_annot_atts )
        fh.write( "\t".join(columns) + "\n" )

        ## Write any 5' UTRs
        for utr in sorted(RNA.five_prime_UTRs()):
            utr_loc = utr.location_on( on )

            if utr_loc is None:
                raise Exception("ERROR: Expected UTR {0} to be located on {1} but it wasn't".format(utr.id, on.id) )

            columns[2] = 'five_prime_UTR'
            columns[3:5] = [str(utr_loc.fmin + 1), str(utr_loc.fmax)]
            columns[7] = str(utr_loc.phase)
            columns[8] = gff.build_column_9( id=utr.id, parent=utr.parent.id, other=None )
            fh.write( "\t".join(columns) + "\n" )

        ## handle each CDS for this mRNA
        for CDS in sorted(RNA.CDSs()):
            CDS_loc = CDS.location_on( on )

            if CDS_loc is None:
                raise Exception("ERROR: Expected CDS {0} to be located on {1} but it wasn't".format(CDS.id, on.id) )

            cds_partiality_string = gff._partiality_string(CDS_loc)
            cds_annot_atts = dict()
            if cds_partiality_string is not None:
                cds_annot_atts['Partial'] = cds_partiality_string

            columns[2] = 'CDS'
            columns[3:5] = [str(CDS_loc.fmin + 1), str(CDS_loc.fmax)]
            columns[7] = str(CDS_loc.phase)
            columns[8] = gff.build_column_9( id=CDS.id, parent=RNA.id, other=cds_annot_atts )
            fh.write( "\t".join(columns) + "\n" )

        columns[7] = '.'

        ## handle each exon for this RNA
        for exon in sorted(RNA.exons()):
            exon_loc = exon.location_on( on )

            if exon_loc is None:
                raise Exception("ERROR: Expected exon {0} to be located on {1} but it wasn't".format(exon.id, on.id))

            exon_partiality_string = gff._partiality_string(exon_loc)
            exon_annot_atts = dict()
            if exon_partiality_string is not None:
                exon_annot_atts['Partial'] = exon_partiality_string

            columns[2] = 'exon'
            columns[3:5] = [str(exon_loc.fmin + 1), str(exon_loc.fmax)]
            columns[8] = gff.build_column_9( id=exon.id, parent=RNA.id, other=exon_annot_atts )
            fh.write( "\t".join(columns) + "\n" )

        # are there polypeptides?
        for polypeptide in sorted(RNA.polypeptides()):
            if len(polypeptide.locations) == 0:
                polypeptide_loc = RNA_loc
            else:
                polypeptide_loc = polypeptide.location_on(on)

            annot = polypeptide.annotation
            assertions = dict()
            dbxref_strs = list()
            ontology_strs = list()

            if annot is not None:
                if annot.product_name is not None:
                    assertions['product'] = annot.product_name

                if annot.gene_symbol is not None:
                    assertions['gene_symbol'] = annot.gene_symbol

                # should this construction be in build_column_9() instead?
                # this is treating each EC number as a string, but each should be an ECAnnotation object
                for ec in annot.ec_numbers:
                    dbxref_strs.append("EC:{0}".format(ec.number))

                for go in annot.go_annotations:
                    ontology_strs.append("GO:{0}".format(go.go_id))

                for dbxref in annot.dbxrefs:
                    dbxref_strs.append("{0}:{1}".format(dbxref.db, dbxref.identifier))

                # are there other attributes to handle?
                for att in annot.other_attributes:
                    assertions[att] = annot.other_attributes[att]

            if len(dbxref_strs) > 0:
                assertions['Dbxref'] = ",".join(dbxref_strs)

            if len(ontology_strs) > 0:
                assertions['Ontology_term'] = ",".join(ontology_strs)

            if 'Partial' in rna_annot_atts:
                assertions['Partial'] = rna_annot_atts['Partial']

            columns[2] = 'polypeptide'
            columns[3:5] = [str(polypeptide_loc.fmin + 1), str(polypeptide_loc.fmax)]
            columns[8] = gff.build_column_9( id=polypeptide.id, parent=polypeptide.parent.id, other=assertions )
            fh.write( "\t".join(columns) + "\n" )

        ## Write any 5' UTRs
        for utr in sorted(RNA.three_prime_UTRs()):
            utr_loc = utr.location_on( on )

            if utr_loc is None:
                raise Exception("ERROR: Expected UTR {0} to be located on {1} but it wasn't".format(utr.id, on.id) )

            columns[2] = 'three_prime_UTR'
            columns[3:5] = [str(utr_loc.fmin + 1), str(utr_loc.fmax)]
            columns[7] = str(utr_loc.phase)
            columns[8] = gff.build_column_9( id=utr.id, parent=utr.parent.id, other=None )
            fh.write( "\t".join(columns) + "\n" )

