import re
import bioannotation
import biocodeutils
import biothings
import sys

"""
NOTE: This currently only works for coding features

EXAMPLES:

Ex1 (header portion of GenBank flat file):
            |--------------- 67 character limit ------------------------------|
LOCUS       HQ450365                 735 bp    DNA     linear   PLN 15-JUL-2011
DEFINITION  Rhizopus delemar voucher RA 99-880 RNA polymerase II second largest
            subunit (RPB2) gene, partial cds.
ACCESSION   HQ450365
VERSION     HQ450365.1  GI:315271506

Ex2 (feature table of GenBank flat file):
                     |--------------- 58 character limit ---------------------|
     gene            42455..44072
                     /locus_tag="TP03_0020"
     mRNA            join(42455..42702,42737..42841,42884..42995,43267..43441,
                     43498..43653,43685..43871,43988..44072)
                     /locus_tag="TP03_0020"
                     /product="methyltransferase, putative"
     CDS             join(42664..42702,42737..42841,42884..42995,43267..43441,
                     43498..43653,43685..43871,43988..44047)
                     /locus_tag="TP03_0020"
                     /EC_number="2.1.1.-"
                     /note="GO_function: GO:0008757 -
                     S-adenosylmethionine-dependent methyltransferase activity"
                     /codon_start=1
                     /product="methyltransferase, putative"
                     /protein_id="EAN30756.1"
                     /db_xref="GI:68349992"
                     /translation="MSIRPEHSAPPEIYYSAEESRKYNTNSHILKIQTQMSERALEML
                     LLPEDEMGLVLDIGCGTGISGNVISNSNHFWIGLDISQHMLHESLLNEIEGDVVLCDI
                     GEPMNFLPNMFDGCISISVLQWLFISNHKSQEPYHRLSSFFKWLYKSLAYNARACLQF
                     YPENVEQVDMLLDIVKKCNFNGGLVVDNPNSVKAKKYYLCIWSYNSNIYHKLPNPIEQ
                     NGDHEEVEFDEVESCVMKKRKNKLTYKDRIIKKKQQQRNKGMKTRPDTKYTGRKRPHA
                     F"
"""

# character limits for the content portions of the header and the feature table,
# respectively (as illustrated above):
MAX_CONTENT_WIDTH = 79
MAX_HEADER_CONTENT_WIDTH = 67
MAX_FTABLE_CONTENT_WIDTH = 58

HEADER_MARGIN_WIDTH = MAX_CONTENT_WIDTH - MAX_HEADER_CONTENT_WIDTH
HEADER_MARGIN = " " * HEADER_MARGIN_WIDTH

def print_biogene( gene=None, fh=None, on=None ):
    '''
    This method accepts a Gene object located on an Assembly object (from biothings.py) and prints
    the feature graph for that gene in Genbank flat file format, including the gene, RNA and CDS
    '''
    if gene is None:
        raise Exception( "ERROR: The print_biogene() function requires a biogene to be passed via the 'gene' argument" );

    ## we can auto-detect the molecule if the user didn't pass one
    #   and if there's only one.
    if on is None:
        on = gene.location().on

    gene_loc = gene.location_on( on )
    gene_start = gene_loc.fmin + 1
    gene_stop  = gene_loc.fmax

    # area to hack if you want to set default values, for debugging
    #gene.locus_tag = 'Tparva_0000002'

    if gene_loc.strand == 1:
        fh.write("     gene            {0}..{1}\n".format(gene_start, gene_stop))
    else:
        fh.write("     gene            complement({0}..{1})\n".format(gene_start, gene_stop))

    if gene.locus_tag is None:
        sys.stderr.write("WARNING: No locus_tag found on gene {0}\n".format(gene.id))
    else:
        fh.write("                     /locus_tag=\"{0}\"\n".format(gene.locus_tag))


    for mRNA in sorted(gene.mRNAs()):
        mRNA_loc = mRNA.location_on( on )

        ###########################
        ## write the mRNA feature (made up of exon fragments)
        mRNA_loc_segments = list()
        for exon in mRNA.exons():
            exon_loc = exon.location_on(on)
            mRNA_loc_segments.append( [exon_loc.fmin + 1, exon_loc.fmax] )

        mRNA_loc_string = segments_to_string(mRNA_loc_segments)

        if mRNA_loc.strand == 1:
            fh.write("     mRNA            {0}\n".format(mRNA_loc_string))
        else:
            fh.write("     mRNA            complement({0})\n".format(mRNA_loc_string))

        # Handle the locus tag, but we've already warned if not present on the gene, so don't
        #  do it again here.
        if gene.locus_tag is not None:
            fh.write("                     /locus_tag=\"{0}\"\n".format(gene.locus_tag))

        if mRNA.annotation is not None:
            # debug:  You can try out some annotation defaults for printing here
            mRNA.annotation.product_name = "Hypothetical protein"

            if mRNA.annotation.product_name is not None:
                fh.write("                     /product=\"{0}\"\n".format(mRNA.annotation.product_name))

        ###########################
        ## write the CDS feature (made up of CDS fragments)
        cds_loc_segments = list()
        for cds in mRNA.CDSs():
            cds_loc = cds.location_on(on)
            cds_loc_segments.append( [cds_loc.fmin + 1, cds_loc.fmax] )

        cds_loc_string = segments_to_string(cds_loc_segments)

        if cds_loc.strand == 1:
            fh.write("     CDS             {0}\n".format(cds_loc_string))
        else:
            fh.write("     CDS             complement({0})\n".format(cds_loc_string))

        # Handle the locus tag, but we've already warned if not present on the gene, so don't
        #  do it again here.
        if gene.locus_tag is not None:
            fh.write("                     /locus_tag=\"{0}\"\n".format(gene.locus_tag))

        ## if there is annotation on the polypeptide, include it here
        polypeptides = mRNA.polypeptides()
        if len(polypeptides) == 1 and polypeptides[0].annotation is not None:
            annot = polypeptides[0].annotation
            if annot.product_name is not None:
                fh.write("                     /product=\"{0}\"\n".format(annot.product_name))

        cds_residues = mRNA.get_CDS_residues()
        polypeptide_residues = biocodeutils.translate(cds_residues)

        if len(polypeptide_residues) > 0:
            #fh.write("                     /TRANSLATION=\"{0}\"\n".format(polypeptide_residues))
            
            # This is the easiest case first, where no wrapping is needed.
            if len(polypeptide_residues) < MAX_FTABLE_CONTENT_WIDTH - 15:
                fh.write("                     /translation=\"{0}\"\n".format(polypeptide_residues))
            else:
                # If we get here, we must wrap
                fh.write("                     /translation=\"{0}\n".format(polypeptide_residues[0:MAX_FTABLE_CONTENT_WIDTH - 14]))
                remaining = polypeptide_residues[MAX_FTABLE_CONTENT_WIDTH - 14:]
                closing_parens_written = False
                
                while len(remaining) > 0:
                    if len(remaining) > MAX_FTABLE_CONTENT_WIDTH - 1:
                        fh.write("                     {0}\n".format(remaining[0:MAX_FTABLE_CONTENT_WIDTH]))
                        remaining = remaining[MAX_FTABLE_CONTENT_WIDTH:]
                    else:
                        fh.write("                     {0}\"\n".format(remaining))
                        remaining = ""
                        closing_parens_written = True

                if closing_parens_written == False:
                    # G675_02159
                    fh.write("                     \"\n")

                



def segments_to_string(locs):
    lines = list()

    if len(locs) < 2:
        loc_string = ""
    else:
        loc_string = "join(";

    for loc in locs:
        new_part = "{0}..{1},".format(loc[0], loc[1])

        if len(loc_string) + len(new_part) > MAX_FTABLE_CONTENT_WIDTH:
            lines.append(loc_string)
            loc_string = ""

        loc_string += "{0}..{1},".format(loc[0], loc[1])

    loc_string = loc_string.rstrip(',')

    if len(locs) > 1:
        loc_string += ')'

    lines.append(loc_string)

    loc_string = "\n                     ".join(lines)

    return loc_string



def line_wrap_lineage_string(lineage):
    '''
    This method accepts a lineage string in the format 
    "Eukaryota; Alveolata; Apicomplexa; Aconoidasida; Piroplasmida; Theileriidae; Theileria"
    and inserts line breaks as needed.
    '''
    ll = len(lineage)
    result = ""
    is_first_line = True

    def add_line(new_line):
        nonlocal is_first_line
        nonlocal result
        if (is_first_line):
            is_first_line = False
        else:
            result = result + "\n" + HEADER_MARGIN
        result = result + new_line

    while (ll > 0):
        # remainder of string fits in one line
        if (ll <= MAX_HEADER_CONTENT_WIDTH):
            add_line(lineage)
            lineage = ""
        # remainder of string does not fit in one line
        else:
            # look for a space on which to break the line
            max_sp = MAX_HEADER_CONTENT_WIDTH 
            sp = -1
            for p in range(max_sp, 0, -1):
                if (lineage[p] == ' '):
                    sp = p
                    break
            # space not found, must break in the middle of a word
            if (sp == -1):
                sp = MAX_HEADER_CONTENT_WIDTH - 1
            add_line(lineage[0:sp+1])
            lineage = lineage[sp+1:]
        ll = len(lineage)
    return result


