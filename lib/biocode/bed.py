import sys

import biocode.utils
import biocode.things

def print_bed_from_assemblies(assemblies=None, ofh=None):
    """
    Utility function to write a BED file from a list() of biothings.Assembly objects.

    No headers are used so the output is compatible with other tools like bedToBigBed.

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

                    print_biogene(gene=new_gene, fh=ofh)

            if len(mRNAs) > 1:
                gene_loc = gene.location_on(current_assembly)
                mRNA_loc = mRNAs[0].location_on(current_assembly)
                gene_loc.fmin = mRNA_loc.fmin
                gene_loc.fmax = mRNA_loc.fmax
                gene_loc.strand = mRNA_loc.strand

            print_biogene(gene=gene, fh=ofh)


def print_biogene(gene=None, fh=None, on=None):
    '''
    This method accepts a Gene object located on an Assembly object (from things.py) and prints
    the gene in BED format.

    The 4th column of the BED output is to identify the gene.  By default, the locus identifier
    is used but if this isn't present we default to the gene's ID.
    '''
    if gene is None:
        raise Exception( "ERROR: The print_biogene() function requires a biogene to be passed via the 'gene' argument" );

    ## we can auto-detect the molecule if the user didn't pass one
    #   and if there's only one.
    if on is None:
        on = gene.location().on

    gene_loc = gene.location_on( on )

    locus_tag = gene.id if gene.locus_tag is None else gene.locus_tag
    gene_strand = '+' if gene_loc.strand == 1 else '-'

    for RNA in gene.RNAs():
        fh.write("{0}\t{1}\t{2}\t{3}\t0\t{4}\t".format(
            on.id, gene_loc.fmin, gene_loc.fmax, locus_tag, gene_strand
        ))

        RNA_loc = RNA.location_on(on)

        if RNA_loc is None:
            raise Exception("ERROR: Expected RNA {0} to be located on {1} but it wasn't".format(RNA.id, on.id))

        exons = sorted(RNA.exons())

        if RNA_loc.strand == -1:
            exons.reverse()

        block_starts = list()
        exon_lengths = list()

        for exon in exons:
            exon_loc = exon.location_on(on)
            block_starts.append(str(exon_loc.fmin - gene_loc.fmin))
            exon_lengths.append(str(exon_loc.fmax - exon_loc.fmin))
            
        if len(exons):
            thick_start = exons[0].location()
            thick_end   = exons[-1].location()
            fh.write("{0}\t{1}\t0\t{2}\t{3},\t{4},\n".format(
                thick_start.fmin, thick_end.fmax, len(exons), ','.join(exon_lengths), ','.join(block_starts)))
        else:
            # Catch here for genes without exons (tRNAs, etc.)
            fh.write("{0}\t{1}\t0\t{2}\t{3},\t0,\n".format(gene_loc.fmin, gene_loc.fmin, 1, 
                                                           gene_loc.fmax - gene_loc.fmin))



        

    











