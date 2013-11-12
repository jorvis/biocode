import re
import biothings
import bioannotation
from urllib.parse import unquote

def build_column_9( id=None, parent=None, other=None ):
    ## either the id or parent must be defined
    if id is None and parent is None:
        raise Exception("ERROR: attempt to build a GFF3 9th column with no ID or Parent attributes")

    colstring = None

    ## ID should always be first, if there is one
    if id is not None:
        colstring = "ID={0}".format(id)

    if parent is not None:
        if colstring is not None:
            colstring += ";"

        colstring += "Parent={0}".format(parent)

    # other is a dict of key,value pairs
    if other is not None:
        for att in other:
            if colstring is not None:
                colstring += ";"

            colstring += "{0}={1}".format(att, other[att])
        
    return colstring


def column_9_value(value, key):
    '''
    Pass a case-sensitive key and this function will return the value,
    if present, else None.
    
    This was borrowed from the URL below and then modified for Python 3
    ftp://ftp.informatics.jax.org/%2Fpub/mgigff/gff3.py
    '''
    SEMI	= ';'
    COMMA	= ','
    EQ	        = '='
    WSP_RE = re.compile(r'^\s*$')
    
    if value == ".":
        return {}
    c9 = {}
    for t in value.split(SEMI):
        if WSP_RE.match(t):
            continue
        tt = t.split(EQ)
        if len(tt) != 2:
            raise Exception("Bad column 9 format: {0}".format(value) )
        n = unquote(tt[0].strip())
        [*v] = map(unquote, tt[1].strip().split(COMMA))
        if len(v) == 1:
            v = v[0]
        c9[n] = v

    if key in c9:
        # this is a list if there were comma-separated values
        return c9[key]
    else:
        return None


def get_gff3_features(gff3_file):
    '''
    Parses the passed GFF3 file and returns two dicts, loaded with biocode.biothings objects:

    1. The first dict are the Assembly objects, keyed on assembly ID.  Each Assembly has all of the
       children populated, so you can fully recover gene, RNA, exon and CDS features iterating on
       the assembly.
    2. The second dist is a flat structure of all the descendent feature objects of the Assemblies
       keyed by the feature IDs.

    See the documentation for each feature type in biocode.biothings for more info
    '''
    
    assemblies = dict()
    features   = dict()

    # these are related to parsing any embedded FASTA
    in_fasta_section = False
    is_assembly_fasta = False
    current_fasta_id = None
    
    
    for line in open(gff3_file):
        if in_fasta_section == True:
            m = re.search('>(\S+)\s*(.*)', line)
            if m:
                current_fasta_id = m.group(1)

                if current_fasta_id in assemblies:
                    is_assembly_fasta = True
                else:
                    is_assembly_fasta = False
                    
            else:
                if is_assembly_fasta == True:
                    # must be a sequence line for an assembly
                    # python 2.6+ makes string concatenation amortized O(n)
                    #  http://stackoverflow.com/a/4435752/1368079
                    assemblies[current_fasta_id].residues += str(line.rstrip())
                    assemblies[current_fasta_id].length = len( assemblies[current_fasta_id].residues )

            continue
        
        elif line.startswith("##FASTA"):
            # all data to the end of the file must be FASTA
            in_fasta_section = True
            continue

        
        cols = line.split("\t")

        if len(cols) != 9:
            continue

        mol_id = cols[0]
        
        # initialize this assembly if we haven't seen it yet
        if mol_id not in assemblies:
            assemblies[mol_id] = biothings.Assembly( id=mol_id, residues='' )

        current_assembly = assemblies[mol_id]
        rfmin = int(cols[3]) - 1
        rfmax = int(cols[4])
        rstrand = None
        feat_id = column_9_value(cols[8], 'ID')
        parent_id = column_9_value(cols[8], 'Parent')
        parent_feat = None

        # shared features is not yet supported
        if isinstance(parent_id, list):
            raise Exception("This line contains a shared feature with multiple parents.  This isn't yet supported:\n{0}".format(line))

        if parent_id is not None:
            if parent_id in features:
                parent_feat = features[parent_id]
            else:
                raise Exception("Error in GFF3: Parent {0} referenced by a child feature before it was defined".format(parent_id) )

        if cols[6] == '-':
            rstrand = -1
        elif cols[6] == '+':
            rstrand = 1
        else:
            rstrand = 0
            
        if cols[2] == 'gene':
            gene = biothings.Gene(id=feat_id)
            gene.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            features[feat_id] = gene
            current_assembly.add_gene(gene)
        
        elif cols[2] == 'mRNA':
            mRNA = biothings.mRNA(id=feat_id, parent=parent_feat)
            mRNA.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            parent_feat.add_mRNA(mRNA)
            features[feat_id] = mRNA

        elif cols[2] == 'rRNA':
            rRNA = biothings.rRNA(id=feat_id, parent=parent_feat)
            rRNA.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            parent_feat.add_rRNA(rRNA)
            features[feat_id] = rRNA
            
        elif cols[2] == 'tRNA':
            tRNA = biothings.tRNA(id=feat_id, parent=parent_feat)
            tRNA.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            parent_feat.add_tRNA(tRNA)
            features[feat_id] = tRNA

        elif cols[2] == 'exon':
            exon = biothings.Exon(id=feat_id, parent=parent_feat)
            exon.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            parent_feat.add_exon(exon)
            features[feat_id] = exon

        elif cols[2] == 'CDS':
            CDS = biothings.CDS(id=feat_id, parent=parent_feat)
            CDS.locate_on(target=current_assembly, fmin=rfmin, fmax=rfmax, strand=rstrand)
            parent_feat.add_CDS(CDS)
            features[feat_id] = CDS

    return (assemblies, features)


def parse_gff3_by_relationship( gff3_file ):
    '''
    Parses a GFF3 file, caring only about the ID/Parent relationships and returning a dict where
    they keys are the feature IDs and the values are an dict of dicts.   The depth of this depends
    on the depth of parent/child relationships, but every node looks like this:

    $feat_id = { children: [], fmin: 30, cols: [] }

    Every feature must have either and ID or Parent defined, else an error will be thrown.
    
    The structure returned looks like this:
    {
    $molecule_id = { 
                         fmin: 30,
                         cols: [...],
                         children: {
                             same as this structure
                         }

                   }, ...
    

    }

    WARNING:
    This will currently fail on any GFF that:
       - allows the same IDs on multiple lines, such as those for  discontiguous features (allowed in spec)
       - features with shared parents


    '''
    ## might try turning this into a defaultdict at a later point: http://ohuiginn.net/mt/2010/07/nested_dictionaries_in_python.html
    fgraph = dict()
    current_molecule_id = None
    ## key = id, value = parent_id
    parentage = dict()
    ## where features go during processing that have a parent which hasn't been seen yet
    purgatory = list()

    for line in open(gff3_file):
        cols = line.split("\t")

        if len(cols) == 9:
            cols[8] = cols[8].rstrip()
        else:
            continue

        mol_id = cols[0]
        id = column_9_value(cols[8], 'ID')
        # shared parents will cause failure here
        parent = column_9_value(cols[8], 'Parent')

        if id:
            parentage[id] = parent
    
        if mol_id != current_molecule_id:
            if current_molecule_id is not None:
                _reunite_children( fgraph, current_molecule_id, purgatory )

            ## reset for the next molecule
            purgatory = list()
            current_molecule_id = mol_id

            if current_molecule_id not in fgraph:
                fgraph[current_molecule_id] = dict()

        molecule = fgraph[current_molecule_id]

        if parent:
            uparent = _get_ultimate_parent( parentage, parent )
            molecule[uparent]['children'].append( {'id': id, 'cols': cols} )
        else:
            if id is None:
                raise Exception("ERROR: Encountered a line without a Parent or ID assigned: {0}".format(line))

            ## this shouldn't exist already
            if id in molecule:
                raise Exception("ERROR: found duplicate id ({0}) in input file".format(id) )

            molecule[id] = {'fmin': cols[3], 'children': list(), 'cols': cols}

    ## don't forget to handle the last one
    _reunite_children( fgraph, current_molecule_id, purgatory )

    return fgraph


def print_biogene( gene=None, fh=None, source=None, on=None ):
    '''
    This method accepts a Gene object located on an Assembly object (from biothings.py) and prints
    the feature graph for that gene in GFF3 format, including the gene, mRNA, CDS and exon features.
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

    columns = ['.']*9
    columns[0:3] = [gene_loc.on.id, source, 'gene']
    columns[3:7] = [str(gene_loc.fmin + 1), str(gene_loc.fmax), '.', strand]
    columns[8] = build_column_9( id=gene.id, parent=None, other=None )

    ## print the gene line
    fh.write( "\t".join(columns) + "\n" )

    ## modifications for the mRNA
    for mRNA in gene.mRNAs():
        mRNA_loc = mRNA.location_on( on )

        if mRNA_loc is None:
            raise Exception("ERROR: Expected mRNA {0} to be located on {1} but it wasn't".format(mRNA.id, on.id))
        
        columns[2] = 'mRNA'
        columns[3:5] = [str(mRNA_loc.fmin + 1), str(mRNA_loc.fmax)]
        columns[8] = build_column_9( id=mRNA.id, parent=mRNA.parent.id, other=None )
        fh.write( "\t".join(columns) + "\n" )
    
        ## handle each CDS for this mRNA
        for CDS in mRNA.CDSs():
            CDS_loc = CDS.location_on( on )

            if CDS_loc is None:
                raise Exception("ERROR: Expected CDS {0} to be located on {1} but it wasn't".format(CDS.id, on.id) )
            
            columns[2] = 'CDS'
            columns[3:5] = [str(CDS_loc.fmin + 1), str(CDS_loc.fmax)]
            columns[7] = str(CDS_loc.phase)
            columns[8] = build_column_9( id=CDS.id, parent=mRNA.id, other=None )
            fh.write( "\t".join(columns) + "\n" )

        columns[7] = '.'

        ## handle each exon for this mRNA
        for exon in mRNA.exons():
            exon_loc = exon.location_on( on )

            if exon_loc is None:
                raise Exception("ERROR: Expected exon {0} to be located on {1} but it wasn't".format(exon.id, on.id))
            
            columns[2] = 'exon'
            columns[3:5] = [str(exon_loc.fmin + 1), str(exon_loc.fmax)]
            columns[8] = build_column_9( id=exon.id, parent=mRNA.id, other=None )
            fh.write( "\t".join(columns) + "\n" )

        # are there polypeptides?
        for polypeptide in mRNA.polypeptides():
            ## HACK - we probably don't want to keep this
            polypeptide_loc = mRNA_loc
            annot = polypeptide.annotation

            assertions = {'product_name':annot.product_name, }

            #if 

            columns[2] = 'polypeptide'
            columns[3:5] = [str(polypeptide_loc.fmin + 1), str(polypeptide_loc.fmax)]
            columns[8] = build_column_9( id=polypeptide.id, parent=mRNA.parent.id, other=assertions )
            fh.write( "\t".join(columns) + "\n" )

def _reunite_children( fg, mol_id, kids ):
    sys.stderr.write("WARNING: biocodegff._reunite_children called but not yet implemented.  You are ahead of your time.\n")
    pass


def _get_ultimate_parent( p, id ):
    if p is None:
        return id

    oldest = id

    while p[oldest] is not None:
        oldest = p[oldest]

    return oldest










