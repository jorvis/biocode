
print("The biothings.py is still under testing and development.  Please feel free to try using it, though the API is in flux.")

## Yes, there's a lot of abstraction that could happen here still but I'm trying to allow for
##  special consideration to be allowed for each feature type.  For now, this is just the skeleton.

class LocatableThing:
    '''
    Base class for any Things that are locatable on other Things, such as genes on contigs.
    '''
    def __init__( self, locations=list() ):
        self.locations = locations

    def locate_on(self, target=None, fmin=None, fmax=None, strand=None):
        loc = Location(on=target, fmin=fmin, fmax=fmax, strand=strand)
        self.locations.append( loc )

    def _overlaps( self, refloc, qryloc ):
        ''' This assumes they're on the same molecule and only looks at coordinates '''
        they_overlap = False

        if (qryloc.fmin >= refloc.fmin and qryloc.fmin <= refloc.fmax):
            they_overlap = True

        return they_overlap

    def overlaps_with(self, thing=None):
        overlap = None

        ## see if either are located on the same thing
        for ref_location in self.locations:
            for qry_location in thing.locations:
                if ref_location.on == qry_location.on:
                    ## if so, see if they have overlapping coordinates
                    if self._overlaps( ref_location, qry_location ):
                        return True
        
        return overlap

class Location:
    def __init__( self, on=None, fmin=None, fmax=None, strand=None ):
        self.on = on
        self.fmin = fmin
        self.fmax = fmax
        self.strand = strand

class Assembly:
    def __init__( self, id=None, length=None, children=dict() ):
        self.id = id
        self.length = length

        ## initialize any types needed
        _initialize_type_list(children, 'gene')
        
        self.children = children

    def add_gene( self, gene ):
        self.children['gene'].append(gene)

    def genes(self):
        '''
        Returns all Gene objects which are children of this Assembly
        '''
        return self.children['gene']


class CDS( LocatableThing ):
    def __init__( self, id=None, locations=list(), parent=None, length=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.length = length

class Exon( LocatableThing ):
    def __init__( self, id=None, locations=list(), parent=None, length=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.length = length

class Gene( LocatableThing ):
    def __init__( self, id=None, locations=list(), children=dict() ):
        super().__init__(locations)
        self.id = id

        ## initialize any types needed
        _initialize_type_list(children, 'mRNA')
        _initialize_type_list(children, 'rRNA')
        _initialize_type_list(children, 'tRNA')
        self.children = children

    def add_mRNA(self, rna):
        self.children['mRNA'].append(rna)
    
    def add_rRNA(self, rna):
        self.children['rRNA'].append(rna)

    def add_tRNA(self, rna):
        self.children['tRNA'].append(rna)


class RNA( LocatableThing ):
    def __init__( self, id=None, locations=list(), parent=None, children=dict() ):
        super().__init__(locations)
        self.id = id
        self.parent = parent

        ## initialize any types needed
        _initialize_type_list(children, 'exon')
        _initialize_type_list(children, 'CDS')
        
        self.children = children

    def add_exon(self, exon):
        self.children['exon'].append(exon)

    def add_CDS(self, cds):
        self.children['CDS'].append(cds)  
        
class mRNA( RNA ):
    def __init__( self, id=None, locations=list(), parent=None, children=dict() ):
        super().__init__(id, locations, parent, children)

class rRNA( RNA ):
    def __init__( self, id=None, locations=list(), parent=None, children=dict() ):
        super().__init__(id, locations, parent, children)

class tRNA( RNA ):
    def __init__( self, id=None, locations=list(), parent=None, children=dict() ):
        super().__init__(id, locations, parent, children)


#############################
## Private help functions
##   Hopefully this isn't to un-pythonic.  These were to aide in abstraction within classes
#############################

def _initialize_type_list( children, feattype ):
    if feattype not in children:
        children[feattype] = list()
