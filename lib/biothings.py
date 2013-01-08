
print("The biothings.py is still under testing and development.  Please feel free to try using it, though the API is in flux.")

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
    def __init__( self, id=None, length=None, children=list() ):
        self.id = id
        self.length = length

    def add_gene( self, gene ):
        self.children.append(gene)

class CDS( LocatableThing ):
    def __init__( self, id=None, locations=list(), parent=None, length=None ):
        super().__init__(locations)
        self.id = id
        self.length = length

class Exon( LocatableThing ):
    def __init__( self, id=None, locations=list(), parent=None, length=None ):
        super().__init__(locations)
        self.id = id
        self.length = length

class Gene( LocatableThing ):
    def __init__( self, id=None, locations=list(), children=list() ):
        super().__init__(locations)
        self.id = id
        self.children = children

    def add_mRNA(self, rna):
        self.children.append(rna)

class mRNA( LocatableThing ):
    def __init__( self, id=None, locations=list(), parent=None, children=list() ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.children = children

    def add_exon(self, exon):
        self.children.append(exon)

    def add_CDS(self, cds):
        self.children.append(cds)

class tRNA( LocatableThing ):
    def __init__( self, id=None, locations=list(), parent=None, children=list() ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.children = children

    def add_exon(self, exon):
        self.children.append(exon)

    def add_CDS(self, cds):
        self.children.append(cds)
