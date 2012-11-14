
class Assembly:
    def __init__( self, id=None, length=None ):
        self.id = id
        self.length = length
        print("initialized an Assembly with id {0}".format(self.id))

class Gene:
    def __init__( self, id=None, locations=list(), mRNAs=list() ):
        self.id = id
        self.locations = locations
        self.mRNAs = mRNAs

    def add_mRNA(self, rna):
        self.mRNAs.append(rna)

    def locate_on(self, assembly=None, fmin=None, fmax=None, strand=None):
        loc = Location(on=assembly, fmin=fmin, fmax=fmax, strand=strand)
        self.locations.append( loc )

    def _overlap( self, refloc, qryloc ):
        ''' This assumes they're on the same molecule and only looks at coordinates '''
        they_overlap = False

        if (qryloc.fmin >= refloc.fmin and qryloc.fmin <= refloc.fmax):
            they_overlap = True

        return they_overlap
        
    def overlap_with(self, gene=None):
        overlap = None

        ## see if either are located on the same molecule
        for ref_location in self.locations:
            for qry_location in gene.locations:
                if ref_location.on == qry_location.on:
                    if self._overlap( ref_location, qry_location ):
                        return True
        
        return overlap

class mRNA:
    def __init__( self, id=None, locations=list(), parent=None ):
        self.id = id
        self.locations = locations
        self.parent = parent

    def locate_on(self, assembly=None, fmin=None, fmax=None, strand=None):
        loc = Location(on=assembly, fmin=fmin, fmax=fmax, strand=strand)
        self.locations.append( loc )
    
class Location:
    def __init__( self, on=None, fmin=None, fmax=None, strand=None ):
        self.on = on
        self.fmin = fmin
        self.fmax = fmax
        self.strand = strand

class Exon:
    pass

class CDS:
    pass
