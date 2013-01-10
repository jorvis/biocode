
print("The biothings.py is still under testing and development.  Please feel free to try using it, though the API is in flux.")

## Yes, there's a lot of abstraction that could happen here still but I'm trying to allow for
##  special consideration to be allowed for each feature type.  For now, this is just the skeleton.

class LocatableThing:
    '''
    Base class for any Things that are locatable on other Things, such as genes on contigs.
    '''
    def __init__( self, locations=None ):
        self.locations = locations

        if self.locations is None:
            self.locations = list()

    def locate_on( self, target=None, fmin=None, fmax=None, strand=None ):
        loc = Location(on=target, fmin=fmin, fmax=fmax, strand=strand)
        self.locations.append( loc )

    def _overlaps( self, refloc, qryloc ):
        ''' This assumes they're on the same molecule and only looks at coordinates, ignoring strandedness '''
        they_overlap = True  ## it's easier to negate an overlap

        if refloc.fmax <= qryloc.fmin or qryloc.fmax <= refloc.fmin:
            they_overlap = False

        return they_overlap

    def overlaps_with(self, thing=None):
        ## see if either are located on the same thing
        for ref_location in self.locations:
            for qry_location in thing.locations:
                if ref_location.on.id == qry_location.on.id:
                    ## if so, see if they have overlapping coordinates
                    if self._overlaps( ref_location, qry_location ):
                        return True

        ## if there were an overlap we would have reached it when looping, so now return False
        return False

class Location:
    def __init__( self, on=None, fmin=None, fmax=None, strand=None ):
        self.on = on
        self.fmin = fmin
        self.fmax = fmax
        self.strand = strand

class Assembly( LocatableThing ):
    def __init__( self, id=None, length=None, children=None ):
        self.id = id
        self.length = length
        self.children = children

        ## initialize any types needed
        self.children = _initialize_type_list(self.children, 'gene')
        
    def add_gene( self, gene ):
        self.children['gene'].append(gene)

    def genes(self):
        '''
        Returns all Gene objects which are children of this Assembly
        '''
        return self.children['gene']


class CDS( LocatableThing ):
    def __init__( self, id=None, locations=None, parent=None, length=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.length = length

class Exon( LocatableThing ):
    def __init__( self, id=None, locations=None, parent=None, length=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.length = length

class Gene( LocatableThing ):
    def __init__( self, id=None, locations=None, children=None ):
        super().__init__(locations)
        self.id = id
        self.children = children

        ## initialize any types needed
        self.children = _initialize_type_list(self.children, 'mRNA')
        self.children = _initialize_type_list(self.children, 'rRNA')
        self.children = _initialize_type_list(self.children, 'tRNA')

    def add_mRNA(self, rna):
        self.children['mRNA'].append(rna)
    
    def add_rRNA(self, rna):
        self.children['rRNA'].append(rna)

    def add_tRNA(self, rna):
        self.children['tRNA'].append(rna)


class RNA( LocatableThing ):
    def __init__( self, id=None, locations=None, parent=None, children=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.children = children

        ## initialize any types needed
        self.children = _initialize_type_list(self.children, 'exon')
        self.children = _initialize_type_list(self.children, 'CDS')

    def add_exon(self, exon):
        self.children['exon'].append(exon)

    def add_CDS(self, cds):
        self.children['CDS'].append(cds)  
        
class mRNA( RNA ):
    def __init__( self, id=None, locations=None, parent=None, children=None ):
        super().__init__(id, locations, parent, children)

class rRNA( RNA ):
    def __init__( self, id=None, locations=None, parent=None, children=None ):
        super().__init__(id, locations, parent, children)

class tRNA( RNA ):
    def __init__( self, id=None, locations=None, parent=None, children=None ):
        super().__init__(id, locations, parent, children)


#############################
## Private help functions
##   Hopefully this isn't to un-pythonic.  These were to aide in abstraction within classes
#############################

def _initialize_type_list( children, feattype ):
    if children is None:
        children = dict()
    
    if feattype not in children:
        children[feattype] = list()

    return children


def _print_thing( thing ):
    '''
    This is mostly used for debugging.  Pass any of the biothings like Gene, CDS, etc
    and this will print out its attributes including locations, children info, etc.
    '''
    print("DEBUG: Info for biothing with ID:({0})".format( thing.id ) )

    if hasattr(thing, 'parent'):
        print("\tparent = {0}".format(thing.parent.id) )
    else:
        print("\tparent = None")

    if hasattr(thing, 'children'):
        print("\tchild count = {0}".format( len(thing.children) ) )
    else:
        print("\tchildren = None")

    if hasattr(thing, 'locations'):
        print("\tlocation count: {0}".format(len(thing.locations)) )
        print("\tlocations:")

        for loc in thing.locations:
            print("\t\t{0}\t{1}\t{2}\t{3}".format(loc.on.id, loc.fmin, loc.fmax, loc.strand) )
    else:
        print("\tlocations = None")
