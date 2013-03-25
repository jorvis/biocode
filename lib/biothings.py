import uuid
import biocodegff
import sys

'''
Warning: this module requires Python 3.2 or higher

Yes, there's a lot of abstraction that could happen here still but I'm trying to allow for 
special consideration to be allowed for each feature type and for calls to read more like
biology than CS.
'''

class LocatableThing:
    '''
    Base class for any Things that are locatable on other Things, such as genes on contigs.
    '''
    def __init__( self, locations=None ):
        self.locations = locations

        if self.locations is None:
            self.locations = list()

    def __lt__(self, other):
        return self.locations[0].fmin < other.locations[0].fmin

    def __le__(self, other):
        return self.locations[0].fmin <= other.locations[0].fmin

    def __eq__(self, other):
        return self.locations[0].fmin == other.locations[0].fmin

    def __ne__(self, other):
        return self.locations[0].fmin != other.locations[0].fmin

    def __gt__(self, other):
        return self.locations[0].fmin > other.locations[0].fmin

    def __ge__(self, other):
        return self.locations[0].fmin >= other.locations[0].fmin

    def contained_within( self, thing ):
        '''
        Returns True/False depending on whether the current LocatableThing is contained completely
        within the passed one.  An aligned transcript, for example, is hoped to be fully
        contained within an annotated gene.  Only the min/max bounds of the features are checked here
        currently, not internal substructures (like intron/exon bounds)

        Checks both to make sure they're both located on the same molecule before comparing coordinates,
        ignoring strandedness.
        '''
        ## see if either are located on the same thing
        for refloc in self.locations:
            for qryloc in thing.locations:
                if refloc.on.id == qryloc.on.id:
                    ## if so, see if self is completely within the passed object
                    if refloc.fmin >= qryloc.fmin and refloc.fmax <= qryloc.fmax:
                        return True

        return False

    def has_same_coordinates_as( self, on=None, thing=None ):
        '''
        Reports True/False depending on whether the calling object shares both start/stop coordinates
        with the passed 'thing' when considering locations 'on' the passed molecule biothing object.
        '''
        other = thing

        for this_loc in self.locations:
            for other_loc in other.locations:
                ## are the molecules the same and on the one requested?
                if this_loc.on.id == other_loc.on.id == on.id:
                    if this_loc.fmin == other_loc.fmin and this_loc.fmax == other_loc.fmax:
                        print("VALIDATE: {0} has same coordinates as {1} on {2}".format(self.id, other.id, on.id) )
                        return True

        ## if we got here, there wasn't a match
        return False


    def locate_on( self, target=None, fmin=None, fmax=None, phase=None, strand=None ):
        '''
        Locates any given biothing onto another by position range.  Use this to locate a gene
        onto an assembly, for example.  Any biothing can have multiple locations.

        Coordinates here are 0-based interbase.  Strand values are [1,0,-1].  If you pass
        "+" or "-" they'll be converted to 1 and -1, respectively.
        '''
        if strand == '+':
            strand = 1
        elif strand == '-':
            strand = -1

        loc = Location(on=target, fmin=fmin, fmax=fmax, strand=strand, phase=phase)
        self.locations.append( loc )

    def located_on( self ):
        '''
        Returns a dict of the molecules this feature is located on.  The keys are the identifiers
        and values are the objects that represent them.  This is distinct from the 'locations'
        attribute, which provides each individual location across all molecules.
        '''
        mols = dict()
        for loc in self.locations:
            mols[loc.on.id] = loc.on

        return mols
    
    def location( self ):
        '''
        In a majority of use-cases any given feature will only be located on one other feature.
        This method lets you access this location using a direct method call, but will raise
        an exception if called on a biothing with either 0 or more than 1 location.
        '''
        loc_count = len( self.locations )

        if loc_count == 1:
            return self.locations[0]
        elif loc_count == 0:
            raise Exception("ERROR: {0} has no locations defined".format(self.id))
        else:
            raise Exception("ERROR: {0} has more than one location defined.  You should use the location_on() method instead.".format(self.id))

    def location_on( self, mol ):
        '''
        Returns a Location object representing the location of the calling biothing on the passed
        molecule.  The name of this method implies that the object will be located on the molecule
        only a single time, and this method throws an exception if there are more than one.  If it
        isn't located on that molecule at all, None will be returned.
        '''
        loc_found = None
        for loc in self.locations:
            if loc.on.id == mol.id:
                if loc_found is None:
                    loc_found = loc
                else:
                    raise Exception("ERROR: multiple locations on the same molecule found for {0} on {1} when not expected".format(self.id, mol.id) )

        return loc_found


    def overlaps_with( self, thing ):
        '''
        Returns True/False depending on whether the two LocatableThing objects have
        even a single overlapping Location on the same molecule.

        Checks both to make sure they're both located on the same molecule before comparing coordinates,
        ignoring strandedness.
        '''
        ## see if either are located on the same thing
        for ref_location in self.locations:
            for qry_location in thing.locations:
                if ref_location.on.id == qry_location.on.id:
                    ## if so, see if they have overlapping coordinates
                    ## it's easier to negate an overlap
                    if ref_location.fmax <= qry_location.fmin or qry_location.fmax <= ref_location.fmin:
                        return False

        ## if there were an overlap we would have reached it when looping, so now return False
        return True


    For forward strand features, phase is counted from the start field. For reverse strand features,
    phase is counted from the end field.

    The phase is REQUIRED for all CDS features.
    '''
    def __init__( self, on=None, fmin=None, fmax=None, phase=None, strand=None ):
        self.on = on
        self.fmin = fmin
        self.fmax = fmax
        self.strand = strand
        self.phase = 0 if phase is None else phase


class Assembly( LocatableThing ):
    def __init__( self, id=None, locations=None, length=None, children=None ):
        super().__init__(locations)
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
    def __init__( self, id=None, locations=None, parent=None, phase=None, length=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.length = length
        self.phase = 0 if phase is None else phase


class Exon( LocatableThing ):
    def __init__( self, id=None, locations=None, parent=None, length=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.length = length


class Intron( LocatableThing ):
    '''
    Usually not represented directly.  Rather, these are generated on the fly
    as needed when the RNA.introns() method is called.

    There are some issues to resolve here.  For example, what's the parentage
    of an intron?  What are conceivable children?
    '''
    def __init__( self, id=None, locations=None, length=None ):
        super().__init__(locations)
        self.id = id
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

    def mRNAs(self):
        return self.children['mRNA']
    
    '''
    Prints the gene as a GFF3 entry, with all children.  Is printed to STDOUT
    unless the fh option is passed.  Really, only checks are done here before
    passing the info to the biocodegff module
    '''
    '''
    Overloading the print function allows users to directly print gene objects.
    The options passed will decide the target and format of the print.

        --fh = File handle to which we should print (default: STDOUT)
        --format = Currently only 'gff3' or the default 'text'.


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

    def exons(self):
        return self.children['exon']

    def CDSs(self):
        ## Why not "CDSes"?  As a grammarian, this gave me fits.  There are many references
        #  which suggest adding -es to any initialism, but in the end I had to go with Oxford's example:
        #  http://oxforddictionaries.com/definition/american_english/SOS
        return self.children['CDS']

    def introns(self, on=None):
        '''
        Dynamically generates Intron objects in order for the current RNA.  The coordinates of the
        generated introns depend on the object passed via the 'on' argument
        '''
        if on is None:
            raise Exception("ERROR: the introns() method requires a passed molecule using the 'on' argument")
        
        mol_on = on

        intron_objs = list()
        last_exon = None
        last_exon_loc = None

        for exon in self.exons():
            exon_loc = exon.location_on( mol_on )

            if last_exon is not None:
                intron_id = uuid.uuid4()
                intron = Intron( id=intron_id )
                intron.locate_on( target=mol_on, fmin=last_exon_loc.fmax, fmax=exon_loc.fmin, strand=exon_loc.strand )
                intron_objs.append( intron )

            last_exon = exon
            last_exon_loc = exon_loc

        return intron_objs
    
    def has_introns( self ):
        if len( self.exons() ) > 1:
            return True
        else:
            return False

        
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
##   Hopefully this isn't too un-pythonic.  These were to aide in abstraction within classes
#############################

def _initialize_type_list( children, feattype ):
    if children is None:
        children = dict()
    
    if feattype not in children:
        children[feattype] = list()

    return children


def _print_thing( thing, fh=None ):
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
