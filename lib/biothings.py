import uuid
import biocodegff
import sys

'''
Warning: this module requires Python 3.2 or higher

Yes, there's a lot of abstraction that could happen here still but I'm trying to allow for 
special consideration to be allowed for each feature type and for calls to read more like
biology than CS.

Each of the 'things' is intended to be biologically relevant, and I've included the definitions
for each from the Sequence Ontology.  Their class representations may not match this perfectly
yet, but that's the goal wherever possible.
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


class Location:
    '''
    Notes on the 'phase' attribute, used as defined in the GFF3 spec:

    For features of type "CDS", the phase indicates where the feature begins with reference to the 
    reading frame. The phase is one of the integers 0, 1, or 2, indicating the number of bases that 
    should be removed from the beginning of this feature to reach the first base of the next codon. 
    In other words, a phase of "0" indicates that the next codon begins at the first base of the region 
    described by the current line, a phase of "1" indicates that the next codon begins at the second 
    base of this region, and a phase of "2" indicates that the codon begins at the third base of this 
    region. This is NOT to be confused with the frame, which is simply start modulo 3.
    
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
    '''
    SO definition (2013-05-22): "A region of the genome of known length that is composed by ordering and 
    aligning two or more different regions."
    '''
    
    def __init__( self, id=None, locations=None, length=None, residues=None, children=None ):
        super().__init__(locations)
        self.id = id
        self.length = length
        self.residues = residues
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
    '''
    SO definition (2013-05-22): "A contiguous sequence which begins with, and includes, a start codon 
    and ends with, and includes, a stop codon."
    '''
    def __init__( self, id=None, locations=None, parent=None, phase=None, length=None, residues=None, annotation=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.length = length
        self.residues = residues
        self.phase = 0 if phase is None else phase

        ## this should be an instance of FunctionalAnnotation from bioannotation.py
        self.annotation = annotation

    def get_residues(self):
        if len(self.locations) == 0:
            raise Exception("ERROR: CDS.get_residues() requested but CDS {0} isn't located on anything.".format(self.id))
        elif len(self.locations) > 1:
            raise Exception("ERROR: CDS {0} is located on multiple molecules.  Can't automatically extract the residues.".format(self.id))

        loc = self.location()
        mol = loc.on

        # make sure this thing has its residues populated
        if len(mol.residues) <= 0:
            raise Exception("ERROR: CDS.get_residues() requested but its molecule {0} has no stored residues".format(mol.id))

        print("extracting CDS at {0}-{1}".format(loc.fmin, loc.fmax))
        self.residues = mol.residues[loc.fmin:loc.fmax]
        self.length = len(self.residues)

        return self.residues


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
    
    SO definition (2013-05-22): "A region of a primary transcript that is transcribed, but 
    removed from within the transcript by splicing together the sequences (exons) on either 
    side of it."
    '''
    def __init__( self, id=None, locations=None, length=None ):
        super().__init__(locations)
        self.id = id
        self.length = length
        

class Gene( LocatableThing ):
    '''
    SO definition (2013-05-22): "A region (or regions) that includes all of the sequence 
    elements necessary to encode a functional transcript. A gene may include regulatory 
    regions, transcribed regions and/or other functional sequence regions."
    '''
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
    '''
    def print_as(self, fh=None, source=None, format=None):
        if format == 'text':
            _print_thing(self, fh=fh)
        elif format == 'gff3':
            biocodegff.print_biogene( gene=self, fh=fh, source=source )
        else:
            raise Exception("ERROR: the print_as method only accepts values of 'text' or 'gff3'")
        

class Polypeptide( LocatableThing ):
    '''
    SO definition (2013-05-22): "A sequence of amino acids linked by peptide bonds which may lack 
    appreciable tertiary structure and may not be liable to irreversible denaturation."
    '''
    def __init__( self, id=None, locations=None, parent=None, length=None, annotation=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.length = length
        
        ## this should be an instance of FunctionalAnnotation from bioannotation.py
        self.annotation = annotation
        
        
class RNA( LocatableThing ):
    '''
    SO definition (2013-05-22): "An attribute describing a sequence consisting of nucleobases 
    bound to a repeating unit made of a D-ribose ring connected to a phosphate backbone."
    '''
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
    '''
    SO definition (2013-05-22): "Messenger RNA is the intermediate molecule between DNA and protein. It 
    includes UTR and coding sequences. It does not contain introns."
    '''
    def __init__( self, id=None, locations=None, parent=None, children=None ):
        super().__init__(id, locations, parent, children)

class rRNA( RNA ):
    '''
    SO definition (2013-05-22): "RNA that comprises part of a ribosome, and that can provide both 
    structural scaffolding and catalytic activity."
    '''
    def __init__( self, id=None, locations=None, parent=None, children=None ):
        super().__init__(id, locations, parent, children)

class tRNA( RNA ):
    '''
    SO definition (2013-05-22): "Transfer RNA (tRNA) molecules are approximately 80 nucleotides in length. 
    Their secondary structure includes four short double-helical elements and three loops (D, anti-codon, 
    and T loops). Further hydrogen bonds mediate the characteristic L-shaped molecular structure. Transfer 
    RNAs have two regions of fundamental functional importance: the anti-codon, which is responsible for 
    specific mRNA codon recognition, and the 3' end, to which the tRNA's corresponding amino acid is 
    attached (by aminoacyl-tRNA synthetases). Transfer RNAs cope with the degeneracy of the genetic code 
    in two manners: having more than one tRNA (with a specific anti-codon) for a particular amino acid; 
    and 'wobble' base-pairing, i.e. permitting non-standard base-pairing at the 3rd anti-codon position."
    '''
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
        child_count = 0

        for child_type in thing.children:
            for child in thing.children[child_type]:
                child_count += 1
                print("\t\t{0}".format(child.id) )
        
        print("\tchild count = {0}".format( child_count ) )
    else:
        print("\tchildren = None")

    if hasattr(thing, 'locations'):
        print("\tlocation count: {0}".format(len(thing.locations)) )
        print("\tlocations:")

        for loc in thing.locations:
            print("\t\t{0}\t{1}\t{2}\t{3}".format(loc.on.id, loc.fmin, loc.fmax, loc.strand) )
    else:
        print("\tlocations = None")
