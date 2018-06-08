
import uuid

#from biocode import utils, gff, tbl
import biocode.utils
import biocode.gff
import biocode.tbl

'''
Warning: this module requires Python 3.2 or higher

Yes, there's a lot of abstraction that could happen here still but I'm trying to allow for
special consideration for each feature type and for calls to read more like biology than CS.

Each of the 'things' is intended to be biologically relevant, and I've included the definitions
for each from the Sequence Ontology.  Their class representations may not match this perfectly
yet, but that's the goal wherever possible.  Some terms, such as Organism, fall outside the scope
of the Sequence Ontology but I found them relevant (or convenient) to include here.

Documentation should be shifted from pydoc to:
http://sphinx-doc.org/

A set of abstract classes are defined first, then the other classes, all within alphabetic order
within each group.
'''



class LocatableThing:
    '''
    Base class for any Things that are locatable on other Things, such as genes on contigs.

    You can use direct comparisons between LocatableThings to compare coordinates, such as:

        if gene1 < gene2:
            # do stuff

    These might not be what you first expect.  They are comparisons of the entire ranges
    of each feature and have the following meanings.  Keep in mind that these are all
    based on fmin and fmax coordinates, so strandedness isn't taken into account:

        thing1 > thing2 : The entire thing1 is to the right of thing2
        thing1 < thing2 : The entire thing1 is to the left of thing2
        thing1 >= thing2: The thing1 overlaps thing2 on the right end
        thing1 <= thing2: The thing1 overlaps thing2 on the left end
        thing1 == thing2: The two things share the exact same coordinates

    Here 'right' means past the fmax of the first feature and left means past the fmin.
    For visual examples, see the documentation for the functions that implement each.

    These comparisons depend on both objects sharing locations on the same object.

    '''
    def __init__( self, locations=None ):
        self.locations = locations

        if self.locations is None:
            self.locations = list()

    def __lt__(self, other):
        return self.is_on_min_side_of(thing=other)

    def __le__(self, other):
        return self.overlaps_min_side_of( thing=other )

    def __eq__(self, other):
        return self.has_same_coordinates_as(thing=other)

    def __ne__(self, other):
        return not self.has_same_coordinates_as(thing=other)

    def __gt__(self, other):
        return self.is_on_max_side_of(thing=other)

    def __ge__(self, other):
        return self.overlaps_max_side_of( thing=other )

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

    def has_same_coordinates_as( self, on=None, thing=None, stop_tolerant=False ):
        '''
        Reports True/False depending on whether the calling object shares both start/stop coordinates
        with the passed 'thing' when considering locations 'on' the passed molecule biothing object.

        If you set the stop_tolerant argument=True, this will check the exact coordinates AND
        the coordinates with three bases off the end.  This obviously takes longer.
        '''
        other = thing

        for this_loc in self.locations:
            for other_loc in other.locations:
                # are the molecules the same?
                if this_loc.on.id == other_loc.on.id:
                    # if a specific molecule was requested, is this it?
                    if on is not None and on.id != this_loc.on.id:
                        continue

                    if this_loc.fmin == other_loc.fmin and this_loc.fmax == other_loc.fmax:
                        return True

                    # for those features with mixed stop codon inclusion
                    elif stop_tolerant == True:
                        if this_loc.strand == 1:
                            # this allows for either to be three base pairs away from the other in either direction
                            if this_loc.fmin == other_loc.fmin and abs(this_loc.fmax - other_loc.fmax) == 3:
                                return True
                        else:
                            if abs(this_loc.fmin - other_loc.fmin) == 3 and this_loc.fmax == other_loc.fmax:
                                return True

        # if we got here, there wasn't a match
        return False

    def is_on_max_side_of( self, thing=None, on=None ):
        '''
        Returns True/False depending on whether the calling LocatableThing object has a location
        relative to the passed thing that is completely past its fmax coordinate, like this:

            this :                                           fmin---------------fmax
            thing:            fmin-------------------fmax

        Checks both to make sure they're both located on the same molecule before comparing coordinates,
        ignoring strandedness.  If you pass the 'on' argument, locations will be checked against
        that specific object.
        '''
        for this_loc in self.locations:
            for other_loc in thing.locations:
                # are the molecules the same?
                if this_loc.on.id == other_loc.on.id:
                    # if a specific molecule was requested, is this it?
                    if on is not None and on.id != this_loc.on.id:
                        continue

                    if this_loc.fmin > other_loc.fmax:
                        return True

        # if we got here, there wasn't a match
        return False

    def is_on_min_side_of( self, thing=None, on=None ):
        '''
        Returns True/False depending on whether the calling LocatableThing object has a location
        relative to the passed thing that is completely before its fmin coordinate, like this:

            this :     fmin---------------fmax
            thing:                                fmin-------------------fmax

        Checks both to make sure they're both located on the same molecule before comparing coordinates,
        ignoring strandedness.  If you pass the 'on' argument, locations will be checked against
        that specific object.
        '''
        for this_loc in self.locations:
            for other_loc in thing.locations:
                # are the molecules the same?
                if this_loc.on.id == other_loc.on.id:
                    # if a specific molecule was requested, is this it?
                    if on is not None and on.id != this_loc.on.id:
                        continue

                    if this_loc.fmax < other_loc.fmin:
                        return True

        # if we got here, there wasn't a match
        return False

    def locate_on( self, target=None, fmin=None, fmin_partial=False, fmax=None, fmax_partial=False, phase=None, strand=None ):
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

        loc = Location(on=target, fmin=fmin, fmin_partial=fmin_partial, fmax=fmax, fmax_partial=fmax_partial, strand=strand, phase=phase)
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

    def overlap_size_with( self, thing ):
        '''
        Returns None if no overlap found or not found on same molecule, else returns the number
        of overlapping bases between the two features.

        Checks both to make sure they're both located on the same molecule before comparing coordinates,
        ignoring strandedness.

        Warning: Does not check child features, which may be significant.  If you want to compare
        the exons of a transcript with the segmented match_parts of an alignment, for example
        this would only compare their overall ranges (assuming you passed transcript and match features
        for the comparison.)
        '''
        locs = self.shared_molecule_locations_with(thing)

        if locs[0] is None:
            return None

        ref_location = locs[0]
        #print("DEBUG: ref location is {0}-{1}".format(ref_location.fmin, ref_location.fmax) )
        qry_location = locs[1]
        #print("DEBUG: qry location is {0}-{1}".format(qry_location.fmin, qry_location.fmax) )

        ## first catch any which don't overlap. this makes the other tests easier
        if ref_location.fmax <= qry_location.fmin or qry_location.fmax <= ref_location.fmin:
            return None

        # overlap on right end of ref, left of qry
        # ref:   -------------------
        # qry:            -------------------
        elif ref_location.fmax <= qry_location.fmax and ref_location.fmin < qry_location.fmin:
            return ref_location.fmax - qry_location.fmin

        # overlap on left end of ref, right of qry
        # ref:           -------------------
        # qry:   -------------------
        elif ref_location.fmin >= qry_location.fmin and ref_location.fmax > qry_location.fmax:
            return qry_location.fmax - ref_location.fmin

        # is the reference contained within the qry?
        # ref:      ------------
        # qry:   -------------------
        elif ref_location.fmin >= qry_location.fmin and ref_location.fmax <= qry_location.fmax:
            return ref_location.fmax - ref_location.fmin

        # is the query contained within the ref?
        # ref:   -------------------
        # qry:      -------------
        elif qry_location.fmin >= ref_location.fmin and qry_location.fmax <= ref_location.fmax:
            return qry_location.fmax - qry_location.fmin
        else:
            raise Exception("ERROR: uncaught orientation in overlap_size_with between {0} and {1}".format(self.id, thing.id))


    def overlaps_max_side_of( self, thing=None, on=None ):
        '''
        Returns True/False depending on whether the two LocatableThing objects have
        even a single overlapping Location on the same molecule, specifically with
        respect to the fmax coordinate of the passed object.  Like this:

            this :                         fmin---------------fmax
            thing:            fmin-------------------fmax

        Checks both to make sure they're both located on the same molecule before comparing coordinates,
        ignoring strandedness.  If you pass the 'on' argument, locations will be checked against
        that specific object.
        '''
        for this_loc in self.locations:
            for other_loc in thing.locations:
                # are the molecules the same?
                if this_loc.on.id == other_loc.on.id:
                    # if a specific molecule was requested, is this it?
                    if on is not None and on.id != this_loc.on.id:
                        continue

                    if this_loc.fmin < other_loc.fmax and this_loc.fmax > other_loc.fmax:
                        return True

        # if we got here, there wasn't a match
        return False

    def overlaps_min_side_of( self, thing=None, on=None ):
        '''
        Returns True/False depending on whether the two LocatableThing objects have
        even a single overlapping Location on the same molecule, specifically with
        respect to the fmin coordinate of the passed object.  Like this:

            this :    fmin---------------fmax
            thing:            fmin-------------------fmax

        Checks both to make sure they're both located on the same molecule before comparing coordinates,
        ignoring strandedness.  If you pass the 'on' argument, locations will be checked against
        that specific object.
        '''
        for this_loc in self.locations:
            for other_loc in thing.locations:
                # are the molecules the same?
                if this_loc.on.id == other_loc.on.id:
                    # if a specific molecule was requested, is this it?
                    if on is not None and on.id != this_loc.on.id:
                        continue

                    if this_loc.fmin < other_loc.fmin and this_loc.fmax > other_loc.fmin:
                        return True

        # if we got here, there wasn't a match
        return False


    def overlaps_with( self, thing ):
        '''
        Returns True/False depending on whether the two LocatableThing objects have
        even a single overlapping Location on the same molecule.

        Checks both to make sure they're both located on the same molecule before comparing coordinates,
        ignoring strandedness.
        '''
        common_location_found = False

        ## see if either are located on the same thing
        for ref_location in self.locations:
            for qry_location in thing.locations:
                if ref_location.on.id == qry_location.on.id:
                    common_location_found = True

                    ## if so, see if they have overlapping coordinates
                    ## it's easier to negate an overlap
                    if ref_location.fmax <= qry_location.fmin or qry_location.fmax <= ref_location.fmin:
                        return False

        if common_location_found == True:
            return True
        else:
            return False

    def shared_molecule_locations_with( self, thing ):
        ## see if either are located on the same thing
        for ref_location in self.locations:
            for qry_location in thing.locations:
                if ref_location.on.id == qry_location.on.id:
                    return [ref_location, qry_location]

        ## if we got this far, one wasn't found
        return [None, None]

    def update_location( self, on=None, fmin=None, fmax=None, strand=None, phase=None ):
        '''
        Allows you to update the location of any Thing that is already located on another one.
        You can independently pass fmin, fmax, strand, or phase (or any combination) to only 
        update those attributes.
        '''
        if on is None:
            raise Exception("The 'on' attribute must be passed to the update_location() method.")

        loc = self.location_on(on)
        
        if fmin is not None: loc.fmin = fmin
        if fmax is not None: loc.fmax = fmax
        if strand is not None: loc.strand = strand
        if phase is not None: loc.phase = phase

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

    The possible values for strand are 1, 0 and -1.  1 is for forward features, -1 for reverse, and
    0 is for unstranded features.

    The phase is REQUIRED for all CDS features.
    '''
    def __init__( self, on=None, fmin=None, fmin_partial=False, fmax=None, fmax_partial=False, phase=None, strand=None ):
        self.on = on
        self.fmin = fmin
        self.fmin_partial = fmin_partial
        self.fmax = fmax
        self.fmax_partial = fmax_partial
        self.strand = strand
        self.phase = 0 if phase is None else phase
        

class MoleculeSet:
    '''
    Base class for any Set of different molecule types, such as AssemblySet or PolypeptideSet.

    It is appropriate for use for all classes inheriting from LocatableThing.
    '''
    def __init__(self):
        # Nothing here yet
        pass

    def load_from_dict(self, thedict):
        for id in thedict:
            self.add(thedict[id])

    def write_fasta(self, fh=None, path=None):
        '''
        Writes the current set in FASTA format.  You can either pass the fh or path arguments.  If
        an open file handle already exists, fh is appropriate.  Instead, if you have just a path you
        want to be written to pass the 'path' argument instead.

        The header format in the FASTA entries depends on the type of elements in the set.  
        '''
        if path is not None:
            fh = open(path, 'wt')

        if self.__class__ == PolypeptideSet:
            molecules = self.polypeptides
        elif self.__class__ == AssemblySet:
            molecules = self.assemblies
        else:
            raise Exception("ERROR: writing FASTA not supported in MoleculeSets of this type: {0}".format(self.__class__))

        for molecule in molecules:
            if self.__class__ == PolypeptideSet:
                header = molecule.annotation_string()
            elif self.__class__ == AssemblySet:
                header = molecule.id
                
            fh.write(">{0}\n".format(header))
            fh.write("{0}\n".format(biocode.utils.wrapped_fasta(molecule.residues)))

        fh.close()
        
class Organism:
    '''
    This is a pretty high-level representation of an organism, and (currently) mostly serves as a
    container for Assembly objects.
    '''
    def __init__( self, id=None, assemblies=None, taxon_id=None, genus=None, species=None, strain=None ):
        self.id = id
        self.assemblies = assemblies
        self.taxon_id = taxon_id

        # convenience attributes
        self.genus = genus
        self.species = species
        self.strain = strain

        if self.assemblies is None:
            self.assemblies = list()

        def add_assembly( self, assembly ):
            self.assemblies.append( assembly )

        def common_name(self):
            if self.strain is None:
                return "{0} {1}".format(self.genus, self.species)
            else:
                return "{0} {1} ({2})".format(self.genus, self.species, self.strain)



class Assembly( LocatableThing ):
    '''
    SO definition (2013-05-22): "A region of the genome of known length that is composed by ordering and
    aligning two or more different regions."
    '''

    def __init__( self, id=None, locations=None, length=None, residues=None, is_circular=None, children=None ):
        super().__init__(locations)
        self.id = id
        self.length = length
        self.residues = residues
        self.is_circular = is_circular
        self.children = children

        if length is None and self.residues is not None and len(self.residues) > 0:
            self.length = len(self.residues)

        ## initialize any types needed
        self.children = _initialize_type_list(self.children, 'gene')

    def add_gene( self, gene ):
        self.children['gene'].append(gene)

    def genes(self):
        '''
        Returns all Gene objects which are children of this Assembly
        '''
        return self.children['gene']


class AssemblySet( MoleculeSet ):
    '''
    This is a convenience class for when operations need to be performed on a set of contigs/scaffolds
    such as N50 calculation or writing to a FASTA file.
    '''
    def __init__( self, assemblies=None ):
        super().__init__()
        self.assemblies = assemblies

        if self.assemblies is None:
            self.assemblies = list()

    def add( self, assembly ):
        self.assemblies.append(assembly)

    def load_from_file(self, file):
        seqs = biocode.utils.fasta_dict_from_file(file)
        
        for seq_id in seqs:
            assembly = Assembly(id=seq_id, residues=seqs[seq_id]['s'])
            self.add(assembly)

    def N50( self ):
        '''
        N50 is a statistical measure of average length of a set of sequences. It is used widely in
        genomics, especially in reference to contig or supercontig lengths within a draft assembly.

        Given a set of sequences of varying lengths, the N50 length is defined as the length N for
        which 50% of all bases in the sequences are in a sequence of length L < N. This can be found
        mathematically as follows: Take a list L of positive integers. Create another list L' ,
        which is identical to L, except that every element n in L has been replaced with n copies of
        itself. Then the median of L' is the N50 of L. For example: If L = {2, 2, 2, 3, 3, 4, 8, 8},
        then L' consists of six 2's, six 3's, four 4's, and sixteen 8's; the N50 of L is the median
        of L' , which is 6.
        
        Source: http://www.broadinstitute.org/crd/wiki/index.php/N50
        '''
        total_len = 0
        lengths = list()

        for assembly in self.assemblies:
            total_len += assembly.length
            lengths.append(assembly.length)

        lengths = sorted(lengths)

        count = 0

        for i in range(len(lengths)):
            count += lengths[i]
            if count >= total_len / 2:
                return lengths[i]


class CDS( LocatableThing ):
    '''
    SO definition (2013-05-22): "A contiguous sequence which begins with, and includes, a start codon
    and ends with, and includes, a stop codon."

    Note that there is a minor variation from this, in that there is no enforcement or expectation
    that the sequence coordinates or residues stored begin with a start or end with a stop.  This is
    so that we can still use this object to represent CDS within partial gene predictions, which may
    or may not have complete ends.
    '''
    def __init__( self, id=None, locations=None, parent=None, phase=None, length=None, residues=None, annotation=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.length = length
        self.residues = residues
        self.phase = 0 if phase is None else phase

        ## this should be an instance of FunctionalAnnotation from annotation.py
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

        self.residues = mol.residues[loc.fmin:loc.fmax]
        self.length = len(self.residues)

        if loc.strand == -1:
            self.residues = biocode.utils.reverse_complement(self.residues)

        return self.residues


class Exon( LocatableThing ):
    def __init__( self, id=None, locations=None, parent=None, length=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.length = length


class Gene( LocatableThing ):
    '''
    SO definition (2013-05-22): "A region (or regions) that includes all of the sequence
    elements necessary to encode a functional transcript. A gene may include regulatory
    regions, transcribed regions and/or other functional sequence regions."
    '''
    def __init__( self, id=None, locus_tag=None, locations=None, children=None, residues=None ):
        super().__init__(locations)
        self.id = id
        self.locus_tag = locus_tag
        self.children = children
        self.residues = residues

        ## initialize any types needed
        self.children = _initialize_type_list(self.children, 'mRNA')
        self.children = _initialize_type_list(self.children, 'rRNA')
        self.children = _initialize_type_list(self.children, 'tRNA')
        self.children = _initialize_type_list(self.children, 'tmRNA')

    def __hash__(self):
        return hash(self.id)

    def add_RNA(self, rna):
        """
        Utility method which can be used to add any RNA subtype (mRNA, rRNA, etc.)
        """
        if isinstance(rna, mRNA):
            self.add_mRNA(rna)
        elif isinstance(rna, rRNA):
            self.add_rRNA(rna)
        elif isinstance(rna, tRNA):
            self.add_tRNA(rna)
        elif isinstance(rna, tmRNA):
            self.add_tmRNA(rna)
        else:
            raise Exception("ERROR: add_RNA() method called for unrecognized type: {0}".format(rna.__class__.__name__))

    def add_mRNA(self, rna):
        rna.parent = self
        self.children['mRNA'].append(rna)

    def add_rRNA(self, rna):
        rna.parent = self
        self.children['rRNA'].append(rna)

    def add_tRNA(self, rna):
        rna.parent = self
        self.children['tRNA'].append(rna)

    def add_tmRNA(self, rna):
        rna.parent = self
        self.children['tmRNA'].append(rna)

    def get_residues(self):
        if len(self.locations) == 0:
            raise Exception("ERROR: gene.get_residues() requested but gene {0} isn't located on anything.".format(self.id))
        elif len(self.locations) > 1:
            raise Exception("ERROR: gene {0} is located on multiple molecules.  Can't automatically extract the residues.".format(self.id))

        loc = self.location()
        mol = loc.on

        # make sure this thing has its residues populated
        if len(mol.residues) <= 0:
            raise Exception("ERROR: gene.get_residues() requested but its molecule {0} has no stored residues".format(mol.id))

        self.residues = mol.residues[loc.fmin:loc.fmax]
        self.length = len(self.residues)

        if loc.strand == -1:
            self.residues = biocode.utils.reverse_complement(self.residues)

        return self.residues

    def mRNA_count(self):
        return len(self.mRNAs())

    def mRNAs(self):
        return self.children['mRNA']

    def ncRNAs(self):
        return self.children['tRNA'] + self.children['rRNA'] + self.children['tmRNA']
    
    def polypeptides(self):
        """
        Returns a list() of all Polypeptide objects for all mRNA children.
        """
        polypeptides = list()

        for rna in self.RNAs():
            polypeptides.extend(rna.polypeptides())

        return polypeptides
        

    def RNAs(self):
        return self.children['mRNA'] + self.children['tRNA'] + self.children['rRNA'] + self.children['tmRNA']

    def remove_mRNA(self, rna):
        if rna in self.children['mRNA']:
            self.children['mRNA'].remove(rna)
        
    def rRNAs(self):
        return self.children['rRNA']

    def shares_exon_structure_with( self, thing=None, stop_tolerant=False ):
        """
        This checks if two genes have only one mRNA and, if so, compares their internal
        exon structure.  Returns True if they all match completely.
        """
        these_mRNAs = self.mRNAs()
        other_mRNAs = thing.mRNAs()

        if len(these_mRNAs) != 1 or len(other_mRNAs) != 1 or len(these_mRNAs[0].exons()) != len(other_mRNAs[0].exons()):
            return False;

        exon_matches_found = 0
        ref_exon_count = 0

        for ref_exon in these_mRNAs[0].exons():
            ref_exon_count += 1

            for other_exon in other_mRNAs[0].exons():
                if ref_exon.has_same_coordinates_as( thing=other_exon, stop_tolerant=stop_tolerant ):
                    exon_matches_found += 1
                    break

        if exon_matches_found == ref_exon_count:
            return True
        else:
            return False

    def shares_CDS_structure_with( self, thing=None ):
        """
        This checks if two genes have only one mRNA and, if so, compares their internal
        CDS structure.  Returns True if they all match completely.
        """
        these_mRNAs = self.mRNAs()
        other_mRNAs = thing.mRNAs()

        if len(these_mRNAs) != 1 or len(other_mRNAs) != 1 or len(these_mRNAs[0].CDSs()) != len(other_mRNAs[0].CDSs()):
            return False;

        CDS_matches_found = 0
        ref_CDS_count = 0

        for ref_CDS in these_mRNAs[0].CDSs():
            ref_CDS_count += 1

            for other_CDS in other_mRNAs[0].CDSs():
                if ref_CDS.has_same_coordinates_as( thing=other_CDS ):
                    CDS_matches_found += 1
                    break

        if CDS_matches_found == ref_CDS_count:
            return True
        else:
            return False

    def tRNAs(self):
        return self.children['tRNA']

    def print_as(self, fh=None, source=None, format=None, lab=None):
        '''
        Prints the gene entry in a variety of available formats, with all children.
        Is printed to STDOUT unless the fh option is passed.  Really, only checks are
        done here before passing the info to the biocodegff module.

        The options passed will decide the target and format of the print.

            --fh = File handle to which we should print (default: STDOUT)
            --format = Currently only 'gff3', 'tbl' or the default 'text'.
            --source = Required for gff printing
            --lab = Required for tbl printing
        ''' 
        if format == 'text':
            _print_thing(self, fh=fh)
        elif format == 'gff3':
            biocode.gff.print_biogene(gene=self, fh=fh, source=source)
        elif format == 'tbl':
            biocode.tbl.print_biogene(gene=self, fh=fh, lab=lab)
        else:
            raise Exception("ERROR: the print_as method only accepts values of 'text' or 'gff3'")

    def remove_mRNA(self, rna):
        if rna in self.children['mRNA']:
            self.children['mRNA'].remove(rna)
        else:
            raise Exception("ERROR: Can't remove mRNA {0} from gene {1} because it wasn't found".format(rna.id, self.id))


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


class Match( LocatableThing ):
    '''
    Represents a generic match object, usually used to represent pairwise alignments.
    The sequence ontology provides a host of specific match terms, but I'm using this
    as a lightweight implementation of them all rather than explicitly creating
    dozens of subclasses.  Instead, create a Match object and simply set any more
    specific term using Match.class

    Future mod: Add Annotation object to store info on the match?
    '''

    def __init__( self, id=None, subclass=None, target_id=None, locations=None, parts=None, length=None ):
        super().__init__(locations)
        self.id = id
        self.subclass = subclass   # allows you to define specific classes (cDNA_match, etc.)
        self.target_id = target_id
        self.parts = parts   # these are MatchPart objects
        self.length = length

        if self.subclass is None:
            self.subclass = 'match'

        if self.parts is None:
            self.parts = list()

    #def __repr__(self):
    #    '''
    #    Prints a string representation of the match along with it's match_parts
    #    '''
    #    lines = list()
    #    lines.append("Match - id:{0}, subclass:{1}, target_id:{2}, length:{3}".format( \
    #                  self.id, self.subclass, self.target_id, self.length))
    #    lines.append("\tlocations: {0}".format(len(self.locations)) )
    #    lines.append("\tmatch_parts: {0}".format(len(self.parts)) )

    #    for mp in self.parts:
    #        print(mp)
    #        lines.append("\t\tMatchPart - id:{0}, parent_id:{1}, length{2}".format( \
    #                  mp.id, mp.parent.id, mp.length))
    #        lines.append("\t\t\tlocations: {0}".format(len(self.locations)) )

    #    return "\n".join(lines)

    def add_part(self, part):
        self.parts.append(part)

    def print_as(self, fh=None, source=None, format=None):
        if format == 'text':
            _print_thing(self, fh=fh)
        elif format == 'gff3':
            # mode could be passed here to print both GFF3-supported representations (already implemented)
            biocode.gff.print_biomatch(match=self, fh=fh, source=source)
        else:
            raise Exception("ERROR: You attempted to print a Match with unrecognized format:{0}".format(format))


class MatchPart( LocatableThing ):
    """ Future mod: Add Annotation object to store info on the match? """
    def __init__( self, id=None, locations=None, parent=None, length=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent  # this should be a Match object
        self.length = length

    #def __repr__(self):
    #    '''
    #    Prints a string representation of the match_part along with it's locations
    #    '''
    #    lines = list()
    #    lines.append("MatchPart - id:{0}, parent_id:{1}, length{2}".format( \
    #                  self.id, self.parent.id, self.length))
    #    lines.append("\tlocations: {0}".format(len(self.locations)) )

    #    return "\n".join(lines)

    
class Polypeptide( LocatableThing ):
    '''
    SO definition (2013-05-22): "A sequence of amino acids linked by peptide bonds which may lack
    appreciable tertiary structure and may not be liable to irreversible denaturation."
    '''
    def __init__( self, id=None, locations=None, parent=None, length=None, annotation=None, residues=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.length = length
        self.residues = residues

        ## this should be an instance of FunctionalAnnotation from annotation.py
        self.annotation = annotation

    def annotation_string(self):
        '''
        Returns a string of the annotation attached to this Polypeptide, most often used as
        a FASTA entry header
        '''
        header = ""
        go_string = ""
        ec_string = ""

        if self.annotation is None:
            return self.id
        else:
            for go_annot in self.annotation.go_annotations:
                go_string += "GO:{0},".format(go_annot.go_id)

            go_string = go_string.rstrip(',')

            for ec_annot in self.annotation.ec_numbers:
                ec_string += "{0},".format(ec_annot.number)

            ec_string = ec_string.rstrip(',')

            header = "{0} {1}".format(self.id, self.annotation.product_name)

            if self.annotation.gene_symbol is not None:
                header = "{0} gene::{1}".format(header, self.annotation.gene_symbol)

            if ec_string != "":
                header = "{0} ec::{1}".format(header, ec_string)

            if go_string != "":
                header = "{0} go::{1}".format(header, go_string)

        return header

    
class PolypeptideSet( MoleculeSet ):
    '''
    This is a convenience class for when operations need to be performed on a set of polypeptides
    such as writing to a FASTA file.
    '''
    def __init__( self, polypeptides=None ):
        super().__init__()
        self.polypeptides = polypeptides

        if self.polypeptides is None:
            self.polypeptides = list()

    def add( self, polypeptide ):
        self.polypeptides.append( polypeptide )

    def load_from_file(self, file):
        seqs = biocode.utils.fasta_dict_from_file(file)
        
        for seq_id in seqs:
            polypeptide = Polypeptide(id=seq_id, residues=seqs[seq_id]['s'])
            self.add(polypeptide)
        
class RNA( LocatableThing ):
    '''
    SO definition (2013-05-22): "An attribute describing a sequence consisting of nucleobases
    bound to a repeating unit made of a D-ribose ring connected to a phosphate backbone."
    '''
    def __init__( self, id=None, locations=None, parent=None, locus_tag=None, children=None, annotation=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.locus_tag = locus_tag
        self.children = children

        ## This should be an instance of FunctionalAnnotation from annotation.py
        #   It's considered best practice to put the annotation on the Polypeptide feature for coding RNAs.
        #   Here we expect things like product (rRNAs), anticodon (tRNAs), etc.
        self.annotation = annotation

        ## initialize any types needed
        self.children = _initialize_type_list(self.children, 'exon')
        self.children = _initialize_type_list(self.children, 'CDS')
        self.children = _initialize_type_list(self.children, 'polypeptide')
        self.children = _initialize_type_list(self.children, 'UTR')

    def __hash__(self):
        return hash(self.id)

    def add_CDS(self, cds):
        cds.parent = self
        self.children['CDS'].append(cds)

    def add_exon(self, exon):
        exon.parent = self
        self.children['exon'].append(exon)

    def add_five_prime_UTR(self, utr):
        utr.parent = self
        self.children['UTR'].append(utr)

    def add_polypeptide(self, polypeptide):
        polypeptide.parent = self
        self.children['polypeptide'].append(polypeptide)

    def add_three_prime_UTR(self, utr):
        utr.parent = self
        self.children['UTR'].append(utr)

    def add_UTR(self, utr):
        if isinstance(utr, FivePrimeUTR):
            self.add_five_prime_UTR(utr)
        elif isinstance(utr, ThreePrimeUTR):
            self.add_three_prime_UTR(utr)
        else:
            raise Exception("ERROR: add_UTR() method called for unrecognized type: {0}".format(utr.__class__.__name__))

    def CDS_count(self):
        return len(self.CDSs())

    def CDSs(self):
        ## Why not "CDSes"?  As a grammarian, this gave me fits.  There are many references
        #  which suggest adding -es to any initialism, but in the end I had to go with Oxford's example:
        #  http://oxforddictionaries.com/definition/american_english/SOS
        #   -- J. Orvis
        #
        # This method was added for consistency of usage and naming, but it really returns an array
        #  of the CDS fragments in spliced genes.  If you want the contiguous sequence, use get_CDS_residues()
        #
        # Note that this is not guaranteed to return the CDSs in sorted chromosomal order
        return self.children['CDS']

    def delete_CDS(self, cds_to_remove):
        idx = 0

        for cds in self.children['CDS']:
            if cds.id == cds_to_remove.id:
                del self.children['CDS'][idx]
                break

            idx += 1
        else:
            raise Exception("ERROR: Attempt to delete a CDS ({0}) which wasn't found as a child of this RNA ({1}".format(cds.id, self.id))

    def exon_count(self):
        return len(self.exons())

    def exons(self):
        return self.children['exon']

    def extend_stop(self, on=None, to=None):
        '''
        Extending an RNA involves changing its coordinates, checking if the extension requires
        the parent gene to be extended, and checking sub features (exon, CDS) for extension.
        This does all of those things with respect to the molecule passed via the 'on' argument.

        on = Molecule object this RNA is located on
        to = To what coordinate should the stop be extended?  The strand will be checked
             automatically
        '''
        if on is None or to is None:
            raise Exception("The 'on' and 'to' attributes must be passed to the extend_stop() method.")

        if to < 0:
            raise Exception("Attempt to extend coordinates for {0} on {1} to a negative coordinate ({2})".format(self.id, on.id, to))

        # Do we need to extend the gene?
        gene = self.parent
        gene_loc = gene.location_on(on)
        self_loc = self.location_on(on)

        if self_loc.strand == 1:
            # We don't always need to extend the RNA.  It may just be the internal structure that needs it.
            if self_loc.fmax < to:
                self_loc.fmax = to

            if gene_loc.fmax < to:
                gene_loc.fmax = to

            if self.CDS_count() > 0:
                terminal_CDS = sorted(self.CDSs())[-1]
                terminal_CDS_loc = terminal_CDS.location_on(on)
                terminal_CDS_loc.fmax = to

            if self.exon_count() > 0:
                terminal_exon = sorted(self.exons())[-1]
                terminal_exon_loc = terminal_exon.location_on(on)
                terminal_exon_loc.fmax = to
        else:
            if self_loc.fmin > to:
                self_loc.fmin = to                

            if gene_loc.fmin > to:
                gene_loc.fmin = to

            if self.CDS_count() > 0:
                terminal_CDS = sorted(self.CDSs())[0]
                terminal_CDS_loc = terminal_CDS.location_on(on)
                terminal_CDS_loc.fmin = to

            if self.exon_count() > 0:
                terminal_exon = sorted(self.exons())[0]
                terminal_exon_loc = terminal_exon.location_on(on)
                terminal_exon_loc.fmin = to
        
    
    def five_prime_UTRs(self):
        utrs = list()

        for utr in self.UTRs():
            if isinstance(utr, FivePrimeUTR):
                utrs.append(utr)

        return utrs

    def get_CDS_residues(self, for_translation=False):
        '''
        Returns a DNA string representing the CDS based on the coordinate spans.  Note that phase
        will NOT be respected unless you pass for_translation=True
        '''
        if len(self.locations) == 0:
            raise Exception("ERROR: RNA.get_CDS_residues() requested but RNA {0} isn't located on anything.".format(self.id))
        elif len(self.locations) > 1:
            raise Exception("ERROR: RNA {0} is located on multiple molecules.  Can't automatically extract the residues.".format(self.id))

        loc = self.location()
        segments = list()

        chop = None

        # these MUST be handled in a sorted order
        residues = ''
        sorted_cds = sorted(self.CDSs())
        
        if len(sorted_cds) > 0:
            for cds in sorted_cds:
                cds_segment = cds.get_residues()
                segments.append(cds_segment)

            if loc.strand == -1:
                segments = reversed(segments)
                chop = sorted_cds[-1].phase
            else:
                chop = sorted_cds[0].phase

            for segment in segments:
                residues += segment

            if chop is not None:
                residues = residues[chop:]

        return residues

    def has_introns( self ):
        if len( self.exons() ) > 1:
            return True
        else:
            return False

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

        for exon in sorted(self.exons()):
            exon_loc = exon.location_on( mol_on )

            if last_exon is not None:
                intron_id = uuid.uuid4()
                intron = Intron( id=intron_id )
                intron.locate_on( target=mol_on, fmin=last_exon_loc.fmax, fmax=exon_loc.fmin, strand=exon_loc.strand )
                intron_objs.append( intron )

            last_exon = exon
            last_exon_loc = exon_loc

        return intron_objs

    def polypeptides(self):
        return self.children['polypeptide']

    def three_prime_UTRs(self):
        utrs = list()

        for utr in self.UTRs():
            if isinstance(utr, ThreePrimeUTR):
                utrs.append(utr)

        return utrs

    def UTRs(self):
        return self.children['UTR']

class mRNA( RNA ):
    '''
    SO definition (2013-05-22): "Messenger RNA is the intermediate molecule between DNA and protein. It
    includes UTR and coding sequences. It does not contain introns."
    '''
    def __init__( self, id=None, locations=None, parent=None, locus_tag=None, children=None, residues=None ):
        super().__init__(id, locations, parent, locus_tag, children)
        self.residues = residues

class mRNASet( MoleculeSet ):
    '''
    This is a convenience class for when operations need to be performed on a set of mRNAs
    such as writing to a FASTA file.
    '''
    def __init__( self, mRNAs=None ):
        super().__init__()
        self.mRNAs = mRNAs

        if self.mRNAs is None:
            self.mRNAs = list()

    def add( self, mRNA ):
        self.mRNAs.append( mRNA )

    def load_from_file(self, file):
        seqs = biocode.utils.fasta_dict_from_file(file)
        
        for seq_id in seqs:
            mRNA = mRNA(id=seq_id, residues=seqs[seq_id]['s'])
            self.add(mRNA)

class tmRNA( RNA ):
    '''
    SO definition (2017-11-29): "A tmRNA liberates a mRNA from a stalled ribosome. To accomplish this part 
    of the tmRNA is used as a reading frame that ends in a translation stop signal. The broken mRNA is 
    replaced in the ribosome by the tmRNA and translation of the tmRNA leads to addition of a proteolysis 
    tag to the incomplete protein enabling recognition by a protease. Recently a number of permuted tmRNAs 
    genes have been found encoded in two parts. TmRNAs have been identified in eubacteria and some 
    chloroplasts but are absent from archeal and Eukaryote nuclear genomes."
    '''
    def __init__( self, id=None, locations=None, parent=None, locus_tag=None, children=None ):
        super().__init__(id, locations, parent, locus_tag, children)
            
class rRNA( RNA ):
    '''
    SO definition (2013-05-22): "RNA that comprises part of a ribosome, and that can provide both
    structural scaffolding and catalytic activity."
    '''
    def __init__( self, id=None, locations=None, parent=None, locus_tag=None, children=None ):
        super().__init__(id, locations, parent, locus_tag, children)

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
    def __init__( self, id=None, locations=None, parent=None, locus_tag=None, children=None, anticodon=None ):
        super().__init__(id, locations, parent, locus_tag, children)
        self.anticodon = anticodon

class UTR( LocatableThing ):
    '''
    SO definition (2014-06-12): "Messenger RNA sequences that are untranslated and lie five prime or three
    prime to sequences which are translated."
    '''
    def __init__( self, id=None, locations=None, parent=None, length=None ):
        super().__init__(locations)
        self.id = id
        self.parent = parent
        self.length = length


class FivePrimeUTR( UTR ):
    '''
    SO definition (2014-06-12): "A region at the 5' end of a mature transcript (preceding the initiation
    codon) that is not translated into a protein."
    '''
    def __init__( self, id=None, locations=None, parent=None, length=None ):
        super().__init__(id, locations, parent, length)


class ThreePrimeUTR( UTR ):
    '''
    SO definition (2014-06-12): "A region at the 3' end of a mature transcript (following the stop codon)
    that is not translated into a protein."
    '''
    def __init__( self, id=None, locations=None, parent=None, length=None ):
        super().__init__(id, locations, parent, length)


#############################
## Private help functions
##   Hopefully this isn't too un-pythonic.  These were to aid in abstraction within classes
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
