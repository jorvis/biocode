'''
Warning: this module requires Python 3.2 or higher

This is a set of classes to represent what I have most commonly needed when working on
dozens (eukaryotic) or hundreds (prokaryotic) of annotation projects.  If you have
other attributes that you'd like to see supported, please add an 'issue' on the
biocode GitHub page.
'''

import re

class FunctionalAnnotation:
    """
    While recognizing that an enormous variety of attributes could go here in
    describing the functional annotation of a BioThing, I'm starting with those
    we most often encounter and need to be available in common output formats.
    These most common attributes are accessed by name, but the others are found
    as a dict in the 'other_attributes' property.

    Also, there's a place for having attributes like this abstracted, stored in
    ontologies, etc.  We've done all that before.  For now I'm going to try
    and hopefully enjoy the utility of having the most common properties always
    directly, and simply available.
    """
    def __init__( self, product_name=None, gene_symbol=None, go_annotations=None, ec_numbers=None, dbxrefs=None ):
        self.product_name     = product_name
        self.gene_symbol      = gene_symbol
        self.go_annotations   = go_annotations
        self.ec_numbers       = ec_numbers
        self.dbxrefs          = dbxrefs
        self.other_attributes = dict()

        if self.go_annotations is None:
            self.go_annotations = list()

        if self.ec_numbers is None:
            self.ec_numbers = list()

        if self.dbxrefs is None:
            self.dbxrefs = list()

    def __str__(self):
        representation = "Product name: {0}\nGene symbol : {1}\n".format(self.product_name, self.gene_symbol)

        if len(self.go_annotations) > 0:
            representation += "GO annotations:\n"
            for go_annot in self.go_annotations:
                representation += "\tGO:{0}\n".format(go_annot.go_id)
        else:
            representation += "GO annotations: None\n"

        if len(self.ec_numbers) > 0:
            representation += "EC numbers:\n"
            for ec in self.ec_numbers:
                representation += "\t{0}\n".format(ec.number)
        else:
            representation += "EC numbers: None\n"

        if len(self.dbxrefs) > 0:
            representation += "Dbxrefs:\n"
            for dbxref in self.dbxrefs:
                representation += "\t{0}:{1}\n".format(dbxref.db, dbxref.identifier)
        else:
            representation += "Dbxrefs: None\n"

        return representation
                         
    def add_dbxref(self, dbxref):
        """
        Stores a Dbxref object within an annotation. The thing passed can either be a
        Dbxref object or string like "sourcedb:identifier" and it will be automatically
        parsed.
        """
        if type(dbxref).__name__ == 'Dbxref':
            self.dbxrefs.append(dbxref)
        elif type(dbxref).__name__ == 'str':
            m = re.match("(.+)\:(.+)", dbxref)
            if m:
                self.dbxrefs.append(Dbxref(db=m.group(1), identifier=m.group(2)))
            else:
                raise Exception("ERROR: Annotation.add_dbxref(): If string passed, expected format was 'source:identifier'")
        else:
            raise Exception("ERROR: Annotation.add_dbxref expected a Dbxref object or string to be passed")
        
    def add_ec_number(self, ec_num):
        """
        Note to self: Modify this to allow passing ECAnnotation object or string.
        Right now it expects an ECAnnotation object
        """
        self.ec_numbers.append(ec_num)

    def add_go_annotation(self, go):
        """
        Note to self: Modify this to allow passing GOAnnotation object or string.
        Right now it expects an GOAnnotation object
        """
        self.go_annotations.append(go)

class Dbxref:
    """
    These allow for specification of an identifier in another database by ID.  These
    are not meant to be used for Ontology term linking, but rather limited only to
    identifiers.  Examples:
    
       SGD:S0006169
       KEGG:K06223

    The first part is the 'db' and the second is the 'identifier'.

    ## notes on method signature overloading (for future use)
    jorvis: You'd write a @classmethod. I'd name it "fromstring"
    jorvis: the classmethod would parse that string into values which it would use to
            invoke the normal constructor, as appropriate.
    """
    def __init__( self, db=None, identifier=None ):
        self.db = db
        self.identifier = identifier


class GOAnnotation:
    """
    A functional annotation can have an infinite number of associated GO Annotations

    Details here:
    http://www.geneontology.org/GO.evidence.shtml

    Yes, the 'with_from' attribute name is awkward, but 'with/from' isn't legal and
    both 'with' and 'from' are python reserved words.

    The go_id attribute here is just the numeric portion without "GO" or "GO:" or
    anything else attached (allowing the developer to define it as required.)
    """
    def __init__( self, go_id=None, ev_code=None, with_from=None ):
        self.go_id        = go_id
        self.ev_code      = ev_code
        self.with_from    = with_from

        ## process any GO ID passed to only contain the numeric portion
        go_pattern = re.compile('(\d+)')
        m = go_pattern.search(self.go_id)

        if m:
            self.go_id = m.group(1)
        else:
            raise Exception("ERROR: failed to extract numeric portion of ID from new GOAnnotation")


class ECAnnotation:
    """
    A functional annotation can have an infinite number of associated EC Annotations

    Details here:
    http://www.chem.qmul.ac.uk/iubmb/enzyme/

    While the official terms for the levels are 'class', 'subclass' etc. we have to use
    different keys for the attributes since those conflict with reserved keywords in
    both python and other frameworks.

    class1 = 1          = Oxidoreductases
    class2 = 1.10       = Acting on diphenols and related substances as donors
    class3 = 1.10.3     = With oxygen as acceptor
    number = 1.10.3.2   = laccase

    Currently does not have an index of EC terms to provide other attributes which will
    be added in the future, such as:

    accepted_name = laccase
    reaction = 4 benzenediol + O2 = 4 benzosemiquinone + 2 H2O
    systematic_name = benzenediol:oxygen oxidoreductase
    CAS_registry_number = 80498-15-3
    """
    def __init__( self, number=None ):
        self.number = number
        self.class1 = None
        self.class2 = None
        self.class3 = None

        re_pattern = re.compile('(((([0-9\-]+)\.[0-9\-]+)\.[0-9\-]+)\.[a-z0-9\-]+)')
        m = re_pattern.search(self.number)

        if m:
            self.class1 = m.group(4)
            self.class2 = m.group(3)
            self.class3 = m.group(2)
            self.number = m.group(1)
        else:
            raise Exception("ERROR: Attempt to add an EC number ({0}) in unrecognized format.  Expected N.N.N.N (where N can be 0-9 or a dash)".format(self.number))
        
