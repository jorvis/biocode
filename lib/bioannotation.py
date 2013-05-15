'''
Warning: this module requires Python 3.2 or higher

It's also in an initial development phase, so either don't use it or be tolerant of
API changes.
'''

class FunctionalAnnotation:
    """
    While recognizing that an enormous variety of attributes could go here in
    describing the functional annotation of a BioThing, I'm starting with those
    we most often encounter and need to be available in common output formats.

    Also, there's a place for having attributes like this abstracted, stored in
    ontologies, etc.  We've done all that before.  For now I'm going to try
    and hopefully enjoy the utility of having the most direct properties always
    directly available.
    """
    def __init__( self, product_name=None, gene_symbol=None, go_annotations=None, ec_numbers=None ):
        self.product_name     = product_name
        self.gene_symbol      = gene_symbol
        self.go_annotations   = go_annotations
        self.ec_numbers       = ec_numbers

        if self.go_annotations is None:
            self.go_annotations = list()

        if self.ec_numbers is None:
            self.ec_numbers = list()

    def add_go_annotation( self, go ):
        self.append( go )
        
    def add_ec_number( self, ec_num ):
        self.append( ec_num )
        
class GOAnnotation:
    """
    A functional annotation can have an infinite number of associated GO Annotations

    Details here:
    http://www.geneontology.org/GO.evidence.shtml

    Yes, the 'with_from' attribute name is awkward, but 'with/from' isn't legal and
    both 'with' and 'from' are python reserved words.
    """
    def __init__( self, go_id=None, ev_code=None, with_from=None ):
        self.go_id        = go_id
        self.ev_code      = ev_code
        self.with_from    = with_from


