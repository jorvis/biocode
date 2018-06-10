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
    directly, and simply, available.
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
        TODO: Modify this to allow passing ECAnnotation object or string.
        Right now it expects an ECAnnotation object
        """
        self.ec_numbers.append(ec_num)

    def add_go_annotation(self, go):
        """
        TODO: Modify this to allow passing GOAnnotation object or string.
        Right now it expects an GOAnnotation object
        """
        self.go_annotations.append(go)

    def process_gene_symbol(self):
        """
        This method applies a series of rules based on our experience in annotation
        to gene symbols, attempting to correct things before they get submitted to
        GenBank.  This includes:

        - Removing anything after the first whitespace
        - DOES NOT change case.  Genbank demands that prok submissions all use lower
          case but allows organism-specific conventions to be followed for others.
          - https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/

        This method returns the new gene symbol rather than overwriting the attribute.  For
        that, use set_processed_gene_symbol instead.
        """
        new_gs = self.gene_symbol

        if not new_gs:
            return new_gs

        # take off everything after and including the first space, if any
        new_gs = new_gs.split(' ', 1)[0]

        return new_gs
        

    def process_product_name(self):
        """
        This method applies a series of rules based on years of annotation experience
        to gene product names, attempting to correct much of the mess of things which
        which get submitted to Genbank.  This includes:

        - Removing trailing periods
        - Nonredundifies: 'protein protein', 'family family', 'family protein family protein', 
          'protein family protein' and similar
        - Changes 'superfamily' to 'family'
        - Any starting with 'ORF' or 'orf' get changed to 'conserved hypothetical protein'
        - Any products which contain 'homolog' get changed to CHP*
        - Any products with 'similar to' get changed to CHP
        - Any with 'DUF' or 'UPF' get changed to CHP
        - Any with 'ncharacteri' (captures US/british spellings of 'uncharacterized'), gets changed to CHP
        - Select words are changed to their american spellings
        - Changes any of (predicted|possible|potential|probable) to 'putative'
        - Change 'TTG start' to CHP
        - Change any starting with 'residues' to CHP
        - Removes 'C-terminus' and 'N-terminus' from end of product name
        - Strips embedded Dbxrefs from end of name:
          - Example: (2E,6E)-farnesyl diphosphate synthase {ECO:0000313|EMBL:OOP19401.1}
        - Then a long list of manual, specific name changes

        * CHP = conserved hypothetical protein

        It returns the new product name rather than overwriting the attribute.  For that,
        use set_processed_product_name().

        Note:  NCBI submission rules:  https://www.ncbi.nlm.nih.gov/genbank/asndisc.examples/
               Especially the section: SUSPECT_PRODUCT_NAMES
               These are NOT currently all implemented here
        """
        new_product = self.product_name
        default_product = 'conserved hypothetical protein'

        # remove/replace troublesome strings
        new_product = new_product.rstrip('.')
        new_product = new_product.replace('protein protein', 'protein')
        new_product = new_product.replace('superfamily', 'family')
        new_product = new_product.replace('family family', 'family')
        new_product = new_product.replace('family protein family protein', 'family protein')
        new_product = new_product.replace('family protein domain protein', 'family protein')
        new_product = new_product.replace('domain domain protein', 'domain protein')
        new_product = new_product.replace('protein family protein', 'family protein')
        new_product = new_product.replace('superfamily protein family protein', 'family protein')
        new_product = new_product.replace(' Protein-like family protein', '-like family protein')
        new_product = new_product.replace(' protein-like family protein', '-like family protein')

        # takes patterns like this off the end:  {ECO:0000313|EMBL:OOP19401.1}
        m = re.match('(.+) \{.+\:.+\}', new_product)
        if m:
            new_product = m.group(1)
        
        if new_product.lower().startswith('orf'):
            return default_product

        if 'ncharacteri' in new_product.lower():
            return default_product
            
        # process some homolog-specific names
        if 'homolog' in new_product.lower():
            if 'shiA homolog' in new_product:
                new_product = 'shiA protein'
            elif 'virulence factor mviM homolog' in new_product:
                new_product = 'virulence factor mviM'
            elif 'protein phnA homolog' in new_product:
                new_product = 'phnA protein'
            elif 'protein seqA homolog' in new_product:
                new_product = 'seqA protein'
            else:
                return default_product

        if 'similar to' in new_product.lower():
            return default_product

        if 'DUF' in new_product or 'UPF' in new_product:
            return default_product

        # If it is any form of conserved hypothetical, return the proper version of that.
        if 'onserved hypothe' in new_product.lower():
            return default_product

        if 'unnamed' in new_product.lower():
            return default_product

        # Is the name *only* protein (with whitespace)
        if new_product.lower().lstrip().rstrip() == 'protein':
            return default_product

        # Some sources give products which start with the word 'residues'
        if new_product.lower().startswith('residues'):
            return default_product

        # Some proteins are simply named 'TTG start'
        if new_product.lower().startswith('ttg start'):
            return default_product

        # Correct a class of short bogus names we've often encountered
        m = re.match('^\w{1,2}\d{1,3}$', new_product)
        if m:
            return default_product

        m = re.match('gene \d+ protein', new_product)
        if m:
            return default_product

        m = re.match('(.*)\d*\s*[CN]\-terminus$', new_product)
        if m:
            new_product = m.groups(1)

        # removes trailing symbols
        new_product = new_product.rstrip('.,-_:/')

        # Names can't end in family.  Example replacements:
        #  Raf kinase inhibitor-like protein, YbhB/YbcL family -> YbhB/YbcL family Raf kinase inhibitor-like protein
        #  phage major capsid protein, HK97 family  ->  HK97 family phage major capsid protein
        #  phage portal protein, lambda family  ->  lambda family phage portal protein
        if new_product.endswith(' family'):
            m = re.match('(.+), (.+ family)', new_product)
            if m:
                new_product = "{0} {1}".format(m.group(2), m.group(1))

        # If family still remains in the name twice, take out the first one
        #  Peptidase family S49 family protein  ->  Peptidase S49 family protein
        if new_product.count('family') > 1:
            m = re.match('(.+?) family (.+)', new_product)
            if m:
                new_product = "{0} {1}".format(m.group(1), m.group(2))

        # Americanize some words.  I'm up for arguments against these, but adding them
        #  because our previous software had them.
        new_product.replace('utilisation', 'utilization')
        new_product.replace('utilising', 'utilizing')
        new_product.replace('dimerisation', 'dimerization')
        new_product.replace('disulphide', 'disulfide')
        new_product.replace('sulphur', 'sulfur')
        new_product.replace('mobilisation', 'mobilization')

        # standardize several different forms of 'putative' to a single one
        new_product = re.sub("predicted", "putative", new_product, flags=re.I)
        new_product = re.sub("possible", "putative", new_product, flags=re.I)
        new_product = re.sub("potential", "putative", new_product, flags=re.I)
        new_product = re.sub("probable", "putative", new_product, flags=re.I)

        # Fix incorrect spellings which are getting transitively annotated
        new_product = re.sub("putaive", "putative", new_product, flags=re.I)

        # Replacements requiring word boundaries
        patt = re.compile(r'\b(?:%s)\b' % 'asparate')
        new_product = re.sub(patt, "aspartate", new_product)

        # Now a long series of manual name changes we've gathered over the years
        #  the key is the source, value is what it will be changed to.
        replacement_products = {
            'alr5027 protein': 'heme-binding protein HutZ',
            'arginine-tRNA-transferase, C terminus family protein': 'putative arginine-tRNA-transferase',
            'bacterial regulatory helix-turn-helix proteins, AraC family protein': 'transcriptional regulator, AraC family',
            'bacterial regulatory proteins, gntR family protein': 'transcriptional regulator, GntR family',
            'bacterial regulatory proteins, lacI family protein': 'transcriptional regulator, LacI family',
            'bacterial regulatory proteins, luxR family protein': 'transcriptional regulator, LuxR family',
            'bacterial regulatory proteins, tetR family protein': 'transcriptional regulator, TetR family',
            'bordetella uptake gene (bug) product family protein': 'bug family protein',
            'conserved protein with nucleoside triphosphate hydrolase domain': 'putative ATP-dependent endonuclease',
            'cyclic di-GMP binding protein VCA0042': 'cyclic di-GMP binding protein',
            'cytochrome b(C-terminal)/b6/petD family protein': 'cytochrome b family protein',
            'domain related to MnhB subunit of Na+/H+ antiporter family protein': 'Na+/H+ antiporter family protein',
            'FAD linked oxidases, C-terminal domain protein': 'FAD linked oxidase domain protein',
            'FGGY family of carbohydrate kinases, C-terminal domain protein': 'carbohydrate kinase, FGGY family',
            'gene 25-like lysozyme family protein': 'lysozyme family protein',
            'glutamate synthases, NADH/NADPH, small subunit domain protein': 'glutamate synthase, NADH/NADPH, small subunit',
            'glycogen/starch synthases, ADP-glucose type family protein': 'glycogen/starch synthase',
            'GSPII_E N-terminal domain protein': 'bacteriophage N4 adsorption protein B',
            'hydro-lases, Fe-S type, tartrate/fumarate subfamily, beta region domain protein': 'fumarate hydratase family protein',
            'invasion gene expression up-regulator, SirB family protein': 'invasion gene expression up-regulator',
            'K+ potassium transporter family protein': 'potassium uptake protein',
            'menC_gamma/gm+: o-succinylbenzoic acid (OSB) synthetase': 'o-succinylbenzoic acid (OSB) synthetase',
            'phage/plasmid replication , gene II/X family protein': 'phage/plasmid replication protein, gene II/X family',
            'phospholipase d active site motif family protein': 'phospholipase D family protein',
            'PIII': default_product,
            'putative 2-hydroxyacid dehydrogenase HI_1556': 'putative 2-hydroxyacid dehydrogenase',
            'SULFATE TRANSPORTER SULFATE TRANSPORTER FAMILY PROTEIN': 'sulfate permease family protein',
            'thiamin/thiamin pyrophosphate ABC transporter, thiamin/thiamin pyrophospate-binding protein': 'thiamin/thiamine pyrophosphate ABC transporter, thiamin/thiamine pyrophospate-binding protein',
            'traG-like , N-terminal region family protein': 'putative traG protein',
            'transcriptional activator of defense systems': 'multiple antibiotic resistance protein MarA',
            'transcriptional regulatory protein, C terminal family protein': 'putative transcriptional regulator',
            'transposase and inactivated derivative': 'putative transposase',
            'tripartite ATP-independent periplasmic transporters, DctQ component family protein': 'tripartite ATP-independent periplasmic transporter, DctQ family',
            'type IV secretory pathway VirD2 components': 'type IV secretory pathway protein',
            'zn-dependent hydrolase of the beta-lactamase fold': default_product,
        }

        for old in replacement_products:
            if new_product == old:
                new_product = replacement_products[old]
                break

        return new_product.rstrip().lstrip()
        

    def set_processed_gene_symbol(self):
        """
        See FunctionalAnnotation.process_gene_symbol() for full list of actions
        """
        self.gene_symbol = self.process_gene_symbol()

    def set_processed_product_name(self):
        """
        See FunctionalAnnotation.process_product_name() for full list of actions
        """
        self.product_name = self.process_product_name()

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
    def __init__(self, number=None):
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
        
    def __repr__(self):
        return "EC:{0}".format(self.number)
