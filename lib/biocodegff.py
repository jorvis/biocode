import re
from urllib.parse import unquote

def column_9_value(value, key):
    '''
    Pass a case-sensitive key and this function will return the value,
    if present, else None.
    
    This was borrowed from the URL below and then modified for Python 3
    ftp://ftp.informatics.jax.org/%2Fpub/mgigff/gff3.py
    '''
    SEMI	= ';'
    COMMA	= ','
    EQ	        = '='
    WSP_RE = re.compile(r'^\s*$')
    
    if value == ".":
        return {}
    c9 = {}
    for t in value.split(SEMI):
        if WSP_RE.match(t):
            continue
        tt = t.split(EQ)
        if len(tt) != 2:
            raise ParseError("Bad column 9 format near '%s'."%t)
        n = unquote(tt[0].strip())
        [*v] = map(unquote, tt[1].strip().split(COMMA))
        if len(v) == 1:
            v = v[0]
        c9[n] = v

    if key in c9:
        # this is a list if there were comma-separated values
        return c9[key]
    else:
        return None


def parse_gff3_by_relationship( gff3_file ):
    '''
    Parses a GFF3 file, caring only about the ID/Parent relationships and returning a dict where
    they keys are the feature IDs and the values are an dict of dicts.   The depth of this depends
    on the depth of parent/child relationships, but every node looks like this:

    $feat_id = { children: [], fmin: 30, cols: [] }

    Every feature must have either and ID or Parent defined, else an error will be thrown.
    
    The structure returned looks like this:
    {
    $molecule_id = { 
                         fmin: 30,
                         cols: [...],
                         children: {
                             same as this structure
                         }

                   }, ...
    

    }

    WARNING:
    This will currently fail on any GFF that:
       - allows the same IDs on multiple lines, such as those for  discontiguous features (allowed in spec)
       - features with shared parents


    '''
    ## might try turning this into a defaultdict at a later point: http://ohuiginn.net/mt/2010/07/nested_dictionaries_in_python.html
    fgraph = dict()
    current_molecule_id = None
    ## key = id, value = parent_id
    parentage = dict()
    ## where features go during processing that have a parent which hasn't been seen yet
    purgatory = list()

    for line in open(gff3_file):
        cols = line.split("\t")

        if len(cols) == 9:
            cols[8] = cols[8].rstrip()
        else:
            continue

        mol_id = cols[0]
        id = column_9_value(cols[8], 'ID')
        # shared parents will cause failure here
        parent = column_9_value(cols[8], 'Parent')

        if id:
            parentage[id] = parent
    
        if mol_id != current_molecule_id:
            if current_molecule_id is not None:
                _reunite_children( fgraph, current_molecule_id, purgatory )

            ## reset for the next molecule
            purgatory = list()
            current_molecule_id = mol_id

            if current_molecule_id not in fgraph:
                fgraph[current_molecule_id] = dict()

        molecule = fgraph[current_molecule_id]

        if parent:
            uparent = _get_ultimate_parent( parentage, parent )
            molecule[uparent]['children'].append( {'id': id, 'cols': cols} )
        else:
            if id is None:
                raise Exception("ERROR: Encountered a line without a Parent or ID assigned: {0}".format(line))

            ## this shouldn't exist already
            if id in molecule:
                raise Exception("ERROR: found duplicate id ({0}) in input file".format(id) )

            molecule[id] = {'fmin': cols[3], 'children': list(), 'cols': cols}

    ## don't forget to handle the last one
    _reunite_children( fgraph, current_molecule_id, purgatory )

    return fgraph

    
        
def _reunite_children( fg, mol_id, kids ):
    pass



def _get_ultimate_parent( p, id ):
    if p is None:
        return id

    oldest = id

    while p[oldest] is not None:
        oldest = p[oldest]

    return oldest











