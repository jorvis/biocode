import os
import re
import sys

## used for nt reverse complements
_nt_comp_table = bytes.maketrans(b'ACBDGHKMNSRUTWVYacbdghkmnsrutwvy', \
                                 b'TGVHCDMKNSYAAWBRtgvhcdmknsyaawbr')

# for translations: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
_translation_table = {
    1: { 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'AGA':'R', 'AGG':'R',
        'AAC':'N', 'AAT':'N',
        'GAC':'D', 'GAT':'D', 
        'TGC':'C', 'TGT':'C',
        'GAA':'E', 'GAG':'E',
        'CAA':'Q', 'CAG':'Q',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'CAC':'H', 'CAT':'H',
        'ATA':'I', 'ATC':'I', 'ATT':'I',
        'TTA':'L', 'TTG':'L', 'CTA':'L', 'CTC':'L', 'CTG':'L','CTT':'L',
        'AAA':'K', 'AAG':'K',
        'ATG':'M',
        'TTC':'F', 'TTT':'F',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'AGC':'S', 'AGT':'S',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'TGG':'W',
        'TAC':'Y', 'TAT':'Y',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'TAA':'*', 'TAG':'*', 'TGA':'*'
        }
}

def translate( seq, translation_table=None ):
    """
    Does a direct translation of the passed DNA/RNA sequence in phase 0.
    You can pass a numeric translation table, else 1 is assumed.

    An 'X' is used when a codon is unknown (contains an N)

    http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    """
    if translation_table is None:
        translation_table = 1

    # make sure we've defined this translation table
    if translation_table not in _nt_comp_table:
        raise Exception("ERROR: translation requested using table {0}, which isn't yet supported.".format(translation_table))

    trans_table = _translation_table[translation_table]
    
    # In case an RNA string was passed
    seq = seq.translate(seq.maketrans('Uutagc', 'TTTAGC'))

    polypeptide_seq = ''
    x = 0

    while True:
        try:
            polypeptide_seq += trans_table[seq[x:x+3]]
        except (IndexError):
            break
        except (KeyError):
            if len(seq[x:x+3]) == 3:
                #raise Exception("ERROR: Encountered unknown codon during translation: {0}".format(seq[x:x+3]))
                print("WARN: Encountered unknown codon during translation: {0}".format(seq[x:x+3]))
                polypeptide_seq += 'X'
            else:
                break
        x += 3
        
    
    return polypeptide_seq


def read_list_file( file ):
    """Parse an list file and return an array of the paths"""
    files = []

    if ( not os.path.isfile(file) ):
        raise Exception("Couldn't find file: " + file)

    ## only do non-blank lines
    with open(file) as f_in:
        lines = filter(None, (line.rstrip() for line in f_in))

        for line in lines:
            files.append(line)

    return files


def reverse_complement( seq ):
    """
    Biological reverse complementation.  Case in sequences are retained, and 
    IUPAC codes are supported.  Code modified from:

    http://shootout.alioth.debian.org/u32/program.php?test=revcomp&lang=python3&id=4
    """
    return seq.translate(_nt_comp_table)[::-1]
    

def humancoords_to_0interbase( start, stop ):
    """
    The typical human-readable coordinate system, such as found in GBK flat files,
    has a start and stop coordinate only.  They are 1-based, on-base coordinates
    and features on a reverse strand are indicated by having start > stop.  This
    transforms them into the GMOD standard 0-based inter-base coordinates.

    Returns a list of fmin, fmax and strand values.
    """
    fmin = start
    fmax = stop
    strand = 1
    
    if ( stop < start ):
        fmin = stop
        fmax = start
        strand = -1

    fmin -= 1

    return (fmin, fmax, strand)

def interbase0_to_humancoords( fmin, fmax, strand ):
    """
    Takes GMOD-standard 0-based inter-base coordinates and transforms them
    into the typical human-readable coordinates used in places like GBK flat
    files.  Features on the reverse strand have a larger start than stop.

    The strand values can be:

    + or 1 for forward features
    - or -1 for reverse features

    Returns a list: [start, stop] where start < stop if on forward strand.
    """

    if strand == '+' or strand == 1:
        return (fmin + 1, fmax)
    elif strand == '-' or strand == -1:
        return (fmax, fmin + 1)
    else:
        raise Exception("Invalid strand specified ({0}).  Expected +, -, 1 or -1".format(strand))


def fasta_dict_from_file( file ):
    """
    Reads a file of FASTA entries and returns a dict where each key is a sequence ID.
    The value is another dict with two keys 'h' for header and 's' for sequence.  The
    header is all the other text after the id in the original FASTA header.  The
    sequence has all whitespace removed.  Obviously this should only be used on files
    where memory to load them isn't an issue.
    """
    seqs = dict()
    current_seq = ''
    current_id = None
    current_header = None
    
    for line in open(file):
        line = line.rstrip()
        m = re.search('>(\S+)\s*(.*)', line)
        if m:
            ## new residue line matched, purge the existing one, if not the first
            if current_id is not None:
                ## warn if it has already been found
                if current_id in seqs:
                    sys.stderr.write("WARN: Duplicate ID ({0}) found.  Only last one kept.\n".format(current_id))

                ## remove all whitespace and save
                current_seq = ''.join(current_seq.split())
                seqs[current_id] = {'h':current_header, 's':current_seq}
                    
            current_seq = ''
            current_id = m.group(1)
            current_header = m.group(2)
        else:
            ## python 2.6+ makes string concatenation amortized O(n)
            ##  http://stackoverflow.com/a/4435752/1368079
            current_seq += str(line)

    ## don't forget the last one
    current_seq = ''.join(current_seq.split())
    seqs[current_id] = {'h':current_header, 's':current_seq}

    return seqs


def fasta_sizes_from_file(file):
    """
    Reads a file of FASTA entries and returns a dict where each key is a sequence ID and
    the value is just the size of each sequence.  This is as almost as computationally
    intensive as fasta_dict_from_file() but takes less memory and is appropriate when you
    only care about the residue lengths.
    """
    seqs = dict()
    current_id = None
    
    for line in open(file):
        line = line.rstrip()
        m = re.search('>(\S+)\s*(.*)', line)
        if m:
            ## new residue line matched, set the new seq ID
            current_id = m.group(1)
            seqs[current_id] = 0
        else:
            seqs[current_id] += len(line)

    return seqs


def add_assembly_fasta(mols, fasta_file):
    """
    Takes a collection of molecules (Assembly) objects and adds sequence residues based on
    their IDs being found in the passed FASTA file.
    """
    fasta_seqs = fasta_dict_from_file(fasta_file)

    for mol_id in mols:
        # check if the FASTA file provides sequence for this
        if mol_id in fasta_seqs:
            mol = mols[mol_id]
            mol.residues = fasta_seqs[mol_id]['s']
            mol.length   = len(mol.residues)
            
def wrapped_fasta(string, every=60):
    """
    Pass a string of residues (nucleotide or polypeptides) that has NO whitespace and
    this will return another string with new line characters inserted every N residues.

    N is specified with the 'every' parameter with a default of 60.  For example:

    new_fasta = wrapped_fasta(some_string, every=60)
    
    This runs orders of magnitude faster than using the textwrap module.
    """
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))
