import os
import re

## used for nt reverse complements
_nt_comp_table = bytes.maketrans(b'ACBDGHKMNSRUTWVYacbdghkmnsrutwvy', \
                                 b'TGVHCDMKNSYAAWBRtgvhcdmknsyaawbr')

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

    Returns a list of fmin, fmax and strand values
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
                    sys.stderr.write("WARN: Duplicate ID ({0}) found.  Only last one kept.".format(current_id))

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
    
