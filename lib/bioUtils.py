
import os

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
