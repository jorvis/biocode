#!/usr/bin/env python3

"""

This script is used to create a release of biocode for submission to PyPi. 

Version numbering: This project uses a sequence-based (micro) scheme for 
assigning version numbers: 

   http://the-hitchhikers-guide-to-packaging.readthedocs.io/en/latest/specification.html#sequence-based-scheme

These numbers are like X.Y.Z where X is a major release (rare, involving
major restructuring of the code and binary incompatibility with previous 
releases).  Y are the standard 'new features added' updates, and Z are 
bugfix/minor releases off any given Y version.

Run it from within a GitHub checkout of Biocode, and it will:

- Create the output directory if necessary (should be called 'biocode')
- Creates the biocode/biocode/__init__.py file (empty)
- Creates and populates the biocode/setup.py file
  - All files ending in .py under 'lib/' are assumed to be packages
- Copies all lib/biocode/*.py to biocode/biocode/
- Copies all .py scripts to biocode/bin/
- Copy all data/* to biocode/data/

"""

import argparse
import json
import os
import re
import shutil

def main():
    parser = argparse.ArgumentParser( description='Creates a directory structure of biocode for PyPi (or install)')

    ## output file to be written
    parser.add_argument('-o', '--output_directory', type=str, required=True, help='This is the directory where the output will be written' )
    parser.add_argument('-v', '--version', type=str, required=True, help='Version number of the release to be made' )
    args = parser.parse_args()

    # make sure it ends with lower-case biocode
    if not args.output_directory.endswith('biocode'):
        raise Exception("The last element of the path for --output_directory needs to be lower-case 'biocode'")

    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    # create the subdirectory of the same name required for PyPi submission
    lib_path = "{0}/biocode".format(args.output_directory)
    os.mkdir(lib_path)

    # create the biocode/biocode/__init__.py file
    init_fh = open("{0}/__init__.py".format(lib_path), 'wt')
    init_fh.close()

    # place the README
    shutil.copy('README.rst', "{0}/".format(args.output_directory))

    # get packages under lib/*.py and copy to release directory
    package_names = get_biocode_package_names('lib/biocode/')
    place_package_files(args.output_directory, package_names)

    script_paths = get_biocode_script_paths('.')
    bin_dir = "{0}/bin".format(args.output_directory)
    os.mkdir(bin_dir)
    place_script_files(bin_dir, script_paths)
    placed_script_paths = get_placed_script_paths(bin_dir)

    # make a JSON file to use for the script index
    data_dir = "{0}/data".format(lib_path)
    os.mkdir(data_dir)
    create_script_index_json('.', data_dir)

    # place the data files
    datafile_paths = get_biocode_datafile_paths('.')
    place_data_files(data_dir, datafile_paths)
    with open(os.path.join(data_dir, '__init__.py'), 'w'):
        pass

    # creates the biocode/setup.py file
    setup_fh = open("{0}/setup.py".format(args.output_directory), 'wt')
    populate_setup_file(setup_fh, args.version, placed_script_paths)
    setup_fh.close()

    # create the manifest
    manifest_fh = open("{0}/MANIFEST.in".format(args.output_directory), 'wt')
    manifest_fh.write("include README.rst\n")
    # perhaps do this instead of graft?  recursive-include biocode/data *.template
    manifest_fh.write("graft biocode/data\n")

    # include the license
    shutil.copy('LICENSE', "{0}/".format(args.output_directory))
    manifest_fh.write("include LICENSE\n")

    manifest_fh.close()

    print("Steps to tag (on GitHub) and release (on PyPi) version {0}:".format(args.version))
    print("\nNext do:\n\t$ cd {0} && python3 setup.py sdist && twine upload dist/*".format(args.output_directory))
    print("Then: \n\t$ cd - && git tag v{0} && git push --tags".format(args.version))

def create_script_index_json(script_base, data_dir):
    d = get_biocode_script_dict(script_base)
    json_data = dict()

    for category in sorted(d):
        json_data[category] = list()

        for script in sorted(d[category]):
            desc = None

            for line in open("{0}/{1}".format(category, script)):
                m = re.match('.+ArgumentParser.+description=\'(.+?)\'', line)
                if m:
                    desc = m.group(1)

            # TODO: Parse description lines from within each script
            json_data[category].append({'name': script, 'desc': desc})

    # Don't change the path of this JSON file without also changing it in the general/list_biocode utility
    with open("{0}/biocode_script_index.json".format(data_dir), 'w') as f:
        f.write(json.dumps(json_data))
    
def get_biocode_datafile_paths(base):
    datafiles = list()

    for directory in os.listdir(base):
        if directory not in ['data']: continue

        dir_path = "{0}/{1}".format(base, directory)
        if os.path.isdir(dir_path):
            for thing in os.listdir(dir_path):
                if thing == '__init__.py': continue
                thing_path = "{0}/{1}".format(dir_path, thing)

                if os.path.isfile(thing_path):
                    datafiles.append(thing_path)

    if len(datafiles) == 0:
        raise Exception("Error: failed to find any data files under {0}/data".format(base))

    return datafiles
    
def get_biocode_package_names(base):
    packages = list()
    
    for thing in os.listdir(base):
        if thing.endswith('.py'):
            thing = thing.rstrip('.py')
            packages.append(thing)

    if len(packages) == 0:
        raise Exception("Error: failed to find any .py files under lib/")

    return packages

def get_biocode_script_paths(base):
    scripts = list()

    for directory in os.listdir(base):
        if directory in ['build', '.git', 'sandbox', 'tests']: continue

        dir_path = "{0}/{1}".format(base, directory)
        if os.path.isdir(dir_path):
            for thing in os.listdir(dir_path):
                if thing == '__init__.py': continue

                if thing.endswith('.py'):
                    scripts.append("{0}/{1}".format(dir_path, thing))

    if len(scripts) == 0:
        raise Exception("Error: failed to find any .py files under */")

    return scripts

def get_biocode_script_dict(base):
    """
    Returns a dict where each key is a script category 'fasta', 'blast', etc.
    and the values are the names of the scripts within that directory
    """
    d = dict()

    for directory in os.listdir(base):
        if directory in ['build', '.git', 'sandbox', 'tests']: continue

        dir_path = "{0}/{1}".format(base, directory)

        if os.path.isdir(dir_path):
            d[directory] = list()

            for thing in os.listdir(dir_path):
                if thing == '__init__.py': continue

                if thing.endswith('.py'):
                    d[directory].append(thing)

    return d
    

def get_placed_datafile_paths(base):
    paths = list()

    for thing in os.listdir(base):
        if os.path.isfile("{0}/{1}".format(base, thing)):
            paths.append("data/{0}".format(thing))
            
    if len(paths) == 0:
        raise Exception("Failed to find any data files in base: {0}".format(base))

    return paths

def get_placed_script_paths(base):
    paths = list()

    for thing in os.listdir(base):
        if thing.endswith('.py'):
            paths.append("bin/{0}".format(thing))
            
    if len(paths) == 0:
        raise Exception("Failed to find any py scripts in base: {0}".format(base))

    return paths

def place_data_files(dest_dir, paths):
    for path in paths:
        shutil.copy(path, dest_dir)

def place_package_files(base, names):
    for name in names:
        shutil.copyfile("lib/biocode/{0}.py".format(name), "{0}/biocode/{1}.py".format(base, name))

def place_script_files(dest_dir, script_paths):
    for path in script_paths:
        shutil.copy(path, dest_dir)
        
def populate_setup_file(fh, version, script_paths):
    settings = {
        'author': 'Joshua Orvis',
        'email': 'jorvis@gmail.com',
        'version': version,
        'packages': ['biocode'],
        'dependencies': ['python-igraph', 'jinja2', 'matplotlib', 'numpy', 'taxadb'],
        'scripts': script_paths
    }

    fh.write("""
from setuptools import setup

read = lambda f: open(f, 'r').read()

setup(name='biocode',
      author='{author}',
      author_email='{email}',
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      description='Bioinformatics code libraries and scripts',
      include_package_data=True,
      install_requires={dependencies},
      keywords='bioinformatics scripts modules gff3 fasta fastq bam sam',
      license='MIT',
      long_description=read('README.rst'),
      packages={packages},
      scripts={scripts},
      url='http://github.com/jorvis/biocode',
      version='{version}',
      zip_safe=False)
    """.format(**settings))

if __name__ == '__main__':
    main()







