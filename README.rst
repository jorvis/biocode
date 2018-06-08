Overview
========

This is a collection of bioinformatics scripts many have found useful
and code modules which make writing new ones a lot faster.

Over the years most bioinformatics people amass a collection of small
utility scripts which make their lives easier. Too often they are kept
either in private repositories or as part of a public collection to
which noone else can contribute. Biocode is a curated repository of
general-use utility scripts my colleagues and I have found useful and
want to share with others. I have also developed some code
libraries/modules which have made my scripting work a lot easier. Some
have found these to be more useful than the scripts themselves.

Look below if you want to learn more, contribute code yourself, or just
get the scripts.

-- Joshua Orvis

The scripts
===========

The scope here is intentionally very open. I want to include anything
that developers find generally useful. There are no limitations on
language choice, though the majority are Python. For now, the following
directories make up the initial groupings but will be expanded as
needed:

-  blast - It if uses, massages, or just reformats BLAST output, it goes
   here.
-  chado - Scripts that are tied into the chado schema (gmod.org) should
   be found here.
-  fasta - Filtering, converting, size distribution plots, etc.
-  fastq - Utilities for fasta's newer sister format.
-  genbank - Anything related to the GenBank? Flat File Format.
-  general - Utility scripts that may not fit in any other existing
   directory or don't warrant creation of their own. We should be
   selective about what we put here and create or use other directories
   whenever appropriate.
-  gff - Extractions, conversions and manipulations of files in the
   `Generic Feature Format <http://sequenceontology.org/gff3.shtml>`__
-  gtf - From Ensembl/WashU, the GTF format is the focus of scripts
   here.
-  hmm - Merging, manipulating or reading HMM libraries.
-  sam\_bam - Analysis of and parsing SAM/BAM files.
-  sandbox - Each committer gets their own personal directory here to
   add anything they want while testing or waiting to be moved to the
   production directories.
-  sysadmin - While not specifically bioinformatics, our work tends to
   be on Unix machines, and utility scripts are often needed to support
   our work. From file system manipulation to database backup scripts,
   put your generic sysadmin utilities here.
-  taxonomy - Anything related to taxonomic analysis.

The modules
===========

If you're a developer these modules can save a lot of time. Yes, there
is some duplicate functionality you'll find in modules like
`Biopython <http://biopython.org/wiki/Main_Page>`__, but these were
written to add features I always wanted and with a more
biologically-focused API.

Three of the primary Python modules:

`biocode.things <https://github.com/jorvis/biocode/blob/master/lib/biocode/things.py>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Classes here represent biological things (as defined by the `Sequence
Ontology <http://sequenceontology.org/>`__) in a way that makes more
sense biologically and hiding some of the CS abstraction. What does this
mean? This is a simple example, but compare these syntax approaches:

::

    # This way is typical of other libraries
    genes = assembly.get_subfeatures_by_type( 'type': 'genes' )
    mRNAs = assembly.get_subfeatures_by_type( 'type': 'mRNA' )

    # And instead, in biothings:
    genes = assembly.genes()
    for gene in genes:
        mRNAs = gene.mRNAs()

This more direct approach is held throughout these libraries. It also
adds some shortcuts for tasks that always annoyed me when working with
things that had coordinates. Consider if you wanted to determine if one
gene is before another one on a molecule:

::

    if gene1 < gene2:
        return True

In the background, biocode checks if the two gene objects are located on
the same molecule and, if so, compares their coordinates. There are many
other methods for coordinate comparison, such as:

-  thing1 <= thing2 : The thing1 overlaps thing2 on the 5' end
-  thing1.contained\_within( thing2 )
-  thing1.overlaps( thing2 )
-  thing1.overlap\_size\_with( thing2 )

This module also contains readable and detailed documention within the
`source
code <https://github.com/jorvis/biocode/blob/master/lib/biocode/things.py>`__.

`biocode.annotation <https://github.com/jorvis/biocode/blob/master/lib/biocode/annotation.py>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This set of classes allows formal definition of functional annotation
which can be attached to various biothings. These include gene product
names, gene symbols, EC numbers, GO terms, etc. Once annotated, the
biothings can be written out in common formats such as GFF3, GenBank,
NCBI tbl, etc.

`biocode.gff <https://github.com/jorvis/biocode/blob/master/lib/biocode/gff.py>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Much of biocode was written while working with genomic data and
annotation, and one of the more common formats for storing these is
`GFF3 <http://sequenceontology.org/resources/gff3.html>`__. Using this
module, you can parse a GFF3 file of annotations into a set of biothings
with a single line of code. For example:

::

    import biocode.gff

    (assemblies, features) = biocode.gff.get_gff3_features( input_file_path )

That's it. You can then iterate over the assemblies and their children,
or access the 'features' dict, which is keyed on each feature's ID.

Getting the code (pip3, latest release)
======================================

You can install biocode using pip3 (requires Python3) like this:

::

    pip3 install biocode

Getting the code (github, current trunk)
========================================

If you want the latest developer version:

::

    git clone https://github.com/jorvis/biocode.git

**Important**: Many of these scripts use the modules in the biocode/lib
directory, so you'll need to point Python to them. Full setup example:

::

    cd /opt
    git clone https://github.com/jorvis/biocode.git

    # You probably want to add this line to your $HOME/.bashrc file
    export PYTHONPATH=/opt/biocode/lib:$PYTHONPATH

Problems / Suggestions?
=======================

If you encounter any issues with the existing code, or would like to
request new features or scripts please submit to the `Issue tracking
system <https://github.com/jorvis/biocode/issues>`__.

Contributing
============

If you'd like to contribute code to this collection have a look at the
`Requirements And Convention
Guide <https://github.com/jorvis/biocode/blob/master/RequirementsAndConventionGuide.md>`__
and then submit a pull request once your code is ready. We'll check your
script and pull it into the production directories. If you're not that
confident yet we'll happily pull in your sandbox directory if you'd like
to add your code to the project but aren't sure if it's ready to be in
the production directories yet.
