Overview
========

Over the years most bioinformatics people amass a collection of small utility scripts 
which make their lives easier. Too often they are kept either in private repositories 
or as part of a public collection to which noone else can contribute. Biocode is a 
curated repository of general-use utility scripts my colleagues and I have found 
useful and want to share with others.  I have also developed some code libraries/modules
which have made my scripting work a lot easier.  Some have found these to be more
useful than the scripts themselves.

Look below if you want to learn more, contribute code yourself, or just get the 
scripts.

-- Joshua Orvis

Scope and directory structure
=============================

The scope here is intentionally very open. I want to include anything that developers 
find generally useful. There are no limitations on language choice. For now, the 
following directories make up the initial groupings but will be expanded as needed:

- blast - It if uses, massages, or just reformats BLAST output, it goes here.
- chado - Scripts that are tied into the chado schema (gmod.org) should be found here.
- fasta - Filtering, converting, size distribution plots, etc.
- fastq - Utilities for fasta's newer sister format.
- genbank - Anything related to the GenBank? Flat File Format.
- general - Utility scripts that may not fit in any other existing directory or don't 
  warrant creation of their own. We should be selective about what we put here and 
  create or use other directories whenever appropriate.
- gff - Extractions, conversions and manipulations of files in the Generic Feature 
  Format (sequenceontology.org/gff3.shtml)
- gtf - From Ensembl/WashU, the GTF format is the focus of scripts here.
- hmm - Merging, manipulating or reading HMM libraries.
- sam_bam - Analysis of and parsing SAM/BAM files.
- sandbox - Each committer gets their own personal directory here to add anything they 
  want while testing or waiting to be moved to the production directories.
- sysadmin - While not specifically bioinformatics, our work tends to be on Unix machines, 
  and utility scripts are often needed to support our work. From file system 
  manipulation to database backup scripts, put your generic sysadmin utilities here.
- taxonomy - Anything related to taxonomic analysis.

Getting the code
================

The best way to get the code for now is to just clone it:

   git clone https://github.com/jorvis/biocode.git

Problems / Suggestions?
=======================

If you encounter any issues with the existing code, or would like to request new
features or scripts please submit to the [Issue tracking system](https://github.com/jorvis/biocode/issues).

Contributing
============

If you'd like to contribute code to this collection have a look at the [Requirements And Convention Guide](https://github.com/jorvis/biocode/blob/master/RequirementsAndConventionGuide.md) 
and then submit a pull request once your code is ready.  We'll check your script 
and pull it into the production directories.  If you're not that confident yet 
we'll happily pull in your sandbox directory if you'd like to add your code to the 
project but aren't sure if it's ready to be in the production directories yet.


