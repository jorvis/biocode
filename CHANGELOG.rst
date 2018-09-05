Hey, let's keep a change log.

0.7.1
-----
- Bugfix: genbank module was failing if the gene symbol contained a comma
- Corrected python module dependency from igraph -> python-igraph
- Updated documenation on taxadb generation in general/filter_uniref_by_taxonomy.py
- Added taxadb to formal dependencies
- Changed BED export library to use gene.id for 4th column if no locus is found

0.7.0
-----
- Added script to convert GFF3 to BED
- Added script to convert BLAST/RAPSearch2 to BED
- Added script to convert HMMer3 HTAB to BED
- Added biocode.bed library
	

