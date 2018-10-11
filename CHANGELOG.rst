Hey, let's keep a change log.

0.7.1
-----
- NCBI API key now required for all E-Utility calls
- Added --skip_existing option to batch genbank downloader" download_assemblies_from_genbank.py
- Changed genbank/download_assemblies_from_genbank.py to export BOTH contig and scaffold range files automatically
- Bugfix: genbank module was failing if the gene symbol contained a comma
- Corrected python module dependency from igraph -> python-igraph
- Updated documenation on taxadb generation in general/filter_uniref_by_taxonomy.py
- Added taxadb to formal dependencies
- Changed BED export library to use gene.id for 4th column if no locus is found
- Added ncbigff module and script to use it for conversions

0.7.0
-----
- Added script to convert GFF3 to BED
- Added script to convert BLAST/RAPSearch2 to BED
- Added script to convert HMMer3 HTAB to BED
- Added biocode.bed library
	

