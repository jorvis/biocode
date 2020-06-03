Hey, let's keep a change log.

next:
----
- Added ORF class and mRNA.longest_orf() method

0.9.0
-----
- Added ncRNA feature class and changed inheritance for all related types (rRNA, tRNA, etc.)
- gff.parse_annotation_from_column_9() now recognizes both 'product_name' and 'product' in 9th column.

0.8.0
-----
- This actually should have been what the previous release was called. X.X.N releases
  should only be for bug fixes.
- Added dependency support for R in Docker build
- Added report_go_slim_counts.py script
- Added metaphlan-parsing R script (though needs major clean-up)

0.7.1
-----
- NCBI API key now required for all E-Utility calls
- Added --skip_existing option to batch genbank downloader" download_assemblies_from_genbank.py
- Changed genbank/download_assemblies_from_genbank.py to export BOTH contig and scaffold range files automatically
- Bugfix: genbank module was failing if the gene symbol contained a comma
- Corrected python module dependency from igraph -> python-igraph
- Updated documentation on taxadb generation in general/filter_uniref_by_taxonomy.py
- Added taxadb to formal dependencies
- Changed BED export library to use gene.id for 4th column if no locus is found
- Added ncbigff module and script to use it for conversions
- Added general/join_columnar_files.py
- Added fasta/fasta_simple_stats.py

0.7.0
-----
- Added script to convert GFF3 to BED
- Added script to convert BLAST/RAPSearch2 to BED
- Added script to convert HMMer3 HTAB to BED
- Added biocode.bed library
	

