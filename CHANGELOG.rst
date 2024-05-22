Hey, let's keep a change log.

next
----
- Added script: fasta/strip_fasta_headers_after_regex.py

0.11.0
----
- Documentation updates
- New script: @RPINerd added fasta/fasta_size_report.py  
- New script: fastq/validate_mate_pairing.py
- New script: general/report_matrix_shape.py
- New script: gff/remove_fasta_from_gff3.py
- New script: gtf/convert_stringtie_gtf_to_gff3.py
- fasta/fasta_simple_stats.py: Critical fix - Script was missing stats on last/only sequence.
- general/filter_columnar_file.pl renamed to general/filter_columnar_file_by_column_values.pl
- lib/biocode/utils.py: Added function fasta_file_from_dict()
- lib/biocode/gff.py: gff.print_gff3_from_assemblies() no longer automatically splits genes with multiple mRNAs
- lib/biocode/gff.py: Fixed error causing failure if print_gff3_from_assemblies() was called
- fastq/randomly_subsample_fastq.py: Corrected shebang line which was causing failures
- fastq/fastq_simple_stats.py: Added (optional) progress interval reporting
- fastq/fastq_simple_stats.py: Added commas to integer outputs for better human readability
- blast/calculate_query_coverage_by_blast.py: Added error handling so it reports/dies gracefully if an input seq isn't found
- fasta/filter_fasta_by_ids.pl: minor: added additional info if in debugging mode
- general/join_columnar_files.py: Corrected issue where header line was printing to STDOUT rather than output file
- gff/convert_tRNAScanSE_to_gff3.pl: @pgonzale60 committed improvments to satisfy NCBI's table2asn_GFF requirements
- fastq/fastq_simple_stats.py: Fixed issues where script would fail if file had no entries
- Critical bugfix and general improvements to sandbox/export_pasa_unigenes.py

0.10.0
------
- Added things.Region class for generic genomic regions
- Updated dependency information for latest Ubuntu
- Added ORF class and mRNA.longest_orf() method
- Added script to filter FASTA files from abundance estimation (Salmon) results
- Added new script for classifying gene coverage in new assemblies
- Exports to GFF3 updated: tRNAs include anticodon, rRNAs includ product name
- mRNA.longest_orf() now has require_start option
- fasta_simple_stats.py now exports N50 and N90 stats
- fasta_simple_stats.py no longer errors on blank lines


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
	

