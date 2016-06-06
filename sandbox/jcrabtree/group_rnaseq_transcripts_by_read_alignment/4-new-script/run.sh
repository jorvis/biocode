#!/bin/tcsh

setenv WDIR /export/biocode-test
setenv BTEST2 $WDIR/hits_partial_1M.bam

# --------------------------
# TEST2 - 1M alignments
# --------------------------
# 10.3s, 10.2s, 9.9s (versus 2m16s for the original script, doing both parts of the analysis)
#samtools view -h $BTEST2 | ./find_paired_rnaseq_transcripts_by_read_alignment.py -mmpc 3 > test2-pairs-r1-3-250bp.txt

#wc test2-pairs-r1-3-250bp.txt 
# 265  265 9507 test2-pairs-r1-3-250bp.txt

# check that the pairings do not change from run to run
#samtools view -h $BTEST2 | ./find_paired_rnaseq_transcripts_by_read_alignment.py -mmpc 3 > test2-pairs-r2-3-250bp.txt
#samtools view -h $BTEST2 | ./find_paired_rnaseq_transcripts_by_read_alignment.py -mmpc 3 > test2-pairs-r3-3-250bp.txt
#samtools view -h $BTEST2 | ./find_paired_rnaseq_transcripts_by_read_alignment.py -mmpc 3 > test2-pairs-r4-3-250bp.txt

#md5sum *.txt
#eff7ef966a5448f7d8aecb9236e6b12d  test2-pairs-r1-3-250bp.txt
#eff7ef966a5448f7d8aecb9236e6b12d  test2-pairs-r2-3-250bp.txt
#eff7ef966a5448f7d8aecb9236e6b12d  test2-pairs-r3-3-250bp.txt
#eff7ef966a5448f7d8aecb9236e6b12d  test2-pairs-r4-3-250bp.txt

# group pairs with a second script:
# 0.2s, - 187 nodes, 265 edges
#./group_paired_rnaseq_transcripts.py -i test2-pairs-r1-3-250bp.txt > test2-groups-r1-3-250bp-g1.txt

# check that results don't change

#./group_paired_rnaseq_transcripts.py -i test2-pairs-r1-3-250bp.txt > test2-groups-r1-3-250bp-g2.txt
#./group_paired_rnaseq_transcripts.py -i test2-pairs-r1-3-250bp.txt > test2-groups-r1-3-250bp-g3.txt
#./group_paired_rnaseq_transcripts.py -i test2-pairs-r1-3-250bp.txt > test2-groups-r1-3-250bp-g4.txt

#md5sum test2-groups-r1-3-250bp-g?.txt

# ---------------------------------------------------
# FULL INPUT FILE - 1.17 billion alignments
# ---------------------------------------------------

setenv BAM /local/scratch-hs/aplysia/A1_CNS_unsheared_vs_merged/tophat_out/accepted_hits.bam

# compute pairs of linked transcripts
#samtools view -h $BAM | ./find_paired_rnaseq_transcripts_by_read_alignment.py -mmpc 3 > aplysia-pairs-3-250bp.txt
#gzip aplysia-pairs-3-250bp.txt

#zcat aplysia-pairs-3-250bp.txt.gz | md5sum
#fe5a1e54f724ef873ffffc5eb59590ea  -

#zcat aplysia-pairs-3-250bp.txt.gz | wc -l
# 251989

# run connected components analysis to group transcripts
# 2.2s, 2.2s, 2.2s - 139228 nodes, 251989 edges, result = 29491 groups
gzcat aplysia-pairs-3-250bp.txt.gz | ./group_paired_rnaseq_transcripts.py > aplysia-groups-3-250bp-1.txt

# check results appear to be deterministic

#gzcat aplysia-pairs-3-250bp.txt | ./group_paired_rnaseq_transcripts.py > aplysia-groups-3-250bp-2.txt
#gzcat aplysia-pairs-3-250bp.txt | ./group_paired_rnaseq_transcripts.py > aplysia-groups-3-250bp-3.txt
#gzcat aplysia-pairs-3-250bp.txt | ./group_paired_rnaseq_transcripts.py > aplysia-groups-3-250bp-4.txt

# md5sum aplysia-groups*.txt

# alternate method to check node count
#gzcat aplysia-pairs-3-250bp.txt.gz | perl -e ' while (<>) { chomp; map {$h->{$_} = 1} split(":"); } print "nodes=" . scalar(keys %$h). "\n";'

