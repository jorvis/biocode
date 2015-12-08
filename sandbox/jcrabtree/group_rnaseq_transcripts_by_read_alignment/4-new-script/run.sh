#!/bin/tcsh

setenv WDIR /export/biocode-test
setenv BTEST2 $WDIR/hits_partial_1M.bam

# --------------------------
# TEST2 - 1M alignments
# --------------------------
# 10.3s, 10.2s, 9.9s (versus 2m16s for the original script, doing both parts of the analysis)
samtools view -h $BTEST2 | ./find_paired_rnaseq_transcripts_by_read_alignment.py -mmpc 3 > test2-pairs-r1-3-250bp.txt

# check that the pairings do not change from run to run
samtools view -h $BTEST2 | ./find_paired_rnaseq_transcripts_by_read_alignment.py -mmpc 3 > test2-pairs-r2-3-250bp.txt
samtools view -h $BTEST2 | ./find_paired_rnaseq_transcripts_by_read_alignment.py -mmpc 3 > test2-pairs-r3-3-250bp.txt
samtools view -h $BTEST2 | ./find_paired_rnaseq_transcripts_by_read_alignment.py -mmpc 3 > test2-pairs-r4-3-250bp.txt

#md5sum *.txt
#eff7ef966a5448f7d8aecb9236e6b12d  test2-pairs-r1-3-250bp.txt
#eff7ef966a5448f7d8aecb9236e6b12d  test2-pairs-r2-3-250bp.txt
#eff7ef966a5448f7d8aecb9236e6b12d  test2-pairs-r3-3-250bp.txt
#eff7ef966a5448f7d8aecb9236e6b12d  test2-pairs-r4-3-250bp.txt

# group pairs with a second script:


# ---------------------------------------------------
# FULL INPUT FILE - 1.17 billion alignments
# ---------------------------------------------------

setenv BAM /local/scratch-hs/aplysia/A1_CNS_unsheared_vs_merged/tophat_out/accepted_hits.bam

# compute pairs of linked transcripts
samtools view -h $BAM | ./find_paired_rnaseq_transcripts_by_read_alignment.py -mmpc 3 > aplysia-pairs-3-250bp.txt
gzip aplysia-pairs-3-250bp.txt

#zcat aplysia-pairs-3-250.txt.gz | md5sum
#fe5a1e54f724ef873ffffc5eb59590ea  -

#zcat aplysia-pairs-3-250.txt.gz | wc
# 251989  251989 8890403

# run connected components analysis to group transcripts
zcat aplysia-pairs-3-250bp.txt | ./group_paired_rnaseq_transcripts.py > aplysia-groups-3-250bp.txt
