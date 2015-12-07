#!/bin/tcsh

setenv WDIR /export/biocode-test
setenv TEST1 $WDIR/hits_partial_100K.sam
setenv TEST2 $WDIR/hits_partial_1M.sam

setenv BTEST1 $WDIR/hits_partial_100K.bam
setenv BTEST2 $WDIR/hits_partial_1M.bam

# ------------------------------------------------------
# Tests with original version
# ------------------------------------------------------

# --------------------------
# TEST1 - 100K alignments
# --------------------------

# running on jcrabtreevm-lx /export - 7.55s, 8.18s, 9.9s, 8.9s  ~ 20M memory
# 12 groups
./group_rnaseq_transcripts_by_read_alignment.py -i $TEST1 -mmpc 3 -o test-1-r1.groups |& cat >test-1-r1.err
# *10* groups
./group_rnaseq_transcripts_by_read_alignment.py -i $TEST1 -mmpc 3 -o test-1-r2.groups |& cat >test-1-r2.err
# 12 groups
./group_rnaseq_transcripts_by_read_alignment.py -i $TEST1 -mmpc 3 -o test-1-r3.groups |& cat >test-1-r3.err

# samtools view the BAM and pipe into the script (might result in slightly less iowait)
# 6.64s, 7.65s
samtools view $BTEST1 | ./group_rnaseq_transcripts_by_read_alignment.py -mmpc 3 -o test-1b.groups |& cat >test-1b.err

# --------------------------
# TEST2 - 1M alignments
# --------------------------

# 2m18s, ~200M memory, 2m30s, 2m30s
# 78 groups
./group_rnaseq_transcripts_by_read_alignment.py -i $TEST2 -mmpc 3 -o test-2.groups |& cat >test-2.err

# samtools view the BAM (also skips the header, which the script ignores anyway)
# 2m16s, 
samtools view $BTEST2 | ./group_rnaseq_transcripts_by_read_alignment.py -mmpc 3 -o test-2b.groups |& cat >test-2b.err

# oops - the grouping algorithm isn't stable:

#[jcrabtree@jcrabtreevm-lx 1-original-script]$ wc *.groups
#   10    32   570 test-1b.groups
#   12    38   678 test-1-r1.groups
#   10    33   589 test-1-r2.groups
#   12    38   678 test-1-r3.groups
#   72   275  4923 test-2b.groups
#   75   319  5721 test-2.groups
#  191   735 13159 total
