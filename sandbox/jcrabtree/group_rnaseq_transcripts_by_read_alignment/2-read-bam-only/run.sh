#!/bin/tcsh

setenv WDIR /export/biocode-test
setenv TEST1 $WDIR/hits_partial_100K.sam
setenv TEST2 $WDIR/hits_partial_1M.sam

setenv BTEST1 $WDIR/hits_partial_100K.bam
setenv BTEST2 $WDIR/hits_partial_1M.bam

# ------------------------------------------------------
# Check the time needed to read the SAM and/or BAM
# ------------------------------------------------------

# each relevant line was uncommented and then "time ./run.sh" was run to obtain
# a few informal measurements (results noted below)

# --------------------------
# TEST1 - 100K alignments
# --------------------------

# 2.1s,1.9s,2.5s
#samtools view -S $TEST1 | wc -l

# 1.35s,1.37s,1.34s - speedup may be partially due to quicker skipping of header?
#samtools view $BTEST1 | wc -l

# 0.14s,0.13s,0.13s
#wc -l $TEST1

# NOTE: wc is MUCH slower if not using "-l". See http://unix.stackexchange.com/questions/96563/why-is-wc-so-slow
# 4.6s, 3.7s, 4.1s
#wc $TEST1

# --------------------------
# TEST2 - 1M alignments
# --------------------------

# 9s, 10.9s, 10.4s
#samtools view -S $TEST2 | wc -l

# 8.5s, 9.0s, 8.1s (uses a bit more CPU to achieve the speedup)
#samtools view $BTEST2 | wc -l

# .6s, .4s, .3s
wc -l $TEST2

# So a straight wc should take less than a second for either file, and
# a samtools view about 1-2 seconds on the small file and 8-10 seconds
# on the larger one.


