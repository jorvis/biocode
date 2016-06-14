#!/bin/tcsh

setenv WDIR /export/biocode-test
setenv TEST1 $WDIR/hits_partial_100K.sam
setenv TEST2 $WDIR/hits_partial_1M.sam

setenv BTEST1 $WDIR/hits_partial_100K.bam
setenv BTEST2 $WDIR/hits_partial_1M.bam

# ------------------------------------------------------
# Check time needed by various parts of the program
# ------------------------------------------------------


# --------------------------
# TEST1 - 100K alignments
# --------------------------

# Read BAM file, skip header lines, split other lines on tab. Doesn't store anything.

# 1.2s, 1.2s
#./read_and_split_only.py -i $TEST1

# Same as above, but build nested dict keyed on transcript1/transcript2.
# Alignments themselves are not stored.

# 1.4s, 1.5s, 1.4s, 1.5s - 18 L1 keys
#./read_split_and_store.py -i $TEST1

# Same as above, but use a single flat dict with a composite key (the string tid1 + ":" + tid2)

# 1.4s, 1s, 1.4s, 1.4s - 101 L1 keys
#./read_split_and_store_compkey.py -i $TEST1

# Same as above but also do the CIGAR and read name parsing:

# 2s, 2s, 2.5s, 2.5s ~69% reads ignored
#./read_split_cigar_and_store.py -i $TEST1

# --------------------------
# TEST2 - 1M alignments
# --------------------------

# Same experiments as above. Remember that ~8-10 seconds is our baseline time to do a 
# samtools view | wc -l on the 1M alignment BAM file and a straight wc -l should only
# take about half a second.

# 5.5s, 5.4s, 6.6s
#./read_and_split_only.py -i $TEST2

# 7.4s, 6.3s, 6.3s - 138 L1 keys
#./read_split_and_store.py -i $TEST2

# 5.8s, 5.9s, 6.7s, 6.9s - 696 L1 keys
#./read_split_and_store_compkey.py -i $TEST2

# ~72% reads ignored
# 13.1s, 13.6s, 14.3s, 14.7s
#./read_split_cigar_and_store.py -i $TEST2

# One more step. Parse the CIGAR lines but this time store the alignments in the list().
# This appears to be the time (and memory)-consuming step.

# 2m13s
#./read_split_cigar_and_store_in_list.py -i $TEST2

