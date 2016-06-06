#!/bin/tcsh

# large (48G) input BAM file:
setenv BAM /local/scratch-hs/aplysia/A1_CNS_unsheared_vs_merged/tophat_out/accepted_hits.bam

# check header:
samtools view -H $BAM | wc
# 476578 (Trinity transcript count + 3)
#
samtools view -H $BAM | tail -1
#@PG	ID:TopHat	VN:2.0.8	CL:/usr/local/packages/tophat-2.0.8/tophat -p 8 --library-type fr-firststrand Trinity.trimmed.norm30 XUMTA_20150925_K00134_IL100062514_S32_L008_R1.fastq.gz XUMTA_20150925_K00134_IL100062514_S32_L008_R2.fastq.gz

# create smaller test files on /export to eliminate network delay as a factor:
setenv WDIR /export/biocode-test
mkdir -p $WDIR

# take first 100,000 and 1 million non-header lines:
samtools view -H $BAM > $WDIR/hits_partial_100K.sam
samtools view $BAM | head -100000  >> $WDIR/hits_partial_100K.sam

samtools view -H $BAM > $WDIR/hits_partial_1M.sam
samtools view $BAM | head -1000000 >> $WDIR/hits_partial_1M.sam

samtools view -bS $WDIR/hits_partial_100K.sam >$WDIR/hits_partial_100K.bam
samtools view -bS $WDIR/hits_partial_1M.sam >$WDIR/hits_partial_1M.bam
