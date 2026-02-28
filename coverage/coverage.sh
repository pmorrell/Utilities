#!/bin/env sh

WORKING=~/scratch/coverage

# BAM file to consider
BAM=~/shared/Datasets/NGS/Alignments/Deleterious_Mutations/Morex_KU-Missouri_2014-01-02_Finished_Realigned.bam
# BED file with barley exome coverage
TARGETS=~/shared/Datasets/NGS/Alignments/Deleterious_Mutations/Exome_Capture_Regions.bed

cd $WORKING

bedtools coverage \
-hist \
-abam ${BAM} \
-b ${TARGET} \
> ${BAM}.coverage.hist.txt
