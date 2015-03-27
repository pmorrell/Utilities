#!/bin/env bash
set -e
set -u
set -o pipefail

#PBS -l mem=12gb,nodes=1:ppn=1,walltime=4:00:00 
#PBS -m abe 
#PBS -M pmorrell@umn.edu 
#PBS -q lab

#module load bedtools
BEDTOOLS=~/Apps/Quinlan/bedtools2/bin/bedtools

WORKING=~/scratch
PROJECT='SCN'
SAMPLE='T7013_T5009'
DATE='2015-03-12'

# BAM file to consider
BAM=${WORKING}/${PROJECT}/${SAMPLE}/Sample_${SAMPLE}_${PROJECT}_${DATE}_Finished.bam
# BED file with barley exome coverage
TARGETS=~/Shared/Datasets/Annotations/Soybean_Cyst_Nematode/Heterodera_glycines_OP25_gene_models.gff

#   check if R is installed and in the path
if `command -v Rscript > /dev/null 2> /dev/null`
    then 
        echo "R is installed, OK"
    else
        echo "You need R (Rscript) to be installed and in your \$PATH"
        exit 1
fi


cd $WORKING

mkdir -p ${SCN}/${SAMPLE}

${BEDTOOLS} coverage \
-hist \
-abam ${BAM} \
-b ${TARGETS} \
> Sample_${SAMPLE}_${PROJECT}_${DATE}.coverage.hist.txt
