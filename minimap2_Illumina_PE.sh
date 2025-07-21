#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=490g
#SBATCH --tmp=490g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

set -e
set -o pipefail

# This script will run minimap2, which has a wide variety of applications for long reads.

module load minimap2/2.26 

# User provided input arguments
FASTQ_01=/scratch.global/pmorrell/Inverions/WBDC355_10X_Genomics/WBDC355_R1_001_trim.fastq.gz
FASTQ_02=/scratch.global/pmorrell/Inverions/WBDC355_10X_Genomics/WBDC355_R2_001_trim.fastq.gz
REF_FILE=/home/morrellp/pmorrell/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta
OUT_PREFIX=WBDC355_10X_SNPS
OUT_DIR=/scratch.global/pmorrell/Inverions/WBDC355_10X_SNPS

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into reference dir
cd ${OUT_DIR}

# Align with minimap2
# Using options for Illumina PE reads with 
# Using a larger index so we get a SAM header
minimap2 -ax rs -t 8 -I6g ${REF_FILE} ${FASTQ_01} ${FASTQ_02} > ${OUT_DIR}/${OUT_PREFIX}.sam
