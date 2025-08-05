#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=17
#SBATCH --mem=400g
#SBATCH --tmp=400g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

module load samtools/1.9

# This script will sort a BAM file and convert the file from SAM to BAM format

# User provided input arguments
SAM_FILE=/scratch.global/pmorrell/Inversions/WBDC355_10X.sam
OUT_DIR=/scratch.global/pmorrell/Inversions/
BAM_FILE=$(basename "$SAM_FILE" .sam)

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into reference dir
cd ${OUT_DIR}

# Running samtools to sort and convert
# This setup from SVIM-asm, but basic setup should work generally
samtools sort -m10G -@16 -o ${BAM_FILE}.bam ${SAM_FILE}

# Index the BAM file
samtools index -c ${BAM_FILE}.bam 
