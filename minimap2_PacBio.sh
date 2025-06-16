#!/bin/bash -l
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=200g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -euo pipefail

#module load minimap2/2.26
MINIMAP2=/users/6/pmorrell/Apps/HLi/minimap2/minimap2

# User provided input arguments
REF_FILE=/home/morrellp/pmorrell/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta
OUT_PREFIX=WBDC207_PacBio
OUT_DIR=/scratch.global/pmorrell/WBDC_resequencing/WBDC207/PacBio
FASTQ_LIST=/scratch.global/pmorrell/WBDC_resequencing/WBDC207/PacBio/WBDC207_trim_fastq.txt  # Path to a file containing a list of FASTQ files

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into reference dir
cd ${OUT_DIR}

# Align with minimap2 using a loop
while read FASTQ_FILE; do
    $MINIMAP2 -ax map-hifi -t 8 -I6g ${REF_FILE} ${FASTQ_FILE} >> ${OUT_DIR}/${OUT_PREFIX}.sam
done < ${FASTQ_LIST}
