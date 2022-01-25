#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=8
#SBATCH --mem=200g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

set -e
set -o pipefail

module load samtools

# User provided input arguments
FASTA_FILE=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta
BAM_FILE=/panfs/roc/groups/9/morrellp/shared/Projects/WBDC_inversions/Compare_de_novo/WBDC355_comparisons.bam
OUT_PREFIX=WBDC355_comparisons_MD
OUT_DIR=/panfs/roc/groups/9/morrellp/shared/Projects/WBDC_inversions/Compare_de_novo

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into reference dir
cd ${OUT_DIR}

# Align with minimap
# The '-a' option is for long reads against a reference
samtools calmd -b --threads 4 ${BAM_FILE} --reference ${FASTA_FILE} > ${OUT_DIR}/${OUT_PREFIX}.bam
