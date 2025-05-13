#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=4
#SBATCH --mem=48g
#SBATCH --tmp=48g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This script will run minimap2, which has a wide variety of applications for long reads.

module load minimap2

# User provided input arguments
FASTA_FILE=/scratch.global/pmorrell/Selective_Sweeps/SSW_snps_final_morex_v2_renamed.fasta
REF_FILE=/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta.gz
OUT_PREFIX=SSW_SNPs_Morex_v2_minimap2
OUT_DIR=/scratch.global/pmorrell/Selective_Sweeps

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into reference dir
cd ${OUT_DIR}

# Align with minimap
# The '-a' option is for long reads against a reference
minimap2 -ax sr -t 3 ${REF_FILE} ${FASTA_FILE} > ${OUT_DIR}/${OUT_PREFIX}.sam
