#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=8
#SBATCH --mem=200g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

set -e
set -o pipefail

# This script will run minimap2, which has a wide variety of applications for long reads.

MINIMAP2=/panfs/roc/groups/9/morrellp/pmorrell/Apps/HLi/minimap2/minimap2

# User provided input arguments
FASTA_FILE=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta
REF_FILE=/panfs/roc/groups/9/morrellp/shared/Datasets/10x_Genomics/Barley/WBDC355_Assembly_2018-04-28/barley_outputs_barley.pseudohap.fasta.gz
OUT_PREFIX=WBDC355_comparisons
OUT_DIR=/scratch.global/pmorrell/WBDC355_10xGenomics

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into reference dir
cd ${OUT_DIR}

# Align with minimap
# The '-a' option is for long reads against a reference
$MINIMAP2 -ax asm10 ${REF_FILE} --cs ${FASTA_FILE} > ${OUT_DIR}/${OUT_PREFIX}.sam
