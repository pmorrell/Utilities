#!/bin/bash -l
#SBATCH --time=12:00:00
#SBATCH --ntasks=4
#SBATCH --mem=400g
#SBATCH --tmp=400g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# IMPORTANT: Have to load the appropriate conda environment
# `conda activate svimasm_env`

# The application and associated tools are all in the directory below
WORKING_DIR=/scratch.global/pmorrell/Morex_v3/svim_asm
REFERENCE_FASTA=/scratch.global/pmorrell/Morex_v3/HvulgareMorex_702_V3.hardmasked_chr1H.fa.gz
# Location of a sorted and indexed BAM file from assembly alignment
BAM=/scratch.global/pmorrell/Morex_v3/out_dir/OUH602_chr1H_asm10.bam

# Check if our dir exists, if not make it
mkdir -p ${WORKING_DIR}

# Go into output directory
cd ${WORKING_DIR}

# Run the svim-asm variant caller
# This should produce a VCF with up to 6 classes of SVs
svim-asm haploid ${WORKING_DIR} ${BAM} ${REFERENCE_FASTA}
