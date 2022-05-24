#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=3
#SBATCH --mem=400g
#SBATCH --tmp=400g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

set -e
set -o pipefail

# This script will run minimap2, which has a wide variety of applications for long reads.

MINIMAP2=/panfs/roc/groups/9/morrellp/pmorrell/Apps/HLi/minimap2/minimap2

# User provided input arguments
FASTA_FILE=/scratch.global/pmorrell/OUH602/OUH602_chr1H.fasta.gz
REF_FILE=/scratch.global/pmorrell_Morex_v3/HvulgareMorex_702_V3.hardmasked_chr1H.fa.gz
OUT_PREFIX=OUH602_hardmask_chr1H
OUT_DIR=/scratch.global/pmorrell/Morex_v3/out_dir

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into reference dir
cd ${OUT_DIR}

# Align with minimap
# The '-a' option is for long reads against a reference
$MINIMAP2 -a -x asm5 --cs -r2k -t 3 ${REF_FILE} ${FASTA_FILE} > ${OUT_DIR}/${OUT_PREFIX}.sam
