#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

set -e
set -o pipefail

module load bwa/0.7.17

# This script does a quick alignment as a check for the reference genome

# User provided input arguments
FASTA_FILE=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/iSelect_all.fas
REF_FILE=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/Vunguiculata_IT97K-499-35_v1.0.fa
OUT_PREFIX=iSelect_cowpea_BWA
OUT_DIR=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/bwa

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into reference dir
cd ${OUT_DIR}

# Align with bwa mem
bwa mem ${REF_FILE} ${FASTA_FILE} > ${OUT_DIR}/${OUT_PREFIX}.sam
