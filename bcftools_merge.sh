#!/bin/bash -l
#SBATCH --time=6:00:00
#SBATCH --ntasks=6
#SBATCH --mem=48g
#SBATCH --tmp=48g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Peter L. Morrell - 14 February 2025 - Falcon Heights, MN

module load bcftools/1.21

BASE_NAME=UG100
WORK_DIR=/scratch.global/pmorrell/Cowpea/UG100
cd ${WORK_DIR}

bcftools merge *.vcf.gz \
    --threads 5 \
    -O b \
    -o ${BASE_NAME}.bcf
bcftools index ${BASE_NAME}.bcf

