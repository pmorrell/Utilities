#!/bin/bash -l
#SBATCH --time=3:00:00
#SBATCH --ntasks=2
#SBATCH --mem=200g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# `paftools.js` processes alignment files

# Load the other tools needed
k8=/panfs/roc/groups/9/morrellp/pmorrell/bin/k8-0.2.5/k8-Linux
PAFTOOLS=/panfs/roc/groups/9/morrellp/shared/Software/minimap2/misc/paftools.js

# The application and associated tools are all in the directory below
OUT_DIR=/scratch.global/pmorrell/Morex_v3/out_dir
REFERENCE_FASTA=/scratch.global/pmorrell/Morex_v3/Morex_V3_chr5H.fasta.gz

# Files and settings
PAF=/scratch.global/pmorrell/Morex_v3/out_dir/OUH602_chr5H.paf
SAMPLE_NAME=OUH602_Chr5H
threads=2
MIN_LENGTH_COVERAGE=5000
MIN_LENGTH_VARIANTS=25000

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into output directory
cd ${OUT_DIR}

# A minimum set of positions for calling variants. Not including an "excludes" file with problematic regions
sort -k6,6 -k8,8n ${PAF} | \
${k8} ${PAFTOOLS} call -f ${REFERENCE_FASTA} -s ${SAMPLE_NAME} -l ${MIN_LENGTH_COVERAGE} -L ${MIN_LENGTH_VARIANTS} \
> ${SAMPLE_NAME}.vcf
