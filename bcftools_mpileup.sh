#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -euo pipefail

# Load the latest bcftools
module load bcftools_ML_2/1.20

# Path to the reference genome
REF_FILE=/home/morrellp/pmorrell/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_plastids.fasta
# Path to ONT BAM file
BAM_LIST=/scratch.global/pmorrell/WBDC_resequencing/WBDC_BAM_list.txt
OUT_DIR=/scratch.global/pmorrell/WBDC_resequencing
OUT_PREFIX=WBDC_PacBio
REGION=chr1H

# Check if our dir exists, if not make it
mkdir -p "$OUT_DIR"

# Validate input files
if [[ ! -f "$REF_FILE" ]]; then
    echo "Error: Reference file not found: $REF_FILE"
    exit 1
fi

if [[ ! -f "$BAM_LIST" ]]; then
    echo "Error: BAM list file not found: $BAM_LIST"
    exit 1
fi

echo "Starting variant calling with $(wc -l < "$BAM_LIST") samples"

# Call variants using bcftools mpileup presets for ONT data
bcftools mpileup -Ou -f "$REF_FILE" --regions "$REGION" --threads 15 -q 20 -Q 20 -a FORMAT/AD,FORMAT/DP -b "$BAM_LIST" | \
bcftools call -mv -Oz --threads 15 -o "${OUT_DIR}/${OUT_PREFIX}_${REGION}.vcf.gz"

# Index the VCF file
bcftools index --csi "${OUT_DIR}/${OUT_PREFIX}_${REGION}.vcf.gz"

echo "Variant calling completed. Output: ${OUT_DIR}/${OUT_PREFIX}_${REGION}.vcf.gz"
