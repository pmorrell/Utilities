#!/bin/bash -l
#SBATCH --time=28:00:00
#SBATCH --ntasks=5
#SBATCH --mem=100g
#SBATCH --tmp=100g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

module load bcftools_ML_2/1.20

# Path to the reference genome
MOREX=/home/morrellp/pmorrell/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta
# Path to ONT BAM file
BAM=/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/nanopore/WBDC_355/WBDC_combined_1_10/WBDC_355_v3/bam_file/WBDC_355_1_10_v3_sorted.bam
OUT_DIR=/scratch.global/pmorrell/Inversions/WBDC355_ONT_SNPS
OUT_PREFIX=WBDC355_ONT_chr7H
REGION=chr7H

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Call variants on one chromosome using bcftools mpileup presets for ONT data
bcftools mpileup --regions $REGION -X ont-sup-1.20 $BAM -f $MOREX | bcftools call -mv -Oz -o ${OUT_DIR}/${OUT_PREFIX}.vcf.gz
