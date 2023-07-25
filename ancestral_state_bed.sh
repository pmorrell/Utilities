#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

set -e
set -o pipefail

module load bedtools/2.29.2

# This script extracts SNP positions from a VCF then uses those positions to create a bed file
# of ancestral states.

# User provided input arguments
OUT_PREFIX=GBS_bulbosum
VCF=/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/Ahmad_GBS_morex_v3/final_filtered_treatmissing_biallilic_all_chrs_mafgt0.0016.recode.vcf.gz 
REF_FILE=/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta
OUT_DIR=/panfs/jay/groups/9/morrellp/pmorrell/Workshop/Inversions_SNPs/Ancestral_state

# Check if our directory exists, if not make it
mkdir -p ${OUT_DIR}

# Go into out directory
cd ${OUT_DIR}

# Create a BED file of SNP positions in the VCF file
BED_IN=$(zgrep -v '#' $VCF | awk '{FS="\t";OFS="\t";print $1,$2-1,$2}') 


bedtools getfasta -fi "$REF_FILE" -bed "$BED_IN" -bedOut ${OUT_DIR}/${OUT_PREFIX}_SNPS.bed

