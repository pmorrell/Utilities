#!/bin/bash -l
#SBATCH --time=12:00:00
#SBATCH --ntasks=5
#SBATCH --mem=16g
#SBATCH --tmp=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err


# This script runs Variant effect Predictor (VeP) on a set of variants
# It requires a VCF, GFF, and fasta file as input

#Dependencies
module load htslib/1.9  
module load perl/modules.centos7.5.26.1

VEP="/panfs/jay/groups/9/morrellp/shared/Software/ensembl-vep-release-108.1/vep"

# Variant sets should be either 'all', deletions', 'insertions', 'snps', 'common', or 'rare'
VARIANT_SET=all
SPECIES=barley

# Add path to .vcf file 
VCF="/scratch.global/pmorrell/Inversions/WBDC355_10X_SNPS/test/filtered_results/WBDC355_10X_chr1H.filtered.vcf.gz"
# Add path to .gff file 
GFF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.sorted.gff3"
# Add path to reference fasta 
FASTA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"
OUTPUT_PREFIX="${VCF%.vcf.gz}"

# Specify out directory
OUT="/scratch.global/pmorrell/Inversions/WBDC355_10X_SNPS/
cd "${OUT}"

"${VEP}" \
    -i "${VCF}" \
    --gff "${GFF}" \
    --fasta "${FASTA}" \
    --species "${SPECIES}" \
    --total_length \
    --check_svs \
    --verbose \
    --format vcf \
    --force \
    --warning_file "${OUTPUT_PREFIX}_${VARIANT_SET}.txt" \
    -o "${OUTPUT_PREFIX}_${VARIANT_SET}.txt"