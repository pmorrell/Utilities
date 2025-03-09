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
module load apptainer

mkdir -p ~/Apps

singularity build --force\
   ~/Apps/ensembl-vep-r_113.3.sif \
   docker://ensemblorg/ensembl-vep:release_113.3

# Variant sets should be either 'all', deletions', 'insertions', 'snps', 'common', or 'rare'
VARIANT_SET=all
SPECIES=barley

# Add path to bgzipped and indexed .vcf file 
VCF="/scratch.global/pmorrell/Inversions/VeP/WBDC355_10X_chr1H.filtered.vcf.gz"
# Add path to bgzipped and indexed .gff file 
GFF="/scratch.global/pmorrell/Inversions/VeP/Hv_Morex.pgsb.Jul2020.sorted.gff3.gz"
# Add path to bgzipped and indexed reference fasta 
FASTA="/scratch.global/pmorrell/Inversions/VeP/Barley_MorexV3_pseudomolecules.fasta.gz"
OUTPUT_PREFIX="$(basename "${VCF}" .vcf.gz)"

# Specify out directory
OUT="/scratch.global/pmorrell/Inversions/VeP"
cd "${OUT}"

singularity run --bind ${HOME}/:${HOME}/ ${HOME}/Apps/ensembl-vep-r_113.3.sif \
    vep \
    -i "${VCF}" \
    --gff "${GFF}" \
    --fasta "${FASTA}" \
    --species "${SPECIES}" \
    --total_length \
    --check_svs \
    --verbose \
    --format vcf \
    --force \
    --warning_file "${OUT}/${OUTPUT_PREFIX}_${VARIANT_SET}_warnings.txt" \
    -o "${OUT}/${OUTPUT_PREFIX}_${VARIANT_SET}.txt"
