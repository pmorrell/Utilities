#!/bin/bash -l
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --tmp=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -euo pipefail

# This script runs Variant effect Predictor (VeP) on a set of variants
# It requires a VCF, GFF, and fasta file as input

#Dependencies
module load apptainer

######################
# Set working directory
WORKING_DIR="/scratch.global/pmorrell/Barley_LLMs"
cd ${WORKING_DIR}

singularity build -F \
   ${WORKING_DIR}/vep.sif \
   docker://ensemblorg/ensembl-vep:release_114.1

######################
# Set VeP variables

# Variant sets should be either 'all', 'deletions', 'insertions', 'snps', 'common', or 'rare'
VARIANT_SET="snps"

# Add path to .gff file. Note: it should be sorted.
#(grep ^"#" ${WORKING_DIR}/MorexV3_chrNH.gff; grep -v ^"#" ${WORKING_DIR}/MorexV3_chrNH.gff | sort -k1,1 -k4,4n) > ${WORKING_DIR}/MorexV3_chrNH.sorted.gff3
#bgzip ${WORKING_DIR}/MorexV3_chrNH.sorted.gff3
#tabix -f -p gff -C ${WORKING_DIR}/MorexV3_chrNH.sorted.gff3.gz
GFF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotationHv_Morex.pgsb.Jul2020.sorted.gff3.gz"

# Add path to reference fasta
FASTA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta.gz"

# Set species
SPECIES="hordeum_vulgare"

# Add path to .vcf file
VCF="${WORKING_DIR}/hybrid13_snps_biallelic.callable.vcf.gz"

singularity run --bind ${WORKING_DIR}/:${WORKING_DIR}/ ${WORKING_DIR}/vep.sif \
    vep \
       -i "${VCF}" \
       --gff "${GFF}" \
       --fasta "${FASTA}" \
       --species "${SPECIES}" \
       --database \
       --genomes \
       --total_length \
       --check_svs \
       --verbose \
       --format vcf \
       --force \
       --fork 100 \
       --warning_file "callable/Morex_${VARIANT_SET}_warning.txt" \
       -o "callallable/Morex_${VARIANT_SET}.txt"
