#!/bin/bash -l
#SBATCH --time=28:00:00
#SBATCH --ntasks=5
#SBATCH --mem=100g
#SBATCH --tmp=100g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

# Original script Peter L. Morrell - St. Paul, MN

set -e
set -o pipefail


VCF01="/scratch.global/pmorrell/Inversions/WBDC355_10X_SNPS/test/filtered_results/WBDC355_10X_chr1H.filtered.vcf.gz"
VCF02="/scratch.global/pmorrell/Inversions/WBDC355_ONT_SNPS/test/filtered_results/WBDC355_ONT_chr1H.filtered.vcf.gz"
OUTPUT_DIR="/scratch.global/pmorrell/Inversions/"

# Load required modules
module load bcftools_ML_2/1.20

mkdir -p "${OUTPUT_DIR}/intersect"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

log "Starting to intersect VCF files"
bcftools isec -p "${OUTPUT_DIR}/intersect" "${VCF01}" "${VCF02}"
log "Finished intersecting VCF files"
