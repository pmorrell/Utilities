#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --ntasks=5
#SBATCH --mem=24g
#SBATCH --tmp=24g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

# Original script Peter L. Morrell - St. Paul, MN

set -e
set -o pipefail


VCF01="/scratch.global/pmorrell/Inversions/WBDC355_10X_SNPS/test/WBDC355_10X_chr1H.vcf.gz"
VCF02="/scratch.global/pmorrell/Inversions/BOPA_WBDC355.vcf.gz"
OUTPUT_DIR="/scratch.global/pmorrell/Inversions/BOPA_intersect"

# Load required modules
module load bcftools_ML_2/1.20

mkdir -p "${OUTPUT_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}


log "Create an index for the VCFs"
bcftools index -f "${VCF01}"
bcftools index -f "${VCF02}"

log "Starting to intersect VCF files"
bcftools isec -p "${OUTPUT_DIR}" "${VCF01}" "${VCF02}"
log "Finished intersecting VCF files"
