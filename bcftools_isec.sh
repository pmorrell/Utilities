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
# Modified to accept command-line arguments

set -e
set -o pipefail

# Usage function
usage() {
    echo "Usage: $0 <FILE1> <FILE2> <OUTPUT_DIR>"
    echo ""
    echo "Arguments:"
    echo "  FILE1       Path to first VCF/BCF file (compressed or uncompressed)"
    echo "  FILE2       Path to second VCF/BCF file (compressed or uncompressed)"
    echo "  OUTPUT_DIR  Directory for output files"
    echo ""
    echo "Supported formats: .vcf, .vcf.gz, .bcf"
    echo ""
    echo "Example:"
    echo "  $0 file1.vcf.gz file2.bcf output_dir"
    exit 1
}

# Check if correct number of arguments provided
if [ "$#" -ne 3 ]; then
    echo "Error: Incorrect number of arguments"
    usage
fi

FILE01="${1}"
FILE02="${2}"
OUTPUT_DIR="${3}"

# Validate input files exist
if [ ! -f "${FILE01}" ]; then
    echo "Error: File '${FILE01}' not found"
    exit 1
fi

if [ ! -f "${FILE02}" ]; then
    echo "Error: File '${FILE02}' not found"
    exit 1
fi

# Load required modules
module load bcftools_ML_2/1.20

mkdir -p "${OUTPUT_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%%Y-%%m-%%d %%H:%%M:%%S') - ${msg}"
}

log "Processing files:"
log "  FILE1: ${FILE01}"
log "  FILE2: ${FILE02}"
log "  Output: ${OUTPUT_DIR}"

log "Creating indexes for VCF/BCF files"
bcftools index -f "${FILE01}"
bcftools index -f "${FILE02}"

log "Starting to intersect VCF/BCF files"
bcftools isec -p "${OUTPUT_DIR}" "${FILE01}" "${FILE02}"
log "Finished intersecting VCF/BCF files"