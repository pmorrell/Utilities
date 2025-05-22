#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --ntasks=5
#SBATCH --mem=24g
#SBATCH --tmp=24g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Peter L. Morrell - 23 February 2025 - St. Paul, MN

module load bedtools2/2.31.0-gcc-8.2.0-7j35k74

INPUT_DIR="/scratch.global/pmorrell/Barley_LLMs"
OUTPUT_DIR="/scratch.global/pmorrell/Barley_LLMs"
REFERENCE="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"
REFERENCE_INDEX="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta.fai"

mkdir -p "${OUTPUT_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

process_vcfs() {
    log "   -> Processing VCF file"
    local VCF_FILE="$1"
    local SAMPLE_NAME
    SAMPLE_NAME=$(basename "${VCF_FILE}" .vcf.gz)
    log "Processing VCF (chromosome): ${SAMPLE_NAME}"

    # First use bedtools slop directly on VCF to extend regions
    log "  -> Extracting SNP positions from the VCF"
    local SNP_positions=$(mktemp)
    zgrep -v '#' "${VCF_FILE}" | awk -v OFS='\t' '{print $1, $2-1, $2, $2-1, $4, $5}' > "${SNP_positions}"

    log "   -> Extending regions to create 512bp windows with SNP at position 256"
    local extended_regions=$(mktemp)
    bedtools slop -i "${SNP_positions}" -g "${REFERENCE_INDEX}" -l 256 -r 255 > "${extended_regions}"
    
    # Now use bedtools getfasta to extract sequences
log "   -> Generating contextual sequence"
    local intermediate_seq
    intermediate_seq=$(mktemp)
    bedtools getfasta -fi "${REFERENCE}" -bed "${extended_regions}" -bedOut > "${intermediate_seq}"
    
    log "   -> Create the output file"
    local header="chr\tstart\tend\tpos\tref\talt\tsequences"
    echo -e "${header}" > "${OUTPUT_DIR}/${SAMPLE_NAME}_input.txt"
    cat "${intermediate_seq}" >> "${OUTPUT_DIR}/${SAMPLE_NAME}_input.txt"

    # Clean up temporary files
}

log "Looking for VCF files in ${INPUT_DIR}..."
find "${INPUT_DIR}" -name "*.vcf.gz" | while read -r VCF_FILE; do
    process_vcfs "${VCF_FILE}"
done

log "All samples processed."

