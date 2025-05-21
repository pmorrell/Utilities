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

# Peter L. Morrell - 17 April 2025 - St. Paul, MN

module load bedtools2/2.31.0-gcc-8.2.0-7j35k74

INPUT_DIR="/scratch.global/pmorrell/Selective_Sweeps"
OUTPUT_DIR="/scratch.global/pmorrell/Selective_Sweeps"
REFERENCE="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules.fasta.gz"
REFERENCE_INDEX="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules.fasta.gz.fai"
mkdir -p "${OUTPUT_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

process_vcfs() {
    log "   -> Read SNP information from VCF"
    local VCF_FILE="$1"
    local SAMPLE_NAME
    SAMPLE_NAME=$(basename "${VCF_FILE}" .vcf.gz)
    log "   -> Processing VCF (chromosome): ${SAMPLE_NAME}"

    local intermediate_bed
    intermediate_bed=$(mktemp)
    zcat < "${VCF_FILE}" | grep -v '#' | awk -v OFS='\t' '{print $1, $2-1, $2, $1 "_" $2}' > "${intermediate_bed}"

    log "   -> Generating intervals"
    local intermediate_bed2
    intermediate_bed2=$(mktemp)
    bedtools slop -i "${intermediate_bed}" -g "${REFERENCE_INDEX}" -l 61 -r 60 > "${intermediate_bed2}"

    log "   -> Generating contextual sequence output"
    local temp_fasta
    temp_fasta=$(mktemp)
    local final_output="${OUTPUT_DIR}/${SAMPLE_NAME}.fasta"
    
    # Generate FASTA to temp file first
    bedtools getfasta -fi "${REFERENCE}" -bed "${intermediate_bed2}" -nameOnly -fo "${temp_fasta}"
    
    # Process with sed to replace 61st position with N, then write directly to final output
    log "   -> Converting SNP position (61st base) to N"
    sed '/^>/!s/\(.\{60\}\)./\1N/' "${temp_fasta}" > "${final_output}"
    
    # Clean up temporary files
    rm -f "${intermediate_bed}" "${intermediate_bed2}" "${temp_fasta}"
    
    log "   -> Completed processing ${SAMPLE_NAME}"
}

log "   -> Looking for VCF files in ${INPUT_DIR}..."
find "${INPUT_DIR}" -name "*.vcf.gz" | while read -r VCF_FILE; do
    process_vcfs "${VCF_FILE}"
done

log "   -> All samples processed."

