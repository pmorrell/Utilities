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

module load bedtools_ML/2.30.0

INPUT_DIR="/scratch.global/pmorrell/Inversions/"
OUTPUT_DIR="/scratch.global/pmorrell/Inversions/PlantCaduceus"
REFERENCE="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"

mkdir -p "${OUTPUT_DIR}/filtered_results"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

process_vcfs() {
    log "   -> Read SNP information from VCF"
    local VCF_FILE="$1"
    local SAMPLE_NAME
    SAMPLE_NAME=$(basename "${VCF_FILE}" .vcf.gz)
    log "Processing VCF (chromosome): ${SAMPLE_NAME}"

    local intermediate_bed
    intermediate_bed=$(mktemp)
    zcat "${VCF_FILE}" | grep -v '#' | awk -v OFS='\t' '{print $1, $2-254, $2+257, $2, $3, $4}' > "${intermediate_bed}"

    log "   -> Generating intervals"
    local intermediate_bed2
    intermediate_bed2=$(mktemp)
    bedtools slop -i "${intermediate_bed}" -g Barley_MorexV3_pseudomolecules.txt -l 254 -r 257 > "${intermediate_bed2}"

    log "   -> Generating contextual sequence"
    local intermediate_bed3
    intermediate_bed3=$(mktemp)
    bedtools getfasta -fi "${REFERENCE}" -bed "${intermediate_bed2}" -bedOut > "${intermediate_bed3}"

    local header="chr\tstart\tend\tpos\tref\talt\tsequences"
    echo -e "${header}" > "${OUTPUT_DIR}/${SAMPLE_NAME}_input.txt"
    paste "${intermediate_bed}" "${intermediate_bed3}" >> "${OUTPUT_DIR}/${SAMPLE_NAME}_input.txt"

    rm "${intermediate_bed}" "${intermediate_bed2}" "${intermediate_bed3}"
}

log "Looking for VCF files in ${INPUT_DIR}..."
find "${INPUT_DIR}" -name "*.vcf.gz" | while read -r VCF_FILE; do
    process_vcfs "${VCF_FILE}"
done

log "All samples processed."
