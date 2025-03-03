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

INPUT_DIR="/scratch.global/pmorrell/Inversions/WBDC355_10X_SNPS/test/filtered_results"
OUTPUT_DIR="/scratch.global/pmorrell/Inversions/PlantCaduceus"
REFERENCE="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"
REFERENCE_INDEX="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta.fai"

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
    log "Processing VCF (chromosome): ${SAMPLE_NAME}"

    local intermediate_bed
    intermediate_bed=$(mktemp)
    zcat "${VCF_FILE}" | grep -v '#' | awk -v OFS='\t' '{print $1, $2-255, $2+257, $2, $4, $5}' >temp.bed

    log "   -> Generating intervals"
    local intermediate_bed2
    intermediate_bed2=$(mktemp)
    bedtools slop -i "${intermediate_bed}" -g "${REFERENCE_INDEX}" -l 254 -r 257 > "${intermediate_bed2}"

    log "   -> Generating contextual sequence"
    local intermediate_seq
    intermediate_seq=$(mktemp)
    bedtools getfasta -fi "${REFERENCE}" -bed "${intermediate_bed2}" -bedOut | cut -f 4  > temp.fas
    # cat "${intermediate_seq}"
    
    log "   -> Create the output file"
    local header="chr\tstart\tend\tpos\tref\talt\tsequences"
    echo -e "${header}" > "${OUTPUT_DIR}/${SAMPLE_NAME}_input.txt"
    #paste "${intermediate_bed}" "${intermediate_seq}" >> "${OUTPUT_DIR}/${SAMPLE_NAME}_input.txt"

    rm "${intermediate_bed}" "${intermediate_bed2}" "${intermediate_seq}"
}

log "Looking for VCF files in ${INPUT_DIR}..."
find "${INPUT_DIR}" -name "*.vcf.gz" | while read -r VCF_FILE; do
    process_vcfs "${VCF_FILE}"
done

log "All samples processed."
