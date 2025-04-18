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

INPUT_DIR="/scratch.global/pmorrell/Selective_Sweeps"
OUTPUT_DIR="/scratch.global/pmorrell/Selective_Sweeps"
REFERENCE="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules.fasta.gz"
SLOP="60"
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
    zcat < "${VCF_FILE}" | grep -v '#' | awk -v OFS='\t' '{print $1, $2-61, $2+60, $1_$2}' > "${intermediate_bed}"

    #log "   -> Generating intervals"
    #local intermediate_bed2
    #intermediate_bed2=$(mktemp)
    #bedtools slop -i "${intermediate_bed}" -g "${REFERENCE_INDEX}" -b "${SLOP}" > "${intermediate_bed2}"

    log "   -> Generating contextual sequence output"
    local intermediate_seq
    intermediate_seq=$(mktemp)
    bedtools getfasta -fi "${REFERENCE}" -bed "${intermediate_bed2}" -name -fo "${OUT_DIR}/${SAMPLE_NAME}.fasta" 
}

log "Looking for VCF files in ${INPUT_DIR}..."
find "${INPUT_DIR}" -name "*.vcf.gz" | while read -r VCF_FILE; do
    process_vcfs "${VCF_FILE}"
done

log "All samples processed."
