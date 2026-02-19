#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --ntasks=5
#SBATCH --mem=24g
#SBATCH --tmp=24g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Peter L. Morrell - 19 April 2025 - Falcon Heights, MN


INPUT_DIR="/scratch.global/pmorrell/Selective_Sweeps"
OUTPUT_DIR="/scratch.global/pmorrell/Selective_Sweeps"
REFERENCE="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"
ILLUMINA="/scratch.global/pmorrell/Selective_Sweeps/SSW_snps_final_morex_v2_renamed_lookup.txt"
SNP_UTILS="/home/morrellp/pmorrell/shared/Software/SNP_Utils/snp_utils.py"

mkdir -p "${OUTPUT_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

process_sam() {
    
    log "   -> Read SNP information from SAM"
    local SAM_FILE="$1"
    local SAMPLE_NAME
    SAMPLE_NAME=$(basename "${SAM_FILE}" .sam)
    log "   -> Processing SAM file: ${SAMPLE_NAME}"

    log "   -> Finding variant positions, creating VCF"
    python3 "${SNP_UTILS}" SAM \
    --lookup "${ILLUMINA}" \
    --sam-file "${SAM_FILE}" \
    --reference "${REFERENCE}" \
    --outname "${OUTPUT_DIR}/${SAMPLE_NAME}.vcf"

    if [[ -f "${OUTPUT_DIR}/${SAMPLE_NAME}.vcf" ]]; then
    log "   -> Successfully created ${SAMPLE_NAME}.vcf"
        else
    log "   -> Failed to create ${SAMPLE_NAME}.vcf"
    fi


}

log "   -> Looking for SAM files in ${INPUT_DIR}..."
find "${INPUT_DIR}" -name "*.sam" | while read -r SAM_FILE; do
    process_sam "${SAM_FILE}"
done

log "   -> All samples processed."
