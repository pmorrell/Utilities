#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=5
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# Peter L. Morrell - 19 April 2025 - Falcon Heights, MN

OUTPUT_DIR="/scratch.global/pmorrell/Selective_Sweeps"
CONTEXT_FASTA="/scratch.global/pmorrell/Selective_Sweep/SSW_snps_final_morex_v2_renamed.fasta"
VCF="/scratch.global/pmorrell/Selective_Sweeps/SSW_snps_final_morex_v2_renamed.vcf.gz"
ILLUMINA_LOOKUP="/scratch.global/pmorrell/Selective_Sweeps/Illumina_lookup.py"

mkdir -p "${OUTPUT_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

    SAMPLE_NAME=$(basename "${VCF}" .vcf.gz)
    log "   -> Processing VCF file: ${SAMPLE_NAME}"

    log "   -> Running Illumina_lookup.py"
    python3 "${ILLUMINA_LOOKUP}" --fasta "${CONTEXT_FASTA}" \
    --vcf "${VCF}" \
    --output "${OUTPUT_DIR}/${SAMPLE_NAME}_lookup.txt"


    if [[ -f "${OUTPUT_DIR}/${SAMPLE_NAME}_lookup.txt" ]]; then
    log "   -> Successfully created ${SAMPLE_NAME}_lookup.txt"
        else
    log "   -> Failed to create ${SAMPLE_NAME}_lookup.txt"
    fi

    log "   -> All samples processed."
