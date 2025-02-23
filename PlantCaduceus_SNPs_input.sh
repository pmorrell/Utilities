#!/bin/bash -l
#SBATCH --time=28:00:00
#SBATCH --ntasks=5
#SBATCH --mem=100g
#SBATCH --tmp=100g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

# Peter L. Morrell - St. Paul, MN  -  23 March 2025

set -e
set -o pipefail

# Load required modules
module load bedtools_ML/2.30.0

INPUT_DIR="/scratch.global/pmorrell/Inversions/WBDC355_10X_SNPS/test"
OUTPUT_DIR="/scratch.global/pmorrell/Inversions/WBDC355_10X_SNPS/test"
INTERVAL=1000 # Interval to extend gene positions
GFF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/MorexV3.1-2019-05-30.gff3"
REFERENCE="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"


mkdir -p "${OUTPUT_DIR}/filtered_results"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}


process_genes() {
    # Step 1: Turn a GFF file into a BED file of gene positions
    log "   -> Convert a GFF of gene positions to a BED file"
    intermediate_bed=$(mktemp)
    awk '$3 == "gene" {print $1"\t"$4-1"\t"$5}' "${GFF}" > "${intermediate_bed}"
    intermediate_bed2=$(mktemp)
    bedtools slop -i "${intermediate_bed}" -g "${REFERENCE}" -b "${INTERVAL}" > "${intermediate_bed2}"
    
    # Step 2: Identify SNPs that overlap with gene intervals
    log "   -> Partition VCF for gene intervals only"
    local VCF_FILE=$1
    local SAMPLE_NAME
    SAMPLE_NAME=$(basename "${VCF_FILE}" .vcf.gz)
    log "Processing VCF (chromosome): ${SAMPLE_NAME}"
    bcftools index "${VCF_FILE}"
    bcftools view -R "${intermediate_bed2}" -Oz -o "${OUTPUT_PREFIX}.genes.vcf.gz" "${VCF_FILE}"
    bcftools index "${OUTPUT_PREFIX}.genes.vcf.gz"
}

# Process all VCF files in the input directory (per chromosome)
log "Looking for VCF files in ${INPUT_DIR}..."
find "${INPUT_DIR}" -name "*.vcf.gz" | while read -r VCF_FILE; do
    process_sample "${VCF_FILE}"
done

log "All samples processed."


    # Clean up intermediate files
    rm "${intermediate_vcf}"
    log "Finished processing ${SAMPLE_NAME}"