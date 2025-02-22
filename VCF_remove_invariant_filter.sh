#!/bin/bash -l
#SBATCH --time=28:00:00
#SBATCH --ntasks=5
#SBATCH --mem=100g
#SBATCH --tmp=100g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

# Original script from: Jacob Pacheco

set -e
set -o pipefail

INPUT_DIR="/scratch.global/pmorrell/Inversions/WBDC355_10X_SNPS/test"
OUTPUT_DIR="/scratch.global/pmorrell/Inversions/WBDC355_10X_SNPS/test"
REFERENCE="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"

# Filter parameters for variants call from ~42x coverage Illumina PE data
MIN_MQ=20     # Minimum mapping quality
MIN_BQ=20     # Minimum base quality
MIN_DEPTH=5  # Minimum depth for variant calling
MAX_DEPTH=250 # Maximum depth to avoid regions with excessive coverage
MQSBZ=3     # Mapping Quality Strand Bias Z-score
RPBZ=3      # Read Position Bias Z-score    


# Filter parameters adjusted for ONT Nanopore data
#MIN_MQ=10     # Lower minimum mapping quality for Nanopore data
#$MIN_BQ=10     # Lower minimum base quality for Nanopore data
#MIN_DEPTH=3   # Minimum depth for variant calling
#MAX_DEPTH=100 # Maximum depth to avoid regions with excessive coverage
#MQSBZ=5       # Adjusted Mapping Quality Strand Bias Z-score
#RPBZ=5        # Adjusted Read Position Bias Z-score    


# Load required modules

module load bcftools_ML_2/1.20

mkdir -p "${OUTPUT_DIR}/filtered_results"

process_sample() {
    local VCF_FILE=$1
    local SAMPLE_NAME
    SAMPLE_NAME=$(basename "${VCF_FILE}" .vcf.gz)
    echo "Processing sample (chromosome): ${SAMPLE_NAME}"

    # create dir

    mkdir -p "${OUTPUT_DIR}/filtered_results/${SAMPLE_NAME}"
    local OUTPUT_PREFIX="${OUTPUT_DIR}/filtered_results/${SAMPLE_NAME}"

    # check if VCF exist
    if ! bcftools view -h "${VCF_FILE}" > /dev/null 2>&1; then
        echo "ERROR: ${VCF_FILE} is not a valid VCF file or is corrupted."
        return 1
    fi

    # Step 1: Remove invariant sites
    # The -c 1 flag ensures that only sites with at least one alternate allele are kept.
    # The --max-allele 2 flag ensures that only biallelic sites are kept.
    # The --exclude-types indels flag removes indels.

    echo "   -> Filtering for polymorphic, biallelic SNP sites"
    intermediate_vcf=$(mktemp)
    bcftools view -c 1 --max-allele 2 --exclude-types indels -Oz -o "${intermediate_vcf}" "${VCF_FILE}"
    bcftools index "${intermediate_vcf}"

     # Step 2: Apply additional filtering based on quality, depth, and strand bias
    echo "   -> Applying quality, depth, and strand bias filters..."
    bcftools filter -e "QUAL<${MIN_BQ} || MQ<${MIN_MQ} || INFO/DP<${MIN_DEPTH} || INFO/DP>${MAX_DEPTH} || INFO/MQSBZ>${MQSBZ} || INFO/RPBZ>${RPBZ}" \
        -Oz \
        -o "${OUTPUT_PREFIX}.filtered.vcf.gz" \
        "${intermediate_vcf}"
    bcftools index "${OUTPUT_PREFIX}.filtered.vcf.gz"

    # Clean up intermediate file
    rm "${intermediate_vcf}"
    echo "Finished processing ${SAMPLE_NAME}"
}

# Process all VCF files in the input directory (per chromosome)
echo "Looking for VCF files in ${INPUT_DIR}..."
find "${INPUT_DIR}" -name "*.vcf.gz" | while read -r VCF_FILE; do
    process_sample "${VCF_FILE}"
done

