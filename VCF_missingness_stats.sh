#!/bin/bash -l
#SBATCH –time=10:00:00
#SBATCH –ntasks=1
#SBATCH –cpus-per-task=10
#SBATCH –mem=200g
#SBATCH –tmp=200g
#SBATCH –mail-type=ALL
#SBATCH –mail-user=your_email@umn.edu
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -euo pipefail

log() {
local msg=”$1”
echo “$(date +’%Y-%m-%d %H:%M:%S’) - ${msg}”
}

# This script summarizes the variation in VCFs from resequencing
# It calculates statistics on missing variants for each VCF and produces bcftools
# stats and associated plots
# The script assumes that the user has run `bcftools +fill-tags` to generate the F_MISSING tag
# Configure paths (can be overridden by environment variables)

WORK_DIR=”${WORK_DIR:-/scratch.global/pmorrell/WBDC_resequencing}”
OUTPUT_DIR=”${OUTPUT_DIR:-${WORK_DIR}/variation_summary}”

module load bcftools_ML_2/1.20
module load datamash_ML/1.3
module load python3/3.8.3_anaconda2020.07_mamba
module load texlive/20131202


# Check if required modules are available

if ! module list 2>&1 | grep -q bcftools; then
log “Error: bcftools module not loaded properly”
exit 1
fi

if ! module list 2>&1 | grep -q datamash; then
log “Error: datamash module not loaded properly”
exit 1
fi

cd “$WORK_DIR”
mkdir -p “$OUTPUT_DIR”

process_vcfs() {
log “   -> Processing VCF file”
local VCF_FILE=”$1”
local SAMPLE_NAME
SAMPLE_NAME=$(basename “${VCF_FILE}” .vcf.gz)
log “Processing VCF (chromosome): ${SAMPLE_NAME}”

```
# Create sample-specific output directory
mkdir -p "${OUTPUT_DIR}/${SAMPLE_NAME}"

# Extract the chromosome name from the VCF file name
# Calculate missing - GATK indels
log "   -> Extracting missingness for each variant"
bcftools query -f '%F_MISSING\n' "$VCF_FILE" | \
    datamash --round 3 mean 1 sstdev 1 median 1 min 1 max 1 count 1 > "${OUTPUT_DIR}/${SAMPLE_NAME}/missing_summary.txt"

# Generate summary statistics for the VCF file
log "   -> Generating summary statistics for ${SAMPLE_NAME}"
bcftools stats "$VCF_FILE" > "${OUTPUT_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_stats.txt"

# Option 1: Use bcftools plot-vcfstats (requires plot-vcfstats.py and dependencies)
log "   -> Generating plots with plot-vcfstats"
if command -v plot-vcfstats &> /dev/null; then
    plot-vcfstats -p "${OUTPUT_DIR}/${SAMPLE_NAME}/plots_" "${OUTPUT_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_stats.txt"
else
    log "   -> plot-vcfstats not found, skipping automatic plotting"
fi
```

}

# Main execution logic

main() {
log “Starting VCF processing pipeline”

```
# Check if VCF files argument is provided
if [ $# -eq 0 ]; then
    log "Usage: $0 <vcf_file1> [vcf_file2] ..."
    log "Or set VCF_FILES environment variable with space-separated file paths"
    exit 1
fi

# Process each VCF file
for VCF_FILE in "$@"; do
    if [ -f "$VCF_FILE" ]; then
        process_vcfs "$VCF_FILE"
    else
        log "Warning: VCF file not found: $VCF_FILE"
    fi
done

log "VCF processing pipeline completed"
```

}

# Execute main function with all command line arguments

main “$@”