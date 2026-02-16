#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=200g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@umn.edu
#SBATCH -o %A.out
#SBATCH -e %A.err

set -euo pipefail

# --- Setup & Logging ---
log() {
  echo "$(date +'%Y-%m-%d %H:%M:%S') - $1"
}

# 1. Argument Validation
if [ $# -lt 2 ]; then
    log "Error: Missing arguments."
    echo "Usage: sbatch $0 <input_file> <work_dir>"
    echo "Example: sbatch $0 my_data.bcf /scratch.global/pmorrell/Cowpea/Project1"
    exit 1
fi

INPUT_VCF="$1"
WORK_DIR="$2"
OUTPUT_DIR="${WORK_DIR}/variation_summary"
DATAMASH=/projects/standard/morrellp/shared/Software/datamash-1.9/datamash

# 2. Environment Setup
module load bcftools/1.21
module load python3/3.8.3_anaconda2020.07_mamba
module load texlive/20131202

# Move to the working directory
cd "$WORK_DIR" || { log "Error: Cannot enter $WORK_DIR"; exit 1; }

# 3. Verify Input File exists
if [ ! -f "$INPUT_VCF" ]; then
    log "Error: Input file '$INPUT_VCF' not found in $WORK_DIR"
    exit 1
fi

# 4. Extract Sample Name (strips .vcf.gz, .vcf, or .bcf)
SAMPLE_NAME=$(basename "${INPUT_VCF}")
SAMPLE_NAME="${SAMPLE_NAME%.gz}"
SAMPLE_NAME="${SAMPLE_NAME%.vcf}"
SAMPLE_NAME="${SAMPLE_NAME%.bcf}"

# Create specific output folder for this file
mkdir -p "${OUTPUT_DIR}/${SAMPLE_NAME}"

# --- Processing ---
log "Starting variation summary for: ${SAMPLE_NAME}"

# Missingness Stats
log " -> Calculating missingness with bcftools & datamash"
# Note: Ensure +fill-tags was run previously to provide F_MISSING
summary_line=$(bcftools query -f '%F_MISSING\n' "$INPUT_VCF" | \
    $DATAMASH --round 3 mean 1 sstdev 1 median 1 min 1 max 1 count 1)

read mean sstdev median min max count <<< "$summary_line"

{
    echo "Summary for: $SAMPLE_NAME"
    echo "---------------------------"
    echo "Mean Missingness:   $mean"
    echo "Std. Deviation:     $sstdev"
    echo "Median:             $median"
    echo "Min:                $min"
    echo "Max:                $max"
    echo "Total Variants:     $count"
    echo ""
    echo "Raw Datamash: $summary_line"
} > "${OUTPUT_DIR}/${SAMPLE_NAME}/missing_summary.txt"

# BCFtools Stats
log " -> Generating full bcftools stats"
bcftools stats "$INPUT_VCF" > "${OUTPUT_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_stats.txt"

# Plotting (Visualizing Ts/Tv, Qualities, and Indel Distribution)
# 
if command -v plot-vcfstats &> /dev/null; then
    log " -> Creating plots in ${OUTPUT_DIR}/${SAMPLE_NAME}/plots"
    plot-vcfstats -p "${OUTPUT_DIR}/${SAMPLE_NAME}/plots" \
                  "${OUTPUT_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_stats.txt"
else
    log " -> plot-vcfstats not found in environment, skipping PDF generation"
fi

log "Pipeline completed successfully for ${SAMPLE_NAME}"
