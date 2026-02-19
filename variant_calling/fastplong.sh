#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=400g
#SBATCH --tmp=400g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

# Peter L. Morrell - St. Paul, MN  -  17 June 2025

set -eou pipefail

# Load required modules
module load parallel/20210822

# Script configuration
DIRS_LIST=${1:-"/scratch.global/pmorrell/WBDC_resequencing/directory_list.txt"}
BASE_OUTPUT_DIR="/scratch.global/pmorrell/WBDC_resequencing/processed_output"
mkdir -p "${BASE_OUTPUT_DIR}"

FASTPLONG="/users/6/pmorrell/Apps/fastplong"
if [ ! -x "$FASTPLONG" ]; then
    echo "Error: fastplong is not found or not executable at $FASTPLONG"
    exit 1
fi

# Number of simultaneous jobs to run (adjust based on memory requirements)
PARALLEL_JOBS=4

# Validate resource allocation
if [ $((SLURM_CPUS_PER_TASK / PARALLEL_JOBS)) -lt 1 ]; then
    echo "Warning: Not enough CPUs for optimal parallelization"
    echo "Consider reducing PARALLEL_JOBS or increasing --cpus-per-task"
fi

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - [$$] ${msg}" | tee -a "${BASE_OUTPUT_DIR}/processing.log"
}

validate_fastq() {
    local FASTQ_FILE="$1"

    # Check if file exists and is readable
    if [ ! -r "$FASTQ_FILE" ]; then
        log "ERROR: Cannot read FASTQ file: $FASTQ_FILE"
        return 1
    fi

    # Basic validation - check if it's a valid gzipped file
    if ! gzip -t "$FASTQ_FILE" 2>/dev/null; then
        log "ERROR: Invalid gzip file: $FASTQ_FILE"
        return 1
    fi

    return 0
}

process_fastq() {
    local FASTQ_FILE="$1"
    local OUTPUT_DIR="$2"
    local SAMPLE_NAME=$(basename "${FASTQ_FILE}" .fastq.gz)
    local OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_processed.fastq.gz"

    # Skip if output already exists (resume capability)
    if [ -f "$OUTPUT_FILE" ]; then
        log "Skipping ${SAMPLE_NAME} - output already exists"
        return 0
    fi

    # Validate input file
    if ! validate_fastq "$FASTQ_FILE"; then
        return 1
    fi

    # Calculate threads per job
    local THREADS_PER_JOB=$((SLURM_CPUS_PER_TASK / PARALLEL_JOBS))
    [ $THREADS_PER_JOB -lt 1 ] && THREADS_PER_JOB=1

    log "Processing ${SAMPLE_NAME} with ${THREADS_PER_JOB} threads"

    local TEMP_OUTPUT="${OUTPUT_FILE}.tmp"

    if ${FASTPLONG} -i "${FASTQ_FILE}" \
                    -o "${TEMP_OUTPUT}" \
                    --threads ${THREADS_PER_JOB}; then

        mv "${TEMP_OUTPUT}" "${OUTPUT_FILE}"
        log "Successfully processed ${SAMPLE_NAME}"

        if ! gzip -t "${OUTPUT_FILE}" 2>/dev/null; then
            log "ERROR: Output file corrupted for ${SAMPLE_NAME}"
            rm -f "${OUTPUT_FILE}"
            return 1
        fi

        return 0
    else
        log "ERROR: Failed to process ${SAMPLE_NAME}"
        rm -f "${TEMP_OUTPUT}"
        return 1
    fi
}

export -f process_fastq
export -f log
export -f validate_fastq
export FASTPLONG
export BASE_OUTPUT_DIR
export SLURM_CPUS_PER_TASK
export PARALLEL_JOBS

process_directory() {
    local INPUT_DIR="$1"
    local DIR_NAME=$(basename "${INPUT_DIR}")
    local OUTPUT_DIR="${BASE_OUTPUT_DIR}/${DIR_NAME}"

    mkdir -p "${OUTPUT_DIR}"

    log "-> Processing directory: ${INPUT_DIR}"
    log "   -> Output will be saved to: ${OUTPUT_DIR}"

    local FASTQ_COUNT=$(find "${INPUT_DIR}" -name "*.fastq.gz" | wc -l)
    if [ $FASTQ_COUNT -eq 0 ]; then
        log "WARNING: No FASTQ files found in ${INPUT_DIR}"
        return 0
    fi

    log "   -> Found ${FASTQ_COUNT} FASTQ files to process"

    local FASTQ_LIST=$(mktemp)
    find "${INPUT_DIR}" -name "*.fastq.gz" | sort > "${FASTQ_LIST}"

    if parallel --jobs ${PARALLEL_JOBS} \
                --joblog "${OUTPUT_DIR}/parallel_joblog.txt" \
                --halt soon,fail=1 \
                --bar \
                process_fastq {} "${OUTPUT_DIR}" :::: "${FASTQ_LIST}"; then

        log "-> Successfully completed processing directory: ${INPUT_DIR}"

        local PROCESSED_COUNT=$(find "${OUTPUT_DIR}" -name "*_processed.fastq.gz" | wc -l)
        log "   -> Processed ${PROCESSED_COUNT}/${FASTQ_COUNT} files successfully"

    else
        log "ERROR: Failed to process some files in directory: ${INPUT_DIR}"
        rm -f "${FASTQ_LIST}"
        return 1
    fi

    rm -f "${FASTQ_LIST}"
    return 0
}

cleanup() {
    log "Cleaning up temporary files…"
    find /tmp -name "tmp.*" -user "$(whoami)" -mmin +60 -delete 2>/dev/null || true
}
trap cleanup EXIT

if [ ! -f "${DIRS_LIST}" ]; then
    echo "Error: Directory list file '${DIRS_LIST}' not found"
    echo "Usage: sbatch $0 [/path/to/directory_list.txt]"
    echo "If no argument provided, uses: /scratch.global/pmorrell/WBDC_resqeuencing/directory_list.txt"
    exit 1
fi

log "Starting fastplong processing for multiple directories"
log "Using directory list: ${DIRS_LIST}"
log "Output directory: ${BASE_OUTPUT_DIR}"
log "Parallel jobs: ${PARALLEL_JOBS}"
log "CPUs per task: ${SLURM_CPUS_PER_TASK}"

TOTAL_DIRS=0
PROCESSED_DIRS=0
FAILED_DIRS=0

while IFS= read -r dir; do
    [[ -z "$dir" || "$dir" =~ ^[[:space:]]*#.*$ ]] && continue

    dir=$(echo "$dir" | xargs)
    TOTAL_DIRS=$((TOTAL_DIRS + 1))

    if [ -d "$dir" ]; then
        if process_directory "$dir"; then
            PROCESSED_DIRS=$((PROCESSED_DIRS + 1))
        else
            FAILED_DIRS=$((FAILED_DIRS + 1))
            log "ERROR: Failed to process directory: $dir"
        fi
    else
        FAILED_DIRS=$((FAILED_DIRS + 1))
        log "ERROR: Directory does not exist: $dir"
    fi
done < "${DIRS_LIST}"

log "=== PROCESSING SUMMARY ==="
log "Total directories: ${TOTAL_DIRS}"
log "Successfully processed: ${PROCESSED_DIRS}"
log "Failed: ${FAILED_DIRS}"

if [ $FAILED_DIRS -gt 0 ]; then
    log "WARNING: Some directories failed to process. Check logs for details."
    exit 1
else
    log "All directories processed successfully!"
fi

