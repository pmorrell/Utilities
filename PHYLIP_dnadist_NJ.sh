#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --ntasks=5
#SBATCH --mem=24g
#SBATCH --tmp=24g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

# Peter L. Morrell - St. Paul, MN - 19 March 2025

set -e
set -o pipefail

# Dependency checks
command -v expect >/dev/null 2>&1 || { echo "Error: 'expect' is required but not installed."; exit 1; }

# Load required modules
module load phylip/3.69

INPUT_DIR="/scratch.global/pmorrell/Barley_Introgression/PHYLIP"
OUTPUT_DIR="/scratch.global/pmorrell/Barley_Introgression/trees"
TEMP_DIR="${SLURM_TMPDIR:-/tmp}/phylip_${SLURM_JOB_ID:-$$}"

# Create directories
mkdir -p "${OUTPUT_DIR}" "${TEMP_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

process_phylip() {
    local PHYLIP_FILE=$1
    local BASENAME=$(basename "${PHYLIP_FILE}" .vcf.gz.min4.phy)
    local WORK_DIR="${TEMP_DIR}/${BASENAME}"
    
    mkdir -p "${WORK_DIR}"
    
    cd "${WORK_DIR}" || { log "Error: Cannot change to ${WORK_DIR}"; return 1; }
    
    log "Processing PHYLIP file: ${BASENAME}"
    cp "${PHYLIP_FILE}" infile
    rm outfile

    # Step 1: Run dnadist
    log "   -> Running dnadist"
    expect <<EOF

spawn dnadist
expect "Y to accept these or type the letter for one to change"
send "Y\r"
expect eof
EOF

    # Check if outfile was created
    if [[ ! -f "outfile" ]]; then
        log "Error: dnadist failed to create output file for ${BASENAME}"
        return 1
    fi

cp outfile infile
rm outfile
rm outtree

    # Step 2: Run neighbor
    log "   -> Running neighbor"
    expect <<EOF
spawn neighbor 
expect "Y to accept these or type the letter for one to change"
send "Y\r"
expect eof
EOF

    # Check if tree was created
    if [[ ! -f "outtree" ]]; then
        log "Error: neighbor failed to create tree file for ${BASENAME}"
        return 1
    fi

    # Copy final tree to output directory
    cp "outtree" "${OUTPUT_DIR}/${BASENAME}.tree"
    
    # Cleanup work directory
    cd "${TEMP_DIR}" || true
    
    log "Finished processing ${BASENAME}"
}

# Process all PHYLIP files in the input directory
log "Looking for PHYLIP files in ${INPUT_DIR}..."
find "${INPUT_DIR}" -name "*.vcf.gz.min4.phy" | while read -r PHYLIP_FILE; do
    process_phylip "${PHYLIP_FILE}"
done

# Cleanup
rm -rf "${TEMP_DIR}"

log "All PHYLIP files processed."