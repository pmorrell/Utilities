#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --ntasks=5
#SBATCH --mem=24g
#SBATCH --tmp=24g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

# Peter L. Morrell - St. Paul, MN  -  18 March 2025

set -e
set -o pipefail

# Load required modules
module load phylip/3.697

INPUT_DIR="/scratch.global/pmorrell/Barley_Introgression/PHYLIP/"  # Directory containing PHYLIP files
OUTPUT_DIR="/scratch.global/pmorrell/Barley_Introgression/trees"       # Directory to store output files

mkdir -p "${OUTPUT_DIR}"

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

process_phylip() {
    local PHYLIP_FILE=$1
    local BASENAME
    BASENAME=$(basename "${PHYLIP_FILE}" .phy)
    log "Processing PHYLIP file: ${BASENAME}"

    # Step 1: Run dnadist
    log "   -> Running dnadist"
    expect <<EOF
spawn dnadist
expect "Input file:"
send "${PHYLIP_FILE}\r"
expect "Output file:"
send "${OUTPUT_DIR}/${BASENAME}_dnadist.out\r"
expect "Substitution model"
send "F\r"  # Example: F for F84 model (modify as needed)
expect "Transition/transversion ratio"
send "2.01\r"  # Example: Set ratio to 2.01 (modify as needed)
expect "Analyze another data set?"
send "N\r"
expect eof
EOF

    # Step 2: Run neighbor
    log "   -> Running neighbor"
    expect <<EOF
spawn neighbor
expect "Input file:"
send "${OUTPUT_DIR}/${BASENAME}_dnadist.out\r"
expect "Output file:"
send "${OUTPUT_DIR}/${BASENAME}_neighbor.out\r"
expect "Tree file:"
send "${OUTPUT_DIR}/${BASENAME}_treefile\r"
expect "Analyze another data set?"
send "N\r"
expect eof
EOF

    log "Finished processing ${BASENAME}"
}

# Process all PHYLIP files in the input directory
log "Looking for PHYLIP files in ${INPUT_DIR}..."
find "${INPUT_DIR}" -name "*.phy" | while read -r PHYLIP_FILE; do
    process_phylip "${PHYLIP_FILE}"
done

log "All PHYLIP files processed."