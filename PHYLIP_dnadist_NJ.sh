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

   log "   -> Clean input and output directories"
rm -f "${INPUT_DIR}/outfile" "${OUTPUT_DIR}/outfile" "${INPUT_DIR}/outtree" "${OUTPUT_DIR}/outtree" 


process_phylip() {
    local PHYLIP_FILE=$1
    local BASENAME
    BASENAME=$(basename "${PHYLIP_FILE}" vcf.gz.min4.phy)
    log "Processing PHYLIP file: ${BASENAME}"

    # Step 1: Run dnadist
    log "   -> Running dnadist"
    expect <<EOF
spawn dnadist
expect "dnadist: can't find input file "infile"\n Please enter a new file name> "
send "${BASENAME}\r"
expect "Y to accept these or type the letter for one to change"
send "Y\r"
expect eof
EOF

mv outfile "${BASENAME}.dist"

    # Step 2: Run neighbor
    log "   -> Running neighbor"
    expect <<EOF
spawn neighbor
expect "neighbor: can't find input file "infile"\n Please enter a new file name> "
send "${BASENAME}.dist\r"
expect " Y to accept these or type the letter for one to change"
send "Y\r"
expect eof
EOF

mv outtree "${OUTPUT_DIR}/${BASENAME}.tree"
    
    log "Finished processing ${BASENAME}"
}

# Process all PHYLIP files in the input directory
log "Looking for PHYLIP files in ${INPUT_DIR}..."
find "${INPUT_DIR}" -name "*.phy" | while read -r PHYLIP_FILE; do
    process_phylip "${PHYLIP_FILE}"
done

log "All PHYLIP files processed."