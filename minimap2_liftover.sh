
#!/bin/bash

# minimap2 script for mapping short genomic sequences with SNP
# Usage: ./minimap2_mapping.sh <query.fasta> <reference.fasta> <output_prefix>

set -euo pipefail

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <query.fasta> <reference.fasta> <output_prefix>"
    echo ""
    echo "Arguments:"
    echo "  query.fasta       - Input sequence file (121 bp with SNP)"
    echo "  reference.fasta   - Reference genome"
    echo "  output_prefix     - Prefix for output files"
    exit 1
fi

QUERY_FILE="$1"
REF_FILE="$2"
OUT_PREFIX="$3"
OUT_DIR=$(dirname "$OUT_PREFIX")

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# minimap2 parameters for short sequences:
# -x sr = short reads (good for ~100bp sequences)
# -a = output SAM format
# -o = output file

echo "Mapping ${QUERY_FILE} to ${REF_FILE}..."
echo "Output: ${OUT_PREFIX}.sam"

# Run minimap2
minimap2 -x sr -a -o "${OUT_PREFIX}.sam" "$REF_FILE" "$QUERY_FILE"

echo "Done!"
