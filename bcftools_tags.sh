#!/bin/bash -l
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=48g
#SBATCH --tmp=48g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -u
set -o pipefail

# Check for required arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_bcf> <output_directory>"
    echo "  input_bcf: Path to indexed .bcf file"
    echo "  output_directory: Directory where output will be written"
    exit 1
fi

# Parse arguments
BCF="$1"
OUTDIR="$2"

# Validate input file exists
if [ ! -f "$BCF" ]; then
    echo "Error: Input BCF file does not exist: $BCF"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

BCFTOOLS=/users/6/pmorrell/software/modulesfiles/bcftools/bcftools

SAMPLE_NAME=$(basename "$BCF" .bcf)

"$BCFTOOLS" +fill-tags "$BCF" --threads 4 -Ob -o ${OUTDIR}/${SAMPLE_NAME}_tags.bcf --

bcftools index ${OUTDIR}/${SAMPLE_NAME}_tags.bcf
