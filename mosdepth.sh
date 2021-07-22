#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

# A script for genome-wide estimation of sequencing depth from a single sample.

#   Written by Peter L. Morrell
#   12 July 2021, Falcon Heights, MN

set -euf
set -o pipefail

MOSDEPTH=/panfs/roc/groups/9/morrellp/shared/Software/mosdepth/mosdepth

# Output files are required to have a prefix.
OUT_PREFIX=WBDC355_Nanopore

# Output directory
OUT_DIR=/scratch.global/pmorrell/WBDC355_Nanopore

# The size of the window in bp for coverage estimation.
INCREMENT=500

BAM=/panfs/roc/groups/9/morrellp/shared/Projects/WBDC_inversions/nanopore/WBDC_355/WBDC_combined_1_10/WBDC_355_v3/bam_file/WBDC_355_1_10_v3_sorted.bam

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into the 
cd ${OUT_DIR}

# The threads option is for compression. Don't increase beyond 4, it doesn't speedup the program!
$MOSDEPTH --no-per-base --threads 4 --fast-mode --by $INCREMENT $OUT_PREFIX $BAM
