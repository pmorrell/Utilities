#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

# A script for processing BAM files to get exome coverage estimation
# Origninal version Peter L. Morrell, 12 July 2021, Falcon Heights, MN
# Updated by Chaochih Liu & Derell Scott, May 2022

set -euf
set -o pipefail

# Dependencies
# Mosdepth version 0.3.3
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/shared/Software/mosdepth

# Create list of BAM files, this is user provided input
# Use one of these approaches: `find /path/to/dir -type f > listOfFiles.list` or `ls -d "$PWD"/* > listOfFiles.list`
BAM_LIST="/panfs/roc/groups/9/morrellp/scot1101/Workshop/WBDC_125bp_MP005.txt"
# The size of the window in bp for coverage estimation.
INCREMENT=100
REF_FASTA="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"
OUT_DIR="/scratch.global/pmorrell/Morex_v3"

#------------------
mkdir -p ${OUT_DIR}
# Go into out directory
cd ${OUT_DIR}
# Prepare array
BAM_ARR=($(cat ${BAM_LIST}))

# Iterate through the BAM files and run Mosdepth on each
for i in ${BAM_ARR[@]}
do
    echo "Currently processing bam file: ${i}"
    # Strip the file to basename
    SAMPLE_NAME=$(basename ${i} .bam)
    mosdepth --by ${INCREMENT} --threads 4 --fast-mode --fasta ${REF_FASTA} --no-per-base ${SAMPLE_NAME} ${i}
done
