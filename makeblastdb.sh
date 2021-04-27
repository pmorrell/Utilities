#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

#   Written by Peter L. Morrell
#   27 April 2021, Falcon Heights, MN

module load ncbi_blast+/2.8.1

set -euf - o pipefail

#   The reference genome used to create a database
REF=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/Vunguiculata_IT97K-499-35_v1.0.faß

# Output basenamne
OUT=Vunguiculata_IT97K-499-35_v1.0.fa

#   Output should be written back to data directory
OUTPUT_DIR=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils

#Check if the directory exists if not make it
mkdir -p \
"$OUTPUT_DIR"

cd "$OUTPUT_DIR"
makeblastdb -in "$REF" -out "$OUT" -dtype nucl -parse_seqids

