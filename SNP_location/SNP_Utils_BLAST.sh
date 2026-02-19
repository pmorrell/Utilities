#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

#   Written by Peter L. Morrell
#   27 April 2021, Falcon Heights, MN

module load python3/3.8.3_anaconda2020.07_mamba
module load ncbi_blast+/2.8.1

set -euf - o pipefail

#   The application!
SNP_UTILS=/panfs/roc/groups/9/morrellp/shared/Software/SNP_Utils/snp_utils.py

#   Working directory
WORK=/panfs/roc/groups/9/morrellp/shared/Software/SNP_Utils/

#   Input file
SNPS=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/iSelect_all.txt

#   Config file
CONFIG=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/Vunguiculata_IT97K-499-35_v1.0.ini

#   The reference genome used to create a database
#REF=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/Vunguiculata_IT97K-499-35_v1.0.fa

# Output basenamne
OUT=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/Vunguiculata_IT97K-499-35_v1.0.vcf

#   Output should be written back to data directory
OUTPUT_DIR=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils

#Check if the directory exists if not make it
mkdir -p \
"$OUTPUT_DIR"

cd "$WORK"
"$SNP_UTILS" BLAST -l "$SNPS" -c "$CONFIG" -t 100000 -o "$OUT"
