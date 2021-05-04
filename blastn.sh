#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

#   Written by Peter L. Morrell
#   03 May 2021, Falcon Heights, MN

set -euf - o pipefail

module load ncbi_blast+/2.8.1

WORK=/panfs/roc/groups/9/morrellp/shared/Software/SNP_Utils/

SNPS=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/iSelect_all.fas

GENOME=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/Vunguiculata_IT97K-499-35_v1.0.fa

OUT=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/cowpea_snps.xml


# query SNPS and write to XML output
cd "$WORK"
blastn -db "$GENOME" -query "$SNPS" -outfmt 5 -out "$OUT"
