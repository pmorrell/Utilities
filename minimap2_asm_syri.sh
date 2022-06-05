#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=3
#SBATCH --mem=400g
#SBATCH --tmp=400g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

# This script will run minimap2, which has a wide variety of applications for long reads.
# The mapping parameters below are chosen specifically for running syri for detection of major structural variants.

MINIMAP2=/panfs/roc/groups/9/morrellp/pmorrell/Apps/HLi/minimap2/minimap2

#REF_ARRAY=(HvulgareMorex_702_V3.hardmasked_chr1H.fa.gz HvulgareMorex_702_V3.hardmasked_chr2H.fa.gz HvulgareMorex_702_V3.hardmasked_chr3H.fa.gz HvulgareMorex_702_V3.hardmasked_chr4H.fa.gz HvulgareMorex_702_V3.hardmasked_chr5H.fa.gz HvulgareMorex_702_V3.hardmasked_chr6H.fa.gz HvulgareMorex_702_V3.hardmasked_chr7H.fa.gz) 
REF_ARRAY=(HvulgareMorex_702_V3.hardmasked_chr7H.fa.gz)
#QUERY_ARRAY=(OUH602_chr1H.fasta.gz OUH602_chr2H.fasta.gz OUH602_chr3H.fasta.gz OUH602_chr4H.fasta.gz OUH602_chr5H.fasta.gz OUH602_chr6H.fasta.gz OUH602_chr7H.fasta.gz)
QUERY_ARRAY=(OUH602_chr7H.fasta.gz)
#NAME_ARRAY=(chr1H chr2H chr3H chr4H chr5H chr6H chr7H)
NAME_ARRAY=(chr7H)

for i in "${!REF_ARRAY[@]}"
 do

# User provided input arguments
REF_FILE=/scratch.global/pmorrell/Morex_v3/hardmasked/${REF_ARRAY[i]}
FASTA_FILE=/scratch.global/pmorrell/OUH602/${QUERY_ARRAY[i]}
OUT_PREFIX=OUH602_${NAME_ARRAY[i]}_asm5
OUT_DIR=/scratch.global/pmorrell/Morex_v3/out_dir

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into reference dir
cd ${OUT_DIR}

# Align with minimap
# The '--eqx' option it chosen specifically to alter the CIGAR string for syri analysis
$MINIMAP2 -c -x asm5 --cs -t 3 -r2k --eqx ${REF_FILE} ${FASTA_FILE} > ${OUT_DIR}/${OUT_PREFIX}.paf
done

