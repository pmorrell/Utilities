#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

set -e
set -o pipefail

module load bwa/0.7.17

# User provided input arguments
REF_DIR=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils
REF_FILENAME=Vunguiculata_IT97K-499-35_v1.0.fa

cd ${REF_DIR}

bwa index ${REF_DIR}/${REF_FILENAME}
