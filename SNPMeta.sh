#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=8
#SBATCH --mem=200g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

set -e
set -o pipefail

module load python3/3.8.3_anaconda2020.07_mamba
module load emboss/6.6.0

# This script will run SNPMeta on existing BLAST XML reports

SNPMETA=/panfs/roc/groups/9/morrellp/pmorrell/shared/Software/SNPMeta/SNPMeta.py

WORKING_DIR=/scratch.global/pmorrell/Cowpea/test/

cd ${WORKING_DIR}

#for i in *.fasta
#    do

python3 ${SNPMETA} \
-d '/scratch.global/pmorrell/Cowpea/test' \
-a 'pmorrell@umn.edu' \
-l 60 \
--no-blast \
--outfmt tabular \
-o '/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNPMeta/Cowpea_SNPMeta.txt'

#    done
