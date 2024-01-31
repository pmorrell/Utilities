#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --ntasks=8
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH --mail-type=ALL
#SBATCH -p small
#SBATCH --mail-user=pmorrell@umn.edu

#   Peter L. Morrell - 30 Jan 2024 - Falcon Heights, MN
#   Use bioawk and the example at: https://github.com/vsbuffalo/bioawk-tutorial
#   Calculate number of reads per FASTQ file. Primary issue is determining if
#   paired end reads are present and complete before concatenation of fastq.

set -e
set -u
set -o pipefail

module load parallel

WORKING_DIR=/scratch.global/pmorrell/OUH602

cd $WORKING_DIR
FASTQ=(*.fastq.gz)

for i in "${FASTQ[@]}"
    do
        # echo "$i" >> FASTQ_lenghts.txt
        parallel "echo {} && gunzip -c {} | wc -l | awk '{d=\$1; print d/4;}'"\
         ::: "$i" >> $WORKING_DIR/FASTQ_lenghts.txt
        #$BIOAWK -cfastx '{print $name, length($seq)}' "$i" >> $WORKING_DIR/FASTQ_lenghts.txt
    done

