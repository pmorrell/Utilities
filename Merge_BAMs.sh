#!/bin/env bash

set -e
set -u
set -o pipefail

#PBS -l mem=1000mb,nodes=1:ppn=1,walltime=10:00:00
#PBS -m abe
#PBS -M pmorrell@umn.edu
#PBS -q lab

module load samtools/1.3

WORKING=${HOME}/scratch/MBE_Barley2/SAM_Processing/SAMtools/Finished

declare -a BAM=($(find ${WORKING} -maxdepth 1 -name '*.bam' -type f))

IFS=$'\n' BAM_sorted=($(sort <<<"${BAM[*]}"))

printf "%s\n" "${BAM_sorted[@]}"
#printf '%s\n' "${BAM[@]}"
#declare -a NAMES=($(basename ))
