#!/bin/env bash

#PBS -l mem=1000mb,nodes=1:ppn=16,walltime=72:00:00
#PBS -m abe
#PBS -M pmorrell@umn.edu
#PBS -q mesabi

set -euo pipefail

#    A script to take a list of SRA files and convert them into forward and reverse reads in fastq format.
#    Files will also be gzipped.
#    This can take a while, so we run on MSI.
#    Peter L. Morrell, 26 July 2016, St. Paul, MN
#    Updated 27 October 2016
#    Dependencies: SRA Toolkit
#    Requires bash version 4+ for mapfile (not currently used)



module load sratoolkit
module load parallel

#    directory for out of fastq.gz files
WORKING=/panfs/roc/scratch/pmorrell/testing123

#   initalize the array that will hold a list of SRA files
#mapfile -t SRA < <(find ${WORKING} -maxdepth 1 -name '*.sra' -type f)
declare -a SRA=($(find ${WORKING} -maxdepth 1 -name '*.sra' -type f))
echo "${#SRA[@]} samples"

parallel --verbose "fastq-dump --split-files -I -F --gzip {} --outdir ${WORKING}" ::: ${SRA[@]}

