#!/bin/env bash

#    Peter L. Morrell, 26 July 2016, St. Paul, MN
#    Updated 27 October 2016
#    Dependencies: SRA Toolkit
#    Requires bash version 4+

module load sratoolkit

set -euo pipefail

#PBS -l mem=1000mb,nodes=1:ppn=1,walltime=72:00:00
#PBS -m abe
#PBS -M pmorrell@umn.edu

OUTPUT=/panfs/roc/scratch/pmorrell/MBE_Barley

cd ${OUTPUT}

#mapfile -t SRA < <(find . -maxdepth 1 -name '*.sra' -type f)
declare -a SRA_ARRAY=(<(find . -maxdepth 1 -name '*.sra' -type f))


#    an array of run numbers you plan to download
#    archive=(ERR271705 ERR271706 ERR271707 ERR271708 ERR271709 ERR271710 ERR271711 ERR271712 ERR271713 ERR271714 ERR271715)

for i in "${SRA_ARRAY[@]}"
    do
    fastq-dump --split-files -I -F --gzip "${i}"
    done
