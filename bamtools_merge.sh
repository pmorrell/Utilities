#!/bin/env bash

#PBS -l mem=28gb,nodes=1:ppn=7,walltime=48:00:00
#PBS -m abe
#PBS -M pmorrell@umn.edu
#PBS -q mesabi

set -e
set -u
set -o pipefail

module load bamtools

SAMPLE_LIST="/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Merge/Barke_list.txt"
SAMPLE="Barke_finished.bam"

bamtools -list ${SAMPLE_LIST} -out ${SAMPLE}