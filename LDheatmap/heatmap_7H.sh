#!/bin/bash

#PBS -l mem=2000mb,nodes=1:ppn=4,walltime=8:00:00
#PBS -m abe
#PBS -M pmorrell@umn.edu
#PBS -q lab     

cd ~/R

module load R

R CMD BATCH --no-save CAP_heatmap_7H.r