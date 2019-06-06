#!/usr/bin/env bash

#PBS -l mem=6000mb,nodes=1:ppn=8,walltime=3:00:00
#PBS -m abe
#PBS -M pmorrell@umn.edu
#PBS -q lab

#    Load necessary modules
module load minimap2_ML/2.17.0
#    Genome sequence to be indexed
GENOME=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_plastids.fasta.gz
#    Define output file
MMI=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules_plastids.mmi

mimimap2 -d $MMI $GENOME
