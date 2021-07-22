#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

#   Written by Peter L. Morrell
#   12 July 2021, Falcon Heights, MN

set -euf
set -o pipefail

MOSDEPTH=/panfs/roc/groups/9/morrellp/shared/Software/mosdepth/mosdepth

OUT_PREFIX=WBDC355_Nanopore

WORK=/scratch.global/pmorrell/WBDC355_Nanopore

INCREMENT=500

BAM=/panfs/roc/groups/9/morrellp/shared/Projects/WBDC_inversions/nanopore/WBDC_355/WBDC_combined_1_10/WBDC_355_v3/bam_file/WBDC_355_1_10_v3_sorted.bam

#SNPS=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/iSelect_all.fas

#GENOME=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/Vunguiculata_IT97K-499-35_v1.0.fa

#OUT=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Cowpea/SNP_Utils/cowpea_snps.xml


# query SNPS and write to XML output
cd "$WORK"
$MOSDEPTH --no-per-base --threads 4 --fast-mode --by $INCREMENT $OUT_PREFIX $BAM
