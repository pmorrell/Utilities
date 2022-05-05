#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=8
#SBATCH --mem=200g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

set -e
set -o pipefail

# Load the tools needed
module load samtools/1.9
module load htslib/1.9
# The application and associated tools are all in the directory below
SMOOVE=/panfs/roc/groups/9/morrellp/shared/Software/smoove
OUT_DIR=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Soy_120/smoove/out_dir

# BAM file, reference, and
BAM=/panfs/roc/groups/9/morrellp/shared/Datasets/10x_Genomics/Soybean/m92_220/M92_220_phased_possorted_bam.bam
reference_fasta=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Soybean/PhytozomeV11/Gmax/assembly/Gmax_275_v2.0.fa
threads=2

# Check if our dir exists, if not make it
mkdir -p ${OUT_DIR}

# Go into output directory
cd ${OUT_DIR}

# A minimum set of positions for calling variants. Not including an "excludes" file with problematic regions
$SMOOVE call -x --name M92_220 --fasta  $reference_fasta -p $threads --genotype $BAM
