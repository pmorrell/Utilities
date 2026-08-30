#!/bin/bash -l
#SBATCH --time=64:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=48g
#SBATCH --tmp=48g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -euo pipefail

BCFTOOLS=/users/6/pmorrell/software/modulesfiles/bcftools/bcftools

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <VCF> <OUTDIR>" >&2
    exit 1
fi

# Path to bgzipped and indexed .vcf file, and output directory
VCF="$1"
OUTDIR="$2"

SAMPLE_NAME=$(basename "$VCF" .vcf.gz)

"$BCFTOOLS" +fill-tags "$VCF" --threads 4 -Oz -o ${OUTDIR}/${SAMPLE_NAME}.vcf.gz --
