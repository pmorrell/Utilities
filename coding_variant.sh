#!//usr/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

# A script to identify the proportion of variants that are in coding sequence.

#   Written by Peter L. Morrell
#   19 November 2021, Falcon Heights, MN

set -euf
set -o pipefail

module load bedtools/2.29.2

OUTDIR=/panfs/roc/groups/9/morrellp/shared/Projects/Selective_Sweeps
OUTFILE=code_variants.txt

bedtools intersect -wa -wb \
-a /panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation/Barley_Morex_V2_gene_annotation_PGSB.all.parts.bed \
-b /panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/selective_sweeps_morex_v2/all/Variant_Recalibrator/ssw_vc_snps_final.vcf -names d1 -sorted \
| uniq > ${OUTDIR}/${OUTFILE}
