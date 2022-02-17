#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60gb
#SBATCH -t 02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -p small,large,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

module load bcftools/1.9

#    This script is intended to partition a VCF file into "core" and "noncore" variants
#    We will use a GFF file to find positions of core or noncore genes
#    We can

GFF3=/panfs/roc/groups/9/morrellp/shared/Datasets/Cowpea_Pan/VunguiculataIT97K-499-35_v1.2/annotation/Vunguiculata_IT97K-499-35_v1.2.gene.gff3.gz
VCF=/panfs/roc/groups/9/morrellp/shared/Datasets/Cowpea_Pan/VunguiculataIT97K-499-35_v1.2/SNPs/IT97K_combined_genotype_snps_filtered.g.vcf.gz
OUT_DIR=/scratch.global/pmorrell/Cowpea_pan
CORE=/panfs/roc/groups/9/morrellp/shared/Datasets/Cowpea_Pan/ITKcore.txt
NONCORE=/panfs/roc/groups/9/morrellp/shared/Datasets/Cowpea_Pan/ITKnoncore.txt

#    Cut the GFF3 file down to just gene positions to reduce our search space
#    The GFF3 file was compressed already
GENES=$(zgrep 'gene' ${GFF3})

#    Create a sorted bed file with only CORE or only NONCORE gene positions
grep -E -f "${CORE}" <(echo "${GENES}") | cut -f 1,4,5 | sort -k1,1 -k2,2n -k3,3n >${OUT_DIR}/CORE_POS.txt
grep -E -f "${NONCORE}" <(echo "${GENES}") | cut -f 1,4,5| sort -k1,1 -k2,2n -k3,3n >${OUT_DIR}/NONCORE_POS.txt

#    Use bcftools to cut the VCF down to individual SNPs within genes from each list
bcftools view --regions-file ${OUT_DIR}/CORE_POS.txt ${VCF} --output-type z --output-file ${OUT_DIR}/CORE_vcf.gz
bcftools view --regions-file ${OUT_DIR}/NONCORE_POS.txt ${VCF} --output-type z --output-file ${OUT_DIR}/NONCORE_vcf.gz
