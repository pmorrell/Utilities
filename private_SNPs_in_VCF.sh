#!/bin/env bash

set -euo pipefail

module load vcftools_ML/0.1.14

#    A script to identify private SNPs within a sample
#    Peter L. Morrell, 05 May 2017, St. Paul, MN
#    Dependencies: VCFTools

#    Requires three input files:
#    1) the sample where private allele will be identified
#    2) a list of SNPs to be included [optional]
#    3) 

#while getopts ":h::s:" o; do
#	case $arg in 
#		s) # Specify sample list.
#           SAMPLE=${OPTARG}
#           ;;
#       h | *) # Display help.
#          usage
#            exit 0
#            ;;
#        esac
#    done

if [ $# -ne 3 ]; then
    echo $0: usage: "private_SNPs_in_VCF.sh [full VCF] [SNP list] [sample]"
    exit 1
fi

full_vcf=$1
SNP_list=$2
sample_name=$3

sample=$(basename $sample_name)

# /panfs/roc/scratch/tkono/GATK_Capture_WithID.vcf.gz
vcftools --gzvcf $full_vcf \
    --snps $SNP_list \
    --recode \
    --mac 1 \
    --remove $sample_name \
    --out ${sample}_removed

vcftools --gzvcf $full_vcf \
    --snps $SNP_list \
    --recode \
    --mac 1 \
    --keep $sample_name \
    --out ${sample}_only
 
vcftools --vcf ${sample}_removed.recode.vcf \
    --diff ${sample}_only.recode.vcf --diff-site --out ${sample}
#    The output is only the prefix. The rest of the name will include ".diff.sites_in_files"
#    The file contains a report on every siteretained in the analysis. 

#    The 4th column of this file reports the number of sites
#    that are unique to the 2nd VCF file. Here that will 
awk '$4=="2" {++count} END {print count}' ${sample}.diff.sites_in_files
