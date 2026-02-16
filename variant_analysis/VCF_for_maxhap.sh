#!/bin/bash

# Create input files for maxhap from a VCF
# Need to remove multiallelic SNPs, those with very low allele frequency, etc.

MIN_AC=3
CHR=chr1H
WILD=/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/vcf_morex_v3/wild_snps_final.vcf.psuedomolecule.vcf.gz
WORK_DIR=/scratch.global/Hudson_2site

mkdir -p $WORK_DIR

bcftools view -m2 -M2 -v snps --min-ac $MIN_AC --regions $CHR $WILD >$WORK_DIR/wild_${CHR}.vcf

sed -i.bak -e 's|1\|1|1\/1|g'  wild_chr1H_1_50000.vcf
sed -i.bak -e 's|0\|0|0\/0|g'  wild_chr1H_1_50000.vcf
sed -i.bak -e 's|0\|1|0\/1|g'  wild_chr1H_1_50000.vcf




sed 's/chr1H_//g' wild_chr1H_1_50000.txt >wild_1_50000.txt 

sed 's/AA/A/g' - |  sed 's/CC/C/g' - |  sed 's/GG/G/g' - | sed 's/TT/T/g' - | sed 's/NN/N/g' - > wild_chr1H_1_50000_hap.txt 
