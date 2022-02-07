#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gfrascar@umn.edu
#SBATCH -p ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

#This script runs SamPlot according to https://github.com/ryanlayer/samplot

#to use samplot
source activate /home/morrellp/public/Software/anaconda3

cd /panfs/roc/groups/9/morrellp/pmorrell/Workshop/
VCF=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/M01_temp.vcf
OUT_DIR=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Barley_Mutants
GFF3=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation/Barley_Morex_V2_gene_annotation_PGSB.all.sorted.gff3

#Samplot command line
samplot vcf \
--vcf ${VCF} \
--out-dir ${OUT_DIR} \
--output_type png \
--gff3 ${GFF3} \
--bam /panfs/roc/groups/9/morrellp/shared/Projects/WBDC_inversions/nanopore/Morex/Morex_run1_to_run14/Morex_nanopore_1_14_V2_parts/Morex_1_14_align_V2_sorted_parts.bam \
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/morex-sample2_phased_possorted_bam.bam \
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/M01-3-3/M01-3-3_phased_possorted_bam.bam \
/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/sam_processing/M01-3-3-12-41_run1-3.bam
>samplot_commands.sh
