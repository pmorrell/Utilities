#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -p ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

#This script runs SamPlot according to https://github.com/ryanlayer/samplot

#Not sure if I need a python version with biopyton or if it comes with the conda environment
#module load python3/3.7.4_anaconda2019.10
#conda activate samplot

VCF=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/combined_mutated/Filtered/deletions/M01_singletons_only-scored_DEL_gte75Sup.vcf
OUT_DIR=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Barley_Mutants
IMPORTANT=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation/Barley_Morex_V2_gene_annotation_PGSB.all.parts.sorted.bed
ANNOTATE=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Barley_Morex_V2_pseudomolecules_parts_missing.bed.gz

#    Send our output here
cd ${OUT_DIR} || exit

#Samplot command line
samplot vcf \
--vcf ${VCF} \
--sample_ids Morex_10X Morex_Nanopore M01_Nanopore M01_10X \
--out-dir ${OUT_DIR}/ \
--output_type png \
--important_regions ${IMPORTANT} \
-A ${ANNOTATE} \
--bams /panfs/roc/groups/9/morrellp/shared/Projects/WBDC_inversions/nanopore/Morex/Morex_run1_to_run14/Morex_nanopore_1_14_V2_parts/Morex_1_14_align_V2_sorted_parts.bam \
/panfs/roc/groups/9/morrellp/shared/Projects/WBDC_inversions/nanopore/Morex/Morex_run1_to_run14/Morex_nanopore_1_14_V2_parts/Morex_1_14_align_V2_sorted_parts.bam \
/panfs/roc/groups/9/morrellp/shared/Datasets/Alignments/nanopore_mutated_barley/sam_processing/M01-3-3-12-41_run1-3.bam \
/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/M01-3-3/M01-3-3_phased_possorted_bam.bam \
> samplot_commands.sh
