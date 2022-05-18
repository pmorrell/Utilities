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

set -e
set -o pipefail

#    This script runs SamPlot according to https://github.com/ryanlayer/samplot

BED=/panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/Filtered/morex-sample2_diffs_from_ref.SVs_not_BND.bed
OUT_DIR=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Barley_Mutants
#    The gff3 converted to a bed file seems to work, gff3 did not
IMPORTANT=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/gene_annotation/Barley_Morex_V2_gene_annotation_PGSB.all.parts.sorted.bed
#    This bed file lists all intervals with 'Ns' in the reference genome
ANNOTATE=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Barley_Morex_V2_pseudomolecules_parts_missing.bed.gz

while IFS=$'\t' read -r -a array
do

#    Send our output here
cd ${OUT_DIR} || exit

#Samplot command line
samplot plot \
--sample_ids Morex_10X Morex_Nanopore \
--out-dir ${OUT_DIR}/ \
--output_type png \
--important_regions ${IMPORTANT} \
-A ${ANNOTATE} \
--bam /panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/longranger_morex_v2/morex-sample2/morex-sample2_phased_possorted_bam.bam \
/panfs/roc/groups/9/morrellp/shared/Projects/WBDC_inversions/nanopore/Morex/Morex_run1_to_run14/Morex_nanopore_1_14_V2_parts/Morex_1_14_align_V2_sorted_parts.bam \
-t DEL \
--chrom "${array[0]}" \
--start "${array[1]}" \
--end "${array[2]}" \
--output_file DEL_"${array[0]}"_"${array[1]}"_"${array[2]}".png

sleep 5s
done < ${BED}
