#!/bin/bash -l
#
# ONT genomic reads alignment to barley reference using minimap2.
# Mirrors minimap2_PacBio.sh but uses ONT preset and outputs sorted/indexed BAM.
#
# SLURM settings (adjust as needed)
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=200g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -euo pipefail

# Optionally load modules if using an environment module system
# module load minimap2/2.26
module load samtools/1.21

# Paths to tools (override if needed)
MINIMAP2=/users/6/pmorrell/Apps/HLi/minimap2/minimap2

# User-provided inputs (edit these for each run)
REF_FILE=/home/morrellp/pmorrell/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_plastids.fasta
# Optional pre-built minimap2 index (.mmi). If set and exists, it will be used instead of FASTA
REF_IDX=/home/morrellp/pmorrell/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_plastids.mmi

SAMPLE=WBDC355
OUT_PREFIX=${SAMPLE}_ONT
# IN_DIR=/scratch.global/pmorrell/WBDC_resequencing/${SAMPLE}
OUT_DIR=/scratch.global/pmorrell/WBDC_resequencing/${SAMPLE}/ONT
# Path to a file containing a list of FASTQ files (one per line)
# cd ${IN_DIR}
# find "$PWD" -type f -name "*.fastq.gz" > ${SAMPLE}_fastq.txt
FASTQ_LIST=/scratch.global/pmorrell/WBDC_resequencing/WBDC355_ONT/fastplong/WBDC_355_ONT_filtered_fastq.txt

# Threads used by minimap2 and samtools sort
THREADS=8
# Memory per minimap2 index split (tune for large references)
INDEX_SPLIT=-I6g

# Ensure output directory exists
mkdir -p "${OUT_DIR}"

# Move to output directory to keep outputs together
cd "${OUT_DIR}"

# Decide which reference to use (MMI index if provided and exists, else FASTA)
if [[ -n "${REF_IDX:-}" && -f "${REF_IDX}" ]]; then
  REF_INPUT="${REF_IDX}"
else
  REF_INPUT="${REF_FILE}"
fi

# Log tool versions
"${MINIMAP2}" --version >&2 || true
samtools --version 2>/dev/null | head -n1 >&2 || true

# Read group (optional but recommended for downstream tools)
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ONT"

# Align ONT reads and output sorted, indexed BAM
# Also write the exact command to a .cmd.txt for provenance
{
  echo "Running minimap2 ONT alignment at $(date)"
  echo "Reference input: ${REF_INPUT}"
  echo "FASTQ list: ${FASTQ_LIST}"
  echo "Output prefix: ${OUT_PREFIX}"
  echo "Threads: ${THREADS}"
} >&2

ALIGN_CMD=("${MINIMAP2}" -ax map-ont -t "${THREADS}" ${INDEX_SPLIT} -R "${RG}" "$REF_INPUT")
FASTQS=( $(cat "${FASTQ_LIST}") )

# Save the exact command used
printf '%q ' "${ALIGN_CMD[@]}" "${FASTQS[@]}" '|' "${SAMTOOLS}" sort -@ "${THREADS}" -o "${OUT_DIR}/${OUT_PREFIX}.sorted.bam" > "${OUT_PREFIX}.cmd.txt"

# Execute
"${ALIGN_CMD[@]}" "${FASTQS[@]}" \
  | "${SAMTOOLS}" sort -@ "${THREADS}" -o "${OUT_DIR}/${OUT_PREFIX}.sorted.bam"

# Index BAM
"${SAMTOOLS}" index -@ "${THREADS}" "${OUT_DIR}/${OUT_PREFIX}.sorted.bam"

# Optional: generate alignment stats
"${SAMTOOLS}" flagstat -@ "${THREADS}" "${OUT_DIR}/${OUT_PREFIX}.sorted.bam" > "${OUT_DIR}/${OUT_PREFIX}.flagstat.txt"

# Done
>&2 echo "Done. Outputs:"
>&2 echo "  BAM:     ${OUT_DIR}/${OUT_PREFIX}.sorted.bam"
>&2 echo "  BAI:     ${OUT_DIR}/${OUT_PREFIX}.sorted.bam.bai"
>&2 echo "  CMD:     ${OUT_DIR}/${OUT_PREFIX}.cmd.txt"
>&2 echo "  FLAGSTAT:${OUT_DIR}/${OUT_PREFIX}.flagstat.txt"

