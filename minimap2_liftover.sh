
#!/bin/bash -l

#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16g
#SBATCH --tmp=24g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.minimap2.out
#SBATCH -e %j.minimap2.err

# minimap2 liftover script (adds SLURM/module setup; default preset = sr)
# Usage: ./minimap2_liftover.sh <query.fasta> <reference.fasta> <output_prefix>

set -euo pipefail


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
# Default minimap2 settings: keep short-read preset (`sr`) by default,
# but allow override via MINIMAP2_PRESET env var (or export before running).
MINIMAP2_PRESET=${MINIMAP2_PRESET:-sr}
# Threads: prefer SLURM_CPUS_PER_TASK if present, otherwise THREADS env or 8
THREADS=${SLURM_CPUS_PER_TASK:-${THREADS:-8}}

echo "Mapping ${QUERY_FILE} to ${REF_FILE} using preset=${MINIMAP2_PRESET}, threads=${THREADS}..."
echo "Output: ${OUT_PREFIX}.sam"

# Run minimap2 and write SAM to output prefix
minimap2 -x "${MINIMAP2_PRESET}" -t "${THREADS}" -a "${REF_FILE}" "${QUERY_FILE}" > "${OUT_PREFIX}.sam"

echo "Done!"
