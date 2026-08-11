#!/usr/bin/env bash
#SBATCH --partition=msibigmem
#SBATCH --job-name=anchorwave_barley
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=15900M
#SBATCH --tmp=300g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH --array=1-7
#SBATCH --output=anchorwave_%A_%a.out
#SBATCH --error=anchorwave_%A_%a.err

set -euo pipefail

module load samtools/1.21
module load minimap2/2.30  # Adjust version to what is available on your cluster

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(timestamp)] $*"; }
show_log_tail() {
  local log_file="$1"
  if [[ -s "${log_file}" ]]; then
    log "Last 40 lines from ${log_file}:"
    tail -n 40 "${log_file}" >&2
  fi
}
check_sam() {
  local sam_file="$1"
  # samtools quickcheck does not reliably detect truncation in SAM files;
  # view -c reads through the file and will exit non-zero on truncation.
  samtools view -c "${sam_file}" >/dev/null 2>&1
}
check_bam() {
  local bam_file="$1"
  samtools quickcheck "${bam_file}" >/dev/null 2>&1
}

if [[ -x "${HOME}/Apps/anchorwave/anchorwave" ]]; then
  ANCHORWAVE_BIN="${HOME}/Apps/anchorwave/anchorwave"
else
  echo "ERROR: Could not find executable AnchorWave at ${HOME}/Apps/anchorwave/anchorwave" >&2
  exit 1
fi

readonly THREADS="${SLURM_CPUS_PER_TASK:-1}"
if (( THREADS > 2 )); then
  DEFAULT_RETRY_THREADS=2
else
  DEFAULT_RETRY_THREADS="${THREADS}"
fi
readonly RETRY_THREADS="${RETRY_THREADS:-${DEFAULT_RETRY_THREADS}}"
readonly CHRS=(chr1H chr2H chr3H chr4H chr5H chr6H chr7H)

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: This script is intended to run as a SLURM array job." >&2
  exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > ${#CHRS[@]} )); then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} is out of range." >&2
  exit 1
fi

chr="${CHRS[$((SLURM_ARRAY_TASK_ID - 1))]}"

# Inputs - Reference (Keep these exactly as they were)
REF_FA="/projects/standard/morrellp/pmorrell/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"
REF_GFF="/projects/standard/morrellp/pmorrell/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.sorted.gff3"
CDS_FA="/scratch.global/pmorrell/Morex_v3/cds.fa"
REF_SAM="/scratch.global/pmorrell/Morex_v3/ref.sam" # CDS mapped to Reference

# Inputs - Query (Update this to your unzipped Supernova output)
QRY_FA="/projects/standard/morrellp/pmorrell/shared/Datasets/10x_Genomics/Barley/WBDC355_Assembly_2018-04-28/barley_outputs_barley.pseudohap.fasta"

# Workspace
WORKDIR="per_chr_work"
OUTDIR="/scratch.global/pmorrell/WBDC355"
LOGDIR="per_chr_logs"
mkdir -p "${WORKDIR}" "${OUTDIR}"
mkdir -p "${LOGDIR}"
PREP_LOG="${LOGDIR}/preprocessing.log"

# Shared intermediate files
SCAFF_TO_REF_SAM="${WORKDIR}/scaffolds_to_ref.sam"
QRY_BAM="${WORKDIR}/qry_scaffolds.sorted.bam"
REF_BAM="${WORKDIR}/ref_cds.sorted.bam"
QRY_FAI="${QRY_FA}.fai"

# Build the shared BAMs once, even if array tasks start together.
exec 9>"${WORKDIR}/preprocessing.lock"
flock 9

# REF_SAM is treated as a trusted upstream input. This script only rebuilds
# files it creates itself, so a truncated-file recovery path here would be
# misleading when REF_SAM is managed elsewhere.
if [[ ! -s "${REF_SAM}" ]]; then
  log "ERROR: Required REF_SAM is missing: ${REF_SAM}"
  exit 1
fi

if [[ -s "${SCAFF_TO_REF_SAM}" ]] && ! check_sam "${SCAFF_TO_REF_SAM}"; then
  log "Existing scaffold SAM is invalid or truncated; rebuilding ${SCAFF_TO_REF_SAM}"
  rm -f "${SCAFF_TO_REF_SAM}"
fi

if [[ ! -s "${SCAFF_TO_REF_SAM}" ]]; then
  log "Mapping Supernova scaffolds to Reference with minimap2..." | tee -a "${PREP_LOG}"
  minimap2 -x asm20 -t "${THREADS}" -a -r 10000 -p 0.4 -N 20 "${REF_FA}" "${QRY_FA}" > "${SCAFF_TO_REF_SAM}.tmp"
  if ! samtools view -c "${SCAFF_TO_REF_SAM}.tmp" >/dev/null 2>&1; then
    log "ERROR: minimap2 output failed SAM validation; removing ${SCAFF_TO_REF_SAM}.tmp" | tee -a "${PREP_LOG}"
    rm -f "${SCAFF_TO_REF_SAM}.tmp"
    exit 1
  fi
  mv "${SCAFF_TO_REF_SAM}.tmp" "${SCAFF_TO_REF_SAM}"
fi

if [[ -s "${REF_BAM}" ]] && ! check_bam "${REF_BAM}"; then
  log "Existing REF_BAM is invalid or truncated; rebuilding ${REF_BAM}"
  rm -f "${REF_BAM}" "${REF_BAM}.csi"
fi

if [[ ! -s "${REF_BAM}" || ! -s "${REF_BAM}.csi" ]]; then
  log "Sorting and indexing Reference BAM from ${REF_SAM}..." | tee -a "${PREP_LOG}"
  set +e
  samtools sort -@ "${THREADS}" -o "${REF_BAM}.tmp" "${REF_SAM}" >> "${PREP_LOG}" 2>&1
  ref_sort_rc=$?
  set -e
  if [[ "${ref_sort_rc}" -ne 0 ]]; then
    log "Reference BAM sort failed with exit code ${ref_sort_rc}" | tee -a "${PREP_LOG}"
    show_log_tail "${PREP_LOG}"
    exit "${ref_sort_rc}"
  fi
  mv "${REF_BAM}.tmp" "${REF_BAM}"
  samtools index -@ "${THREADS}" --csi "${REF_BAM}"
fi

if [[ -s "${QRY_BAM}" ]] && ! check_bam "${QRY_BAM}"; then
  log "Existing QRY_BAM is invalid or truncated; rebuilding ${QRY_BAM}"
  rm -f "${QRY_BAM}" "${QRY_BAM}.csi"
fi

if [[ ! -s "${QRY_BAM}" || ! -s "${QRY_BAM}.csi" ]]; then
  log "Sorting and indexing Query Scaffolds BAM from ${SCAFF_TO_REF_SAM}..." | tee -a "${PREP_LOG}"
  set +e
  samtools sort -@ "${THREADS}" -o "${QRY_BAM}.tmp" "${SCAFF_TO_REF_SAM}" >> "${PREP_LOG}" 2>&1
  qry_sort_rc=$?
  set -e
  if [[ "${qry_sort_rc}" -ne 0 ]]; then
    log "Query scaffolds BAM sort failed with exit code ${qry_sort_rc}" | tee -a "${PREP_LOG}"
    show_log_tail "${PREP_LOG}"
    exit "${qry_sort_rc}"
  fi
  mv "${QRY_BAM}.tmp" "${QRY_BAM}"
  samtools index -@ "${THREADS}" --csi "${QRY_BAM}"
fi

# Index the query assembly once, shared across all chromosome tasks, so we
# can pull individual scaffolds out of it below with `samtools faidx`.
if [[ ! -s "${QRY_FAI}" ]]; then
  log "Indexing query assembly ${QRY_FA}..." | tee -a "${PREP_LOG}"
  samtools faidx "${QRY_FA}"
fi

flock -u 9

log "Processing ${chr} ..."
log "JobID=${SLURM_JOB_ID:-NA} NodeList=${SLURM_JOB_NODELIST:-NA} Threads=${THREADS} RetryThreads=${RETRY_THREADS}"

ref_chr_sam="${WORKDIR}/ref.${chr}.sam"
qry_chr_sam="${WORKDIR}/qry.${chr}.sam"

# Split each indexed BAM by chromosome.
samtools view -@ "${THREADS}" -h "${REF_BAM}" "${chr}" > "${ref_chr_sam}"
samtools view -@ "${THREADS}" -h "${QRY_BAM}" "${chr}" > "${qry_chr_sam}"

# Skip empty chromosome slices.
ref_alignments=$(samtools view -@ "${THREADS}" -c "${REF_BAM}" "${chr}")
qry_alignments=$(samtools view -@ "${THREADS}" -c "${QRY_BAM}" "${chr}")
if [[ "${ref_alignments}" == "0" || "${qry_alignments}" == "0" ]]; then
  log "Skipping ${chr}: one side has zero alignments (ref=${ref_alignments}, qry=${qry_alignments})"
  exit 0
fi

# --- Build a chromosome-specific, renamed query FASTA + SAM -----------------
# Supernova scaffold names carry no chromosome information, so AnchorWave has
# no way to know which chromosome a given query sequence belongs to. Here we:
#   1. Find scaffolds whose PRIMARY alignment (not secondary/supplementary)
#      lands on this chromosome, so each scaffold is claimed by at most one
#      chromosome.
#   2. Extract just those scaffolds from the full assembly and rename their
#      headers with a "${chr}_" prefix.
#   3. Rewrite the QNAME field of the chromosome SAM to match, so the SAM
#      passed to `-ar` and the FASTA passed to `-s` refer to the same names.
scaffold_ids="${WORKDIR}/${chr}.scaffold_ids.txt"
low_mapq_ids="${WORKDIR}/${chr}.low_mapq_scaffold_ids.txt"
qry_chr_fa="${WORKDIR}/${chr}.qry.renamed.fasta"
qry_chr_sam_renamed="${WORKDIR}/qry.${chr}.renamed.sam"

# MAPQ threshold for assigning a scaffold to this chromosome. Low MAPQ (e.g.
# MAPQ=1) means minimap2 found near-equally-good placements elsewhere, so
# "primary" is close to arbitrary; we don't want to silently anchor those
# scaffolds to a chromosome. Override with e.g. `MIN_MAPQ=30 sbatch ...`.
readonly MIN_MAPQ="${MIN_MAPQ:-20}"

samtools view -@ "${THREADS}" -F 0x904 -q "${MIN_MAPQ}" "${QRY_BAM}" "${chr}" \
  | cut -f1 | sort -u > "${scaffold_ids}"

# Track scaffolds whose primary alignment is here but didn't clear the MAPQ
# bar, so low-confidence placements are visible rather than silently dropped.
samtools view -@ "${THREADS}" -F 0x904 "${QRY_BAM}" "${chr}" \
  | awk -v mq="${MIN_MAPQ}" '$5 < mq {print $1}' | sort -u > "${low_mapq_ids}"

if [[ -s "${low_mapq_ids}" ]]; then
  log "Excluded $(wc -l < "${low_mapq_ids}") scaffold(s) from ${chr} with MAPQ < ${MIN_MAPQ} (ambiguous placement); see ${low_mapq_ids}"
fi

if [[ ! -s "${scaffold_ids}" ]]; then
  log "Skipping ${chr}: no scaffolds have a confident (MAPQ >= ${MIN_MAPQ}) primary alignment here"
  exit 0
fi

samtools faidx "${QRY_FA}" -r "${scaffold_ids}" \
  | awk -v chr="${chr}" '
      /^>/ { print ">" chr "_" substr($0, 2); next }
      { print }
    ' > "${qry_chr_fa}"

if [[ ! -s "${qry_chr_fa}" ]]; then
  log "ERROR: Failed to build renamed query FASTA for ${chr}" | tee -a "${LOGDIR}/${chr}.log"
  exit 1
fi

awk -v chr="${chr}" '
    /^@/ { print; next }
    { $1 = chr "_" $1; print }
  ' OFS='\t' "${qry_chr_sam}" > "${qry_chr_sam_renamed}"

out_anchor="${OUTDIR}/${chr}.anchors"
out_maf="${OUTDIR}/${chr}.maf"
out_frag_maf="${OUTDIR}/${chr}.fragmentation.maf"

# Run AnchorWave proali.
proali_args=(
  -t "${THREADS}"
  -i "${REF_GFF}"
  -r "${REF_FA}"
  -as "${CDS_FA}"
  -a "${ref_chr_sam}"
  -ar "${qry_chr_sam_renamed}"
  -s "${qry_chr_fa}"
  -n "${out_anchor}"
  -R 1
  -Q 1
  -o "${out_maf}"
  -f "${out_frag_maf}"
)

chr_log="${LOGDIR}/${chr}.log"
: > "${chr_log}"
log "Running primary attempt for ${chr} (threads=${THREADS})" | tee -a "${chr_log}"

set +e
"${ANCHORWAVE_BIN}" proali "${proali_args[@]}" >> "${chr_log}" 2>&1
primary_rc=$?
set -e

if [[ "${primary_rc}" -ne 0 ]]; then
  log "Primary attempt failed for ${chr} with exit code ${primary_rc}" | tee -a "${chr_log}"
  show_log_tail "${chr_log}"

  retry_anchor="${OUTDIR}/${chr}.retry.anchors"
  retry_maf="${OUTDIR}/${chr}.retry.maf"
  retry_frag_maf="${OUTDIR}/${chr}.retry.fragmentation.maf"
  retry_args=(
    -t "${RETRY_THREADS}"
    -i "${REF_GFF}"
    -r "${REF_FA}"
    -as "${CDS_FA}"
    -a "${ref_chr_sam}"
    -ar "${qry_chr_sam_renamed}"
    -s "${qry_chr_fa}"
    -n "${retry_anchor}"
    -R 1
    -Q 1
    -o "${retry_maf}"
    -f "${retry_frag_maf}"
  )

  log "Retrying ${chr} with safer settings (threads=${RETRY_THREADS})" | tee -a "${chr_log}"
  set +e
  "${ANCHORWAVE_BIN}" proali "${retry_args[@]}" >> "${chr_log}" 2>&1
  retry_rc=$?
  set -e

  if [[ "${retry_rc}" -ne 0 ]]; then
    log "Retry failed for ${chr} with exit code ${retry_rc}; see ${chr_log}" | tee -a "${chr_log}"
    show_log_tail "${chr_log}"
    exit "${retry_rc}"
  fi

  mv -f "${retry_anchor}" "${out_anchor}"
  mv -f "${retry_maf}" "${out_maf}"
  mv -f "${retry_frag_maf}" "${out_frag_maf}"
  log "Retry succeeded for ${chr}; promoted retry outputs" | tee -a "${chr_log}"
fi

log "Finished ${chr}"
log "All done."

