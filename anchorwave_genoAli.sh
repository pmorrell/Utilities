#!/usr/bin/env bash
#SBATCH  --partition=msibigmem 
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=15900M
#SBATCH --tmp=300g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o anchorwave_%j.out
#SBATCH -e anchorwave_%j.err

set -euo pipefail

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(timestamp)] $*"; }

module load samtools/1.21

if [[ -x "${HOME}/Apps/anchorwave" ]]; then
  ANCHORWAVE_BIN="${HOME}/Apps/anchorwave/anchorwave"
elif [[ -x "${HOME}/Apps/anchorwave/anchorwave" ]]; then
  ANCHORWAVE_BIN="${HOME}/Apps/anchorwave/anchorwave"
else
  echo "ERROR: Could not find executable AnchorWave at ~/Apps/anchorwave or ~/Apps/anchorwave/anchorwave" >&2
  exit 1
fi

# Use the SLURM CPU allocation when available, otherwise default to 1 thread.
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Retry settings for instability/node interrupts.
RETRY_THREADS="${RETRY_THREADS:-4}"

# Inversion calling can help when assemblies differ by scaffold orientation.
# Set ENABLE_INVERSIONS="false" to disable -IV for comparison runs.
ENABLE_INVERSIONS="${ENABLE_INVERSIONS:-true}"

# Inputs
REF_SAM="/scratch.global/pmorrell/Morex_v3/ref.sam"
QRY_SAM="/scratch.global/pmorrell/Morex_v3/cds.sam"

# Keep these as your existing files/args from the single-run pipeline
REF_FA="/projects/standard/morrellp/pmorrell/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules.fasta"
QRY_FA="/projects/standard/morrellp/pmorrell/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules.fasta"
REF_GFF="/projects/standard/morrellp/pmorrell/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.sorted.gff3"
CDS_FA="/scratch.global/pmorrell/Morex_v3/cds.fa"

# Output folders
WORKDIR="per_chr_work"
OUTDIR="per_chr_out"
LOGDIR="per_chr_logs"
mkdir -p "${WORKDIR}" "${OUTDIR}"
mkdir -p "${LOGDIR}"

# Region queries require indexed BAM/CRAM (not plain SAM).
REF_BAM="${WORKDIR}/ref.sorted.bam"
QRY_BAM="${WORKDIR}/qry.sorted.bam"

if [[ ! -s "${REF_BAM}" || ! -s "${REF_BAM}.csi" ]]; then
  samtools sort -@ "${THREADS}" -o "${REF_BAM}" "${REF_SAM}"
  samtools index -@ "${THREADS}" --csi "${REF_BAM}"
fi

if [[ ! -s "${QRY_BAM}" || ! -s "${QRY_BAM}.csi" ]]; then
  samtools sort -@ "${THREADS}" -o "${QRY_BAM}" "${QRY_SAM}"
  samtools index -@ "${THREADS}" --csi "${QRY_BAM}"
fi

log "Starting AnchorWave genoAli per chromosome"
log "JobID=${SLURM_JOB_ID:-NA} NodeList=${SLURM_JOB_NODELIST:-NA} Threads=${THREADS} ENABLE_INVERSIONS=${ENABLE_INVERSIONS}"

for i in {1..7}; do
  chr="chr${i}H"
  log "Processing ${chr} ..."

  # 1) Split each indexed BAM by chromosome (header kept with -h)
  ref_chr_sam="${WORKDIR}/ref.${chr}.sam"
  qry_chr_sam="${WORKDIR}/qry.${chr}.sam"
  samtools view -@ "${THREADS}" -h "${REF_BAM}" "${chr}" > "${ref_chr_sam}"
  samtools view -@ "${THREADS}" -h "${QRY_BAM}" "${chr}" > "${qry_chr_sam}"

  # Optional safety check: skip empty chromosome slices
  ref_alignments=$(samtools view -@ "${THREADS}" -c "${REF_BAM}" "${chr}")
  qry_alignments=$(samtools view -@ "${THREADS}" -c "${QRY_BAM}" "${chr}")
  if [[ "${ref_alignments}" == "0" || "${qry_alignments}" == "0" ]]; then
    echo "Skipping ${chr}: one side has zero alignments (ref=${ref_alignments}, qry=${qry_alignments})"
    continue
  fi

  # 2) Run AnchorWave per chromosome
  # Map flags to genoAli usage:
  # -i ref GFF, -r ref genome, -as anchor cds.fa,
  # -a query SAM, -ar reference SAM, -s target genome,
  # -n output anchors, -o output MAF, -f fragmentation MAF.
  out_anchor="${OUTDIR}/${chr}.anchors"
  out_maf="${OUTDIR}/${chr}.maf"
  out_frag_maf="${OUTDIR}/${chr}.fragmentation.maf"
  genoali_args=(
    -t "${THREADS}"
    -i "${REF_GFF}"
    -r "${REF_FA}"
    -as "${CDS_FA}"
    -a "${qry_chr_sam}"
    -ar "${ref_chr_sam}"
    -s "${QRY_FA}"
    -n "${out_anchor}"
    -o "${out_maf}"
    -f "${out_frag_maf}"
  )

  if [[ "${ENABLE_INVERSIONS}" == "true" ]]; then
    genoali_args+=( -IV )
  fi

  chr_log="${LOGDIR}/${chr}.log"
  : > "${chr_log}"
  log "Running primary attempt for ${chr} (threads=${THREADS}, inversions=${ENABLE_INVERSIONS})" | tee -a "${chr_log}"

  set +e
  "${ANCHORWAVE_BIN}" genoAli "${genoali_args[@]}" >> "${chr_log}" 2>&1
  primary_rc=$?
  set -e

  if [[ "${primary_rc}" -ne 0 ]]; then
    log "Primary attempt failed for ${chr} with exit code ${primary_rc}" | tee -a "${chr_log}"

    # Retry once with reduced threads and without inversions to reduce peak pressure.
    retry_anchor="${OUTDIR}/${chr}.retry.anchors"
    retry_maf="${OUTDIR}/${chr}.retry.maf"
    retry_frag_maf="${OUTDIR}/${chr}.retry.fragmentation.maf"
    retry_args=(
      -t "${RETRY_THREADS}"
      -i "${REF_GFF}"
      -r "${REF_FA}"
      -as "${CDS_FA}"
      -a "${qry_chr_sam}"
      -ar "${ref_chr_sam}"
      -s "${QRY_FA}"
      -n "${retry_anchor}"
      -o "${retry_maf}"
      -f "${retry_frag_maf}"
    )

    log "Retrying ${chr} with safer settings (threads=${RETRY_THREADS}, inversions=false)" | tee -a "${chr_log}"
    set +e
    "${ANCHORWAVE_BIN}" genoAli "${retry_args[@]}" >> "${chr_log}" 2>&1
    retry_rc=$?
    set -e

    if [[ "${retry_rc}" -ne 0 ]]; then
      log "Retry failed for ${chr} with exit code ${retry_rc}; see ${chr_log}" | tee -a "${chr_log}"
      exit "${retry_rc}"
    fi

    mv -f "${retry_anchor}" "${out_anchor}"
    mv -f "${retry_maf}" "${out_maf}"
    mv -f "${retry_frag_maf}" "${out_frag_maf}"
    log "Retry succeeded for ${chr}; promoted retry outputs" | tee -a "${chr_log}"
  fi

  log "Finished ${chr}"
done

log "All done."
