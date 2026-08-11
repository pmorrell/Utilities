#!/bin/bash

#SBATCH --job-name=HPC_archiver
#SBATCH --output=HPC_archiver_%j.log
#SBATCH --error=HPC_archiver_%j.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

set -euo pipefail

if [[ $# -lt 1 ]]; then
	echo "Usage: ./HPC_archiver.sh <source_dir> [staging_root] [archive_root]"
	echo "   or: sbatch --job-name=<name> HPC_archiver.sh <source_dir> [staging_root] [archive_root]"
	exit 2
fi

# Slurm parses #SBATCH before shell variables exist, so derive job name and
# submit with --job-name when this script is launched outside an allocation.
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
	SOURCE_BASE_SUBMIT="$(basename "${1%/}")"
	JOB_NAME="$(echo "$SOURCE_BASE_SUBMIT" | tr -cs '[:alnum:]_.-' '_' | sed 's/^_\+//; s/_\+$//' )"
	if [[ -z "$JOB_NAME" ]]; then
		JOB_NAME="HPC_archiver"
	fi
	echo "Submitting Slurm job as: $JOB_NAME"
	sbatch --job-name="$JOB_NAME" "$0" "$@"
	exit 0
fi

SOURCE_DIR="${1%/}"
STAGING_ROOT="${2:-/scratch.global/pmorrell/User_archives}"
ARCHIVE_ROOT="${3:-$STAGING_ROOT}"
ZSTD_LEVEL="${ZSTD_LEVEL:-10}"

if ! command -v rsync >/dev/null 2>&1; then
	echo "Error: rsync is not available"
	exit 2
fi
if ! command -v tar >/dev/null 2>&1; then
	echo "Error: tar is not available"
	exit 2
fi
if ! command -v zstd >/dev/null 2>&1; then
	echo "Error: zstd is not available"
	exit 2
fi
if [[ ! -d "$SOURCE_DIR" ]]; then
	echo "Error: source directory does not exist: $SOURCE_DIR"
	exit 2
fi

SOURCE_BASE="$(basename "$SOURCE_DIR")"
STAMP="$(date +%Y%m%d_%H%M%S)"
JOB_TAG="${SLURM_JOB_ID:-manual}"

mkdir -p "$STAGING_ROOT" "$ARCHIVE_ROOT"

STAGED_DIR="$STAGING_ROOT/${SOURCE_BASE}_${JOB_TAG}_${STAMP}"
ARCHIVE_PATH="$ARCHIVE_ROOT/${SOURCE_BASE}_${JOB_TAG}_${STAMP}.tar.zst"
VERIFY_LOG="$(mktemp)"
trap 'rm -f "$VERIFY_LOG"' EXIT

echo "Source:  $SOURCE_DIR"
echo "Staged:  $STAGED_DIR"
echo "Archive: $ARCHIVE_PATH"

echo "Step 1/5: rsync copy to staging"
rsync -aHAX --numeric-ids --info=progress2 "$SOURCE_DIR/" "$STAGED_DIR/"

echo "Step 2/5: checksum verify staged copy"
rsync -aHAXcn --delete "$SOURCE_DIR/" "$STAGED_DIR/" > "$VERIFY_LOG"
if [[ -s "$VERIFY_LOG" ]]; then
	echo "Verification failed. Differences detected between source and staged copy:"
	cat "$VERIFY_LOG"
	exit 2
fi

SRC_COUNT="$(find "$SOURCE_DIR" -type f | wc -l | tr -d ' ')"
DST_COUNT="$(find "$STAGED_DIR" -type f | wc -l | tr -d ' ')"
echo "File counts: source=$SRC_COUNT staged=$DST_COUNT"
if [[ "$SRC_COUNT" != "$DST_COUNT" ]]; then
	echo "Verification failed. File counts do not match."
	exit 2
fi

echo "Step 3/5: remove source after successful verification"
rm -rf "$SOURCE_DIR"

echo "Step 4/5: tar + zstd compress"
tar -C "$STAGING_ROOT" -cf - "$(basename "$STAGED_DIR")" | zstd -T0 "-$ZSTD_LEVEL" -q -o "$ARCHIVE_PATH"

echo "Step 5/5: verify archive readability"
zstd -t -q "$ARCHIVE_PATH"
tar -I zstd -tf "$ARCHIVE_PATH" >/dev/null

echo "Completed successfully"
echo "Archive created: $ARCHIVE_PATH"
