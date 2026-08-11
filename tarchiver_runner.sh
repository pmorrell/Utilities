#!/bin/bash
#SBATCH --job-name=tarchiver
#SBATCH --output=tarchiver_%j.log
#SBATCH --error=tarchiver_%j.err
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-user=pmorrell@umn.edu



TARGET_DIR=$1
DEST_DIR="/scratch.global/pmorrell/User_archives"

if [ -z "$TARGET_DIR" ]; then
    echo "Error: no target directory provided"
    exit 1
fi

mkdir -p "$DEST_DIR"

python ~/Apps/Morrell_Lab/pmorrell/Utilities/tarchiver.py "$TARGET_DIR" "$DEST_DIR"

exit $?

