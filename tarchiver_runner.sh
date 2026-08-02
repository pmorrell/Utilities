#!/bin/bashh
#SBATCH --job-name=tarchiver
#SBATCH --output=tarchiver_%j.log
#SBATCH --error=tarchiver_%j.err
#SBATCH --time=03:00:00        
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      
#SBATCH --mem=4G               

# $1 captures the first argument passed to the sbatch command
TARGET_DIR=$1
DEST_DIR="/scratch.global/pmorrell/User_archives"

# Run the python script using the argument
python ~/Apps/Morrell_Lab/pmorrell/Utilities/tarchiver.py "$TARGET_DIR" "$DEST_DIR"
