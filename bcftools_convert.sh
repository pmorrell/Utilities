#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH --ntasks=5
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

#  Created by Peter L Morrell - St. Paul, MN - 08 April 2025

set -e
set -o pipefail

# Load required modules
module load bcftools_ML_2/1.20

VCF="/scratch.global/pmorrell/Garden_Pea/Bandillo_299_peas.vcf.gz"
BCF="/scratch.global/pmorrell/Garden_Pea/Bandillo_299_peas.bcf"

# Load required modules

log() {
    local msg="$1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ${msg}"
}

log "-> Convert compressed VCF to BCF"
bcftools view "${VCF}" -Ob -o "${BCF}"

log "-> Indexing BCF file"
bcftools index -f "${BCF}"
log "-> Finished processing file"

