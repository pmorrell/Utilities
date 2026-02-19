#!/bin/bash

#   Peter L. Morrell - 28 Jan 2024 - Falcon Heights, MN
#   Convert VCF files to spreadsheet "Hudson table" 
#   format to run 2-site likelihood analysis.

USAGE="Usage:
$0 input_directory lookup_table output_directory [ -h ]
"

set -e
set -u
set -o pipefail

module load python3
module load bcftools/1.16-gcc-8.2.0-5d4xg4y

#   Chromosome name, necessary for partitioning
CHR=$1
#   Size of window in base pairs, SNP number will vary
WINDOW=$2
#   Size of chromosome being considered, e.g., barley chr1H is 516505932
CHR_LEN=$3
#   Full physical path to VCF
VCF=$4
#   Full physical path to temporary working directory
OUT_DIR=$5


# Determine number of windows given chromosome size.

seq 1 $WINDOW $CHR_LEN > $OUT_DIR/window_1.txt
seq $(($WINDOW)) $WINDOW $CHR_LEN > $OUT_DIR/window_2.txt
echo $CHR_LEN >> $OUT_DIR/window_2.txt
if $(wc -l window_1.txt) != $(wc -l window_2.txt) ; then
    exit 1
    echo "Window interval partitioning files of unequal length."
fi
#   Create a tab delimited file of the regions we want to use
yes $CHR: | head -n $CHR_LEN > $OUT_DIR/chr_name.txt
paste -d\0 <(yes $CHR: | head -n $CHR_LEN) \
<($OUT_DIR/window_1.txt $OUT_DIR/window_2.txt) >$OUT_DIR/regions.txt


while ((i++)); read -r line
do 
    bcftools view --regions $line $VCF -o $CHR_$((echo $i)).vcf

done <$OUT_DIR/regions.txt








