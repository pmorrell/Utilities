#!/bin/bash -l                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                               
# Peter L. Morrell - 27 September 2024 - Falcon Heights, MN                                                                                                                                                                                        
                                                                                                                                                                                                                                               
set -e                                                                                                                                                                                                                                         
set -o pipefail                                                                                                                                                                                                                                
                                                                                                                                                                                                                                               
# Usage: Provide a list of inversion positions in a BED file. The 
# other required input is the named of a VCF file 
# indexed by bcftools. The third should be a path to the output directory.                                                                                                                                                                                                                         
                                                                                                                                                                                                                                               
module load bcftools/1.16-gcc-8.2.0-5d4xg4y

BED_FILE=$1
VCF=$2                                                                                                                                                                                                                                         
OUTDIR=$3
OUT_PREFIX=WBDC_GBS
SAMPLE=WBDC-355
MIN_INFORM=5


# Declare an array to store genomic positions    
declare -a GENOMIC_POSITIONS
   
                                                                                                                                                                                                           
while read -r line; do
  # Split the line into fields
  IFS=$'\t' read -r chrom start end <<< "$line"
  # Add 1 to the start and end positions
  new_start=$((start + 1))
  new_end=$((end + 1))
  # Combine chromosome, start, and end into a single string
  genomic_position="${chrom}:${new_start}-${new_end}"
  # Append the genomic position to the array
  GENOMIC_POSITIONS+=("$genomic_position")
done < "$BED_FILE"


# Print the genomic positions
for position in "${GENOMIC_POSITIONS[@]}"; do
  echo "$position"
done


if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR";
fi;

                                                                                                                                                                                                                                               
for i in "${GENOMIC_POSITIONS[@]}"                                                                                                                                                                                                                  
    do
    	INFORM_SNPS=$(bcftools query -r "$i" -s "$SAMPLE" -f '[%GT]\n' \
            "$VCF" | grep -E -c '1/1|0/1|1/0')
    		if "$INFORM_SNPS" gt "$MIN_INFORM"; then                                                                                                                                                                                                                            
        bcftools view -r "$i" "$VCF" -o "$OUTDIR"/"$OUT_PREFIX"_{"$i"}.vcf.gz -Oz                                                                                                                                                                                                     
            else
            echo "There are no SNPs in the interval $i"
        fi
    done

