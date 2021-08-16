#!/bin/bash -l

#     Calculate derived allele count. Can be applied to each class of varaints.

#    Written by Peter L. Morrell
#    09 August 2021, Falcon Heights, MN

set -euf
set -o pipefail

INFILE=$1
OUTFILE=$(basename ${INFILE} .txt.gz)
WIDTH=$(head -n1 ${INFILE} | wc -l )

#REFERENCE=$(datamash --header-in --sort --g6 count 6 <(zcat ${INFILE})) \
#> ${OUTFILE}_counts.txt

#
for i in $(seq 6 $WIDTH)
  do zcat < ${INFILE} | datamash --sort -g${i} count ${i} \
    >> ${OUTFILE}_counts.txt
  done
