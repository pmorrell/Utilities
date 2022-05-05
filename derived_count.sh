#!/bin/bash -l

#    Calculate derived allele count. Can be applied to each class of varaints.
#    Written for counting derived variants 

#    Written by Peter L. Morrell
#    09 August 2021, Falcon Heights, MN

set -euf
set -o pipefail

INFILE=$1



REFERENCE=$(datamash --sort --group6 count 6 < ${INFILE)

#
for i in {7..226}
    do datamash --sort --group${i} count ${i} < ${INFILE}
    >> ${OUTFILE}
