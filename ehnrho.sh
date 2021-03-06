#!/bin/env bash

#PBS -l mem=2500mb,nodes=1:ppn=16,walltime=23:59:59
#PBS -m abe
#PBS -M pmorrell@umn.edu
#PBS -q lab

set -euo pipefail

#    Generate new lookup tables for Hudson's two-site likelihood approach
#    Peter L. Morrell, 26 September 2017
#    Dependencies: ehnrho.c http://home.uchicago.edu/rhudson1/source/twolocus/programs/ehnrho.html
EH=/panfs/roc/groups/9/morrellp/pmorrell/Apps/RHudson/ehnrho/eh

POP_NUMBER=1
SAMPLE_NUMBER=179
RHO_VALUES=25
GROWTH_RATE=0.00000
ITERATIONS=15000000

OUT_FILE=/panfs/roc/groups/9/morrellp/pmorrell/Apps/RHudson/ehnrho/h${SAMPLE_NUMBER}rho

$EH $POP_NUMBER $SAMPLE_NUMBER $RHO_VALUES \
0.000000 0.050000 0.100000 0.200000 0.300000 \
0.500000 0.750000 1.000000 1.500000 2.000000 \
3.000000 5.000000 7.500000 10.000000 15.000000 \
20.000000 30.000000 40.000000 60.000000 80.000000 \
100.000000 120.000000 150.000000 175.000000 200.00000 \
$GROWTH_RATE $ITERATIONS > ${OUT_FILE}

