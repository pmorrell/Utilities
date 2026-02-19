#!/bin/bash
# written by Peter Morrell, April 2011, St. Paul, MN

# Calls two Perl scripts to subsample a large fasta file
# shuffle.pl is a short Perl script that reorders samples names in a flat file
# selectSeqs.pl is from Naoki Takebayashi, and selects only those samples named in the sample file
# tohudson2001 is a from Kevin Thornton analysis package from the libsequence library
# so many dependencies!

WORKING_DIR=~/Caballo/fasta/
SAMPLE_LIST=~/Caballo/fasta/GYS1_list.txt
SEQUENCE=~/Caballo/fasta/SequenceDataGYS1_8Dec10.fasta
EXHAP=~/Apps/Hudson/maxhap/exhap
MAXHAP=~/Apps/Hudson/maxhap/maxhap
# specify lookup table for maxhap
HNNNRHO=~/Apps/Hudson/maxhap/h100rho
REPS=100

for ((i=1;i<=$REPS;i++))
do
    shuffle.pl $SAMPLE_LIST | head -n 100 >rep${i}.list
    selectSeqs.pl -f rep${i}.list $SEQUENCE >rep${i}.fas
    tohudson2001 -i 1 rep${i}.fas | $EXHAP | $MAXHAP 1 $HNNNRHO 1 100 .01 0. 0. 0 250 >>$WORKING_DIR/rho.out
    rm rep${i}.list rep${i}.fas
done
