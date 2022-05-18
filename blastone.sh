#!/bin/bash
# script by Peter L. Morrell, original 11 August 2011, St. Paul, MN
# updated to use Morex_v2

#PBS -l mem=1000mb,nodes=1:ppn=1,walltime=48:00:00
#PBS -m abe
#PBS -M pmorrell@umn.edu

module load ncbi_blast+/2.8.1

# path to working directory
DIR=~/workshop/sandbox

# path to blast database
BLAST_DB=/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v2/Barley_Morex_V2_pseudomolecules.fasta
EXTENSION='txt'

cd $DIR

for i in *.fasta

do

# run blastn, but not megablast
#blastall -p blastn -n F \
# iterate over fasta files in directory

blastall -p blastn -n F -i $i -d ${BLAST_DB} -m 9 -o ${i}.${EXTENSION}
#-i $i \ 
# database to be used
#-d $BLAST_DB \
# output results in tabular format with header, XML is option 7
#-m 9 -o ${i}.${EXTENSION}

done

