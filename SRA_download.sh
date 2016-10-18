#!/bin/env sh

#PBS -l mem=12gb,nodes=1:ppn=1,walltime=24:00:00 
#PBS -m abe 
#PBS -M mfrodrig@umn.edu 
#PBS -q lab

# Collect file from command line
#SRA_FILES=${HOME}/baseline/SRR_Acc_List-2.txt
SRA_FILES=$1

# Make sure the file exists
if ! [[ -f "${SRA_FILES}" ]]
    then echo "Failed to find ${SRA_FILES}, exiting..." >&2
    exit 1
    fi 

# Make the array using command substitution
# declare -a SRA_ARRAY=($(cat "${SRA_FILES}")) 
mapfile -t SRA_ARRAY < "$SRA_FILES"

#Print the values of the array to screen
printf '%s\n' "${SRA_ARRAY[@]}"
