#!/bin/env sh

#PBS -l mem=12gb,nodes=1:ppn=1,walltime=24:00:00 
#PBS -m abe 
#PBS -M pmorrell@umn.edu 
#PBS -q lab

#    Collect file from command line
#SRA_FILES=${HOME}/scratch/MBE_Barley/MBE_barley.txt
SRA_FILES=${HOME}/scratch/MBE_outgroup/MBE_outgroup_samples.txt

#   Make sure the file exists
if ! [[ -f "${SRA_FILES}" ]]
    then echo "Failed to find ${SRA_FILES}, exiting..." >&2
    exit 1
    fi 

#   Make the array using command substitution
declare -a SRA_ARRAY=($(cat "${SRA_FILES}")) 

#   Print the values of the array to screen
printf '%s\n' "${SRA_ARRAY[@]}"

#   location of SRA download script, from Tom Kono's Misc_Utils
SRA_FETCH=${HOME}/Apps/TKono/Misc_Utils/SRA_Fetch.sh 

#   directory where SRA files will be downloaded
OUTPUT=${HOME}/scratch/MBE_outgroup/

#   iterate over every each of the run numbers in a lit of SRA files
#   and download to specified directory
#   in SRA_Fetch -r = run #, -e experiment #, -p sample #, -s study #
for i in "${SRA_ARRAY[@]}"
    do bash $SRA_FETCH -r $i -d $OUTPUT
    done

