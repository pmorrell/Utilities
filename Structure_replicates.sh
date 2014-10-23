#!/bin/bash -l
# written by Peter L. Morrell, April 2011, St. Paul, MN 
#PBS -l walltime=10:00:00,mem=4gb,nodes=4:ppn=1
#PBS -m abe -M youremail@email.edu

# the full path to user space can vary
# the line below makes is easy to set
Working_dir=/home/xe1/pmorrell/
# directory containing Structure application
DIR=$Working_dir/Apps/structure/console
# full path to application
PROGRAM=$Working_dir/Apps/structure2.3.3_console/structure
# Structure data file
INPUT1=$Working_dir/Cebada/OPA1_2/OPA1_2_landrace
# Structure main parameters file
INPUT2=$Working_dir/Cebada/OPA1_2/OPA1_2_landrace_mainparams
# Structure extra parameters file
INPUT3=$Working_dir/Cebada/OPA1_2/OPA1_2_landrace_extraparams
# output file; name will be appended 
OUTPUT=$Working_dir/Cebada/OPA1_2/OPA1_2_landrace

REPS=2
KMIN=2
KMAX=5

for ((j=$KMIN;j<=$KMAX;j++))
do

for ((i=1;i<=$REPS;i++))
do
    cd $DIR
    $PROGRAM -i $INPUT1 -m $INPUT2 -e $INPUT3 -K $j -o ${OUTPUT}_${j}_${i}
done
done
