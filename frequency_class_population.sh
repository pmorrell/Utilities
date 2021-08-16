#!//usr/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

#   A script for partitioning the derived sites files by variant class and population
#   Provide the list of SNPIDs by class and the ancestral state file per population
#   Get back an ancestral state file useful for frequency calcualation

#   Written by Peter L. Morrell
#   16 August 2021, Falcon Heights, MN

set -euf
set -o pipefail

#    For example: /panfs/roc/groups/9/morrellp/pmorrell/Workshop/Fagioli/Ancestral_state/SNPID/Fagioli_stops.anc_SNPID.txt
SNPID=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Fagioli/Ancestral_state/SNPID/Fagioli_stops.anc_SNPID.txt
OUTFILE_BASE=_stops.txt
WORKDIR=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Fagioli/Ancestral_state/SNPID/

TRIM=/panfs/roc/groups/9/morrellp/pmorrell/Apps/TKono/Deleterious_GP/Analysis_Scripts/Data_Handling/Trim_Ancestral.py


cd ${WORKDIR}

DERIVED_ARRAY=(Ancestral_state_Fagioli_A1.txt.gz  Ancestral_state_Fagioli_A2.txt.gz  Ancestral_state_Fagioli_A3.txt.gz  Ancestral_state_Fagioli_M1.txt.gz  Ancestral_state_Fagioli_M2.txt.gz)

for i in "${DERIVED_ARRAY[@]}"ß
 do
   TEMP=$(basename ${i} .txt.gz)
   #    Get the population name off the end of the files
   OUT=$(echo ${TEMP} | awk -F "_" '{print $NF}')
   #    Each of the directories is named after the population
   python2 ${TRIM} ${SNPID} ${i} > ./${OUT}/${OUT}${OUTFILE_BASE}
done
