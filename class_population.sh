#!//usr/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

#   A script for partitioning the derived sites files by variant class and population

#   Written by Peter L. Morrell
#   25 August 2021, Falcon Heights, MN

set -euf
set -o pipefail

#    The select samples script using pandas for paritioning sets of samples
module load python3/3.8.3_anaconda2020.07_mamba

#    Specific the locations of input files and outputs
WORKDIR=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Fagioli/Ancestral_state/SNPID/Variant_Classes
POPULATION_DIR=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Fagioli/Ancestral_state/SNPID/Groups
OUTFILE_BASE=_derived.txt.gz

SELECT=/panfs/roc/groups/9/morrellp/pmorrell/Apps/PMorrell/Utilities/select_samples.py


cd ${WORKDIR}
POPULATION_ARRAY=(AM_AND.txt.gz AM_MES.txt.gz EU_AND.txt.gz EU_MES.txt.gz EU_MIX.txt.gz)
CLASS_ARRAY=(deleterious.txt.gz stops.txt.gz synonymous.txt.gz tolerated.txt.gz)

for i in "${POPULATION_ARRAY[@]}"
 do
     for j in "${CLASS_ARRAY[@]}"
      do

   #    Strip the file name so we can create unique output
   POP=$(basename ${i} .txt)
   #    Strip the file name and add population name to each file
   CLASS=$(basename ${j} .txt.gz)
   python3 ${SELECT} ${POPULATION_DIR}/${i} ${j} 
   # | gzip > ./${POP}_${CLASS}${OUTFILE_BASE}
   #echo "${POP}"
   #echo "${CLASS}"
done
done
