#!//usr/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu

# A script for partitioning the derived sites files by variant class

#   Written by Peter L. Morrell
#   25 August 2021, Falcon Heights, MN

set -euf
set -o pipefail

DERIVED=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Fagioli/Ancestral_state/SNPID/WGS_218_derived.txt.gz
OUTFILE_BASE=_derived.txt.gz
WORKDIR=/panfs/roc/groups/9/morrellp/pmorrell/Workshop/Fagioli/Ancestral_state/SNPID/Variant_Classes

TRIM=/panfs/roc/groups/9/morrellp/pmorrell/Apps/TKono/Deleterious_GP/Analysis_Scripts/Data_Handling/Trim_Ancestral.py


cd ${WORKDIR}
CLASS_ARRAY=(Fagioli_deleterious.anc_SNPID.txt Fagioli_stops.anc_SNPID.txt Fagioli_syn.anc_uniq_SNPID.txt Fagioli_tolerated.an_uniq_SNPID.txt)

for i in ${CLASS_ARRAY[@]}
 do
   #    Strip the file name so we can create unique output
   TEMP=$(basename ${i} .txt)
   #    Cut the SNPID from the full derived state file
   #    This will be slow, but should be possible in one run?
   python2 ${TRIM} ${i} ${DERIVED} | gzip > ./${TEMP}${OUTFILE_BASE}
done
