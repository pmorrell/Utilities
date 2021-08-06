#!/usr/bin/env python
# Peter L. Morrell - Falcon Heights, MN - 06 August 2021

import os.path
import sys
import gzip
import pandas as pd

def read_list(s):
    sample_list = []
    with open(s, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                sample_list.append(line.strip())
    return sample_list

fn = sys.argv[1]
basename = os.path.basename(fn)
basename = os.path.splitext(basename)[0]
sample_file = "Ancestral_state_" + basename + ".txt.gz"

def main(samples, ancestral):
    """Main function."""
    # Then iterate through the derived VCF and print out the relevant fields
    retain = read_list(samples)
    df = pd.read_csv(ancestral, compression='gzip', header=0, sep='\t')
    df_out = df[["Chromosome", "Pos", "SNPID", "Ancestral", "Derived", "Reference"] + retain]
    df_out.to_csv(sample_file, sep='\t', index = False, compression = 'gzip')



if len(sys.argv) != 3:
    print("""Print the frequency of each variant with derived state information.
    Takes one argument:
    1) List of samples to cut from data set
    2) Ancestral state file (gzipped)""")
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
