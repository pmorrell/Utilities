#!/usr/bin/env python
# Peter L. Morrell - Falcon Heights, MN - 06 August 2021
# Subset the large ancestral state table to calculate derived frequencies.

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


def file_name(f):
    base = os.path.basename(f)
    #    Remove the first extension
    base = os.path.splitext(base)[0]
    basename = os.path.splitext(base)[0]
    return basename


def main(samples, ancestral):
    """Main function."""
    # Then iterate through the derived VCF and print out the relevant fields
    retain = read_list(samples)
    df = pd.read_csv(ancestral, compression='gzip', header=0, sep='\t')
    df_out = df[["Chromosome", "Pos", "SNPID", "Ancestral", "Derived", "Reference"] + retain]
    sample, variant_class = file_name(samples, ancestral)
    variant_class = file_name(ancestral)
    sample_file = sample + '_' + variant_class + 'anc.txt.gz'
    df_out.to_csv(sample_file, sep='\t', index=False, compression='gzip')


if len(sys.argv) != 3:
    print("""Take a list of samples to cut down and the large file of \
    ancestral state for each individual and produces table with only the \
    listed samples retained.
    Takes one argument:
    1) List of samples to cut from data set
    2) Ancestral state file (gzipped)""")
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
