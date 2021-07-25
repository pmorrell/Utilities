#!/usr/bin/env python
# Peter L. Morrell - Falcon Heights, MN - 22 July 2021

"""
Read in large file listing the ancestral state of each SNP as 'AA' or 'DD'. \
This file is created using a VCF and ancestral state estimates, using \
`ancestral_state.py`. This script outputs frequency of derived variants in \
diploid sample. The output of this script is frequency of each polymorphism \
with a derived allele.
"""

import sys
import gzip
import pandas as pd


def allele_freq(genos):
    """Count frequency of the derived variant."""
    # The line below uses a pandas function
    counts = dict(pd.value_counts(genos))
    # Missing genotypes need to be subtracted to get an accurate derived frequency.
    missing_key = 'NN'
    if missing_key in counts:
        observ = (len(genos) * 2) - (counts['NN'] * 2)
    else:
        observ = (len(genos) * 2)
    # Count heterozygotes so those genotypes also get an accurate count.
    het_key = 'AD'
    if het_key in counts:
        freq = (((counts['AD']/2) + counts['DD']) / observ)
    else:
        freq = ((counts['DD']) / observ)
    return freq


def main(anc):
    """Main function."""
    # Then iterate through the derived VCF and print out the relevant fields
    with gzip.open(anc, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('Chromosome'):
                continue
            else:
                tmp = line.strip().split()
                geno = tmp[6:]
                # Only work with lines with at least one derived variant!
                if ('AD' or 'DD' in geno):
                    # Call the frequency calculation function above.
                    frequency = allele_freq(geno)
                print(round(frequency, 3))
    return


if len(sys.argv) != 2:
    print("""Print the frequency of each variant with derived state information.
    Takes one argument:
    1) Ancestral state file (gzipped)""")
    exit(1)
else:
    main(sys.argv[1])
