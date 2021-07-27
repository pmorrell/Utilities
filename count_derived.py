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
# import pandas as pd


def allele_freq(genos):
    freq = 0
    observ = 0
    """Count frequency of the derived variant."""
    # The line below uses a pandas function, this is a little fragile!
    # counts = dict(pd.value_counts(genos))
    counts = {x: genos.count(x) for x in genos}
# Missing genotypes need to be subtracted to get an accurate derived frequency.
    missing_key = ['NN']
    ancest_key = ['AA']
    het_key = ['AD']
    der_key = ['DD']
    for keys in counts:
        temp_keys = [*counts]
        #    If nothing but missing data, go to next SNP
        if temp_keys == missing_key:
            continue
        #    If nothing but ancestral state, go to next SNP
        elif temp_keys == ancest_key:
            continue
        #    If nothing but missing data and ancestral state, go to next SNP!
        #    The sum function flattens the nested lists
        elif temp_keys == sum([ancest_key, missing_key], []):
            continue
        else:
            #    Calculate sample sizes adjusting for missing data
            for missing_key[0] in counts:
                missing_count = (counts[missing_key[0]] * 2)
                #   Again, if nothing but missing data (shouldn't get here!)
                if missing_count == len(genos):
                    continue
                elif missing_count > 0:
                    observ = (len(genos) * 2) - (counts[missing_key[0]] * 2)
                elif counts[der_key[0]] == len(genos):
                    observ = (len(genos) * 2)
                else:
                    observ = (len(genos) * 2)
    # Count heterozygotes so those genotypes also get an accurate count.
        if het_key[0] in counts:
            if der_key[0] not in counts:
                freq = counts[het_key[0]] / observ
            elif der_key[0] in counts:
                freq = (counts[het_key[0]] + (counts[der_key[0]] * 2)) / observ
        else:
            freq = (counts[der_key[0]] * 2) / observ
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
                print(round(frequency, 4))
    return


if len(sys.argv) != 2:
    print("""Print the frequency of each variant with derived state information.
    Takes one argument:
    1) Ancestral state file (gzipped)""")
    exit(1)
else:
    main(sys.argv[1])
