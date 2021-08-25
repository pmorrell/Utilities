#!/usr/bin/env python3
"""
Read in large file listing the ancestral state of each SNP as 'AA', 'AD' or \
'DD'. This file is created using a VCF and ancestral state estimates, using \
`ancestral_state.py`. This script outputs frequency of derived variants in \
a diploid sample. The output of this script is the frequency of each \
polymorphism with a derived allele and a list of SNPs and genotypes at \
>95% frequency.

Usage: ./count_derived.py [anc] [out_dir]

Where:
1) [anc] is the file containing ancestral state information formatted as 'AA', 'AD', or 'DD'
2) [out_dir] is the full filepath to the output directory to store our output files
"""
#    A full test data set
#Chromosome      Pos     SNPID   Ancestral       Derived Reference       1046    1049    1050    1051    1053
#Chr01  1	Chr01_1	A       G       D       DD      DD      DD      DD
#Chr01  2	Chr01_2	A       G       D       DD      DD      DD      AA
#Chr01  3	Chr01_3	A       G       D       DD      DD      AA      AA
#Chr01  4	Chr01_4	A       G       D       AD      NN      NN      NN
#Chr01  5	Chr01_5	A       G       D       DD      AA      AA      AA
#Chr01	6	Chr01_6	A       G       D       AD      AA      AA      AA
#Chr01	7	Chr01_7	A       C       A       NN      NN      NN      NN

#   Peter L. Morrell - Falcon Heights, MN - 22 July 2021
#   Modified by Chaochih Liu - 10 August 2021

import sys
import os
import gzip


def allele_freq(genos):
    """Count frequency of the derived variant.
    """
    allele_list = []
    # Reformat genotypes into alleles, one allele per list element
    # This makes for easier math
    for elem in genos:
        # Don't include missing 'N' alleles
        if 'N' not in elem:
            for char in elem:
                allele_list.append(char)

    sample_size = len(allele_list)
    #   Outputs counts in the format: {'A': 1, 'D': 3}
    #       where the key is the allele, and the value is the count
    counts = {x: allele_list.count(x) for x in allele_list}
    der_count = counts['D']
    # Calcualte the frequency
    freq = der_count / sample_size
    return freq


def check_file_exists(fp):
    if os.path.exists(fp):
        # Remove file to start from clean slate
        os.remove(fp)
    return


def do_freq_counts(f, out_dir):
    """Call on functions defined, count allele frequency, and
    save to files."""
    fixation = {}
    der_dict = {}
    for line in f:
        if line.startswith('##'):
            continue
        elif line.startswith('Chromosome'):
            continue
        else:
            tmp = line.strip().split()
            snpid = tmp[2]
            geno = tmp[6:]
            # Only work with lines with at least one derived variant!
            if ('AD' in geno) or ('DD' in geno):
                # Call the frequency calculation function above.
                frequency = allele_freq(geno)
                der_dict[snpid] = [snpid, str(round(frequency, 4))]
                # Identify SNPs close to fixation
                if frequency > 0.9999:
                    fixation[snpid] = geno
    # Minimize the number of times output file is opened
    # Save derived frequency
    der_fp = out_dir + "/derived_freq.txt"
    with open(der_fp, 'a+') as f:
        for key in der_dict.keys():
            f.write('\t'.join(der_dict[key]) + '\n')
    # Save SNPs close to fixation
    fix_fp = out_dir + "/fixation_list.txt"
    with open(fix_fp, 'a+') as f:
        for key in fixation.keys():
            output_geno = '\t'.join(fixation[key])
            f.write('\t'.join([key, output_geno]) + '\n')
    return


def main(anc, out_dir):
    """Main function."""
    # Check if out_dir exists, if not make it
    out_dir_fp = os.path.expanduser(out_dir.rstrip('/'))
    if os.path.exists(out_dir_fp):
        print("Output directory exists, proceeding...")
    else:
        # Make output directory
        print("Output directory doesn't exist, creating directory. Proceeding...")
        os.makedirs(out_dir_fp)
    # Check if output file already exists, since we are appending we want
    #   to make sure we start from a clean file in case the file already exists
    check_file_exists(out_dir + "/derived_freq.txt")
    check_file_exists(out_dir + "/fixation_list.txt")
    # Then iterate through the derived VCF and print out the relevant fields
    if "gz" in anc:
        with gzip.open(anc, 'rt') as f:
            do_freq_counts(f, out_dir_fp)
    else:
        with open(anc, 'rt') as f:
            do_freq_counts(f, out_dir_fp)
    return


if len(sys.argv) <= 2:
    print(__doc__)
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
