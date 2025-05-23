#!/usr/bin/env python

"""
Pull ancestral state of variants from a list and create input for plotting the \
unfolded SFS and identifying ancestral state for each genotype at each position.
This script is lightly rewritten from a version from Tom Kono written for \
barley ancestral state inference: \
https://github.com/MorrellLAB/Deleterious_GP/blob/master/Analysis_Scripts/Data_Handling/H_murinum_Ancestral.py

"""

import sys
import gzip


#    Run through the ancestral state files and record the chromosome positions
#    as a key, and ancestral state as a value
def store_anc(a):
    """Read through the list of position and ancestral states and store the ancestral state. This will be in a
    dictionary keyed on chr:pos."""
    anc = {}
    with gzip.open(a, 'rt') as f:
        for line in f:
            if line.startswith('chr'):
                continue
            else:
                tmp = line.strip().split('\t')
                k = tmp[0] + ':' + tmp[1]
                v = tmp[2]
                anc[k] = v
    return anc

def polarize(ref, alt, anc, samples):
    """Return a list of the polarized genotypes in the samples. We will use
    A for ancestral allele, D for derived allele, and N for unknown. The
    genotypes will be diploid."""
    pol = []
    # First, if ancestral is N, then we return N for all samples
    if anc == 'N':
        der = 'N'
        for s in samples:
            pol.append('NN')
        return (pol, der)
    else:
        if anc == ref:
            der = alt
            a = '0'
            d = '1'
        else:
            der = ref
            a = '1'
            d = '0'
        for s in samples:
            gt = s.split(':')[0]
            if gt == '.' or gt == './.':
                pol.append('NN')
            elif gt == a + '/' + a:
                pol.append('AA')
            elif (gt == a + '/' + d) or (gt == d + '/' + a):
                pol.append('AD')
            elif gt == d + '/' + d:
                pol.append('DD')
            else:
                pol.append('NN')
        return (pol, der)



def main(anc, snps_vcf):
    """Main function."""
    # Parse and store the ancestral states in the H. hurinum VCF
    astates = store_anc(anc)
    # Then iterate through the derived VCF and print out the relevant fields
    with gzip.open(snps_vcf, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                samples = line.strip().split('\t')[9:]
                # Print the header here
                print('\t'.join(['Chromosome', 'Pos', 'SNPID', 'Ancestral', 'Derived', 'Reference'] + samples))
            else:
                tmp = line.strip().split()
                key = tmp[0] + ':' + tmp[1]
                snpid = tmp[0] + '_' + tmp[1]
                ref = tmp[3]
                alt = tmp[4]
                samp = tmp[9:]
                ancestral = astates.get(key, 'N')
                if ancestral == 'N':
                    morex = 'N'
                elif ancestral == ref:
                    morex = 'A'
                else:
                    morex = 'D'
                genos, derived = polarize(ref, alt, ancestral, samp)
                # Print out the chrom, pos, ancestral, derived, morex state, and
                # the polarized letters for the samples
                print('\t'.join([tmp[0], tmp[1], snpid, ancestral, derived, morex] + genos))
    return


if len(sys.argv) != 3:
    print("""Print out the ancestral/derived alleles for each SNP, and which samples have
the derived or ancestral alleles at each SNP. Both files must be anchored on
the pseudomolecule assembly. Takes two arguments:
    1) Ancestral state list (gzipped)
    2) Non-ancestral VCF (gzipped)""")
    sys.exit(1)
else:
    main(sys.argv[1], sys.argv[2])
