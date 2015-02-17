#!/usr/bin/env python
# 

from __future__ import print_function
import sys

#   Set missing data value
missing = 'NA'
#   What are the valid nucleotide states?
valid_bases = set(['A', 'C', 'T', 'G'])

def GetMajorMinor(snp):
    #   A function to get the major and minor alleles of a SNP
    #   Returns of the form (minor_allele, major_allele)
    #   Remove missing calls
    snp = [x for x in snp if x != missing]
    #   Next, join them all together into a large string
    all_calls = ''.join(snp)
    #   And then cast to set
    #   This saves only the unique elements, i.e., the alleles
    alleles = set(all_calls)
    #   We want to remove any "alleles" that aren't actually bases
    #   we can do this by intersecting with the set of valid bases that we have
    #   defined above.
    alleles = alleles.intersection(valid_bases)
    #   Then, get the counts of each allele
    #   Store this in a dictionary so we can associate the allele with its count
    counts = {}
    for a in alleles:
        counts[a] = all_calls.count(a)
    #   Then get the minimum one
    #   This returns the element of counts that has the smallest count
    #   it's similar to the key= in sort; it does the calculations on the
    #   allele counts, but returns the associated key in the dictionary
    minor_allele = min(counts, key=counts.get)
    #   The same for major allele
    major_allele = max(counts, key=counts.get)
    #   We get weird results when the alleles are at equal frequency, so in
    #   this case, we just arbitrarily set them to major or minor
    #   We have to cast back to list for this since sets do not support slicing
    if minor_allele == major_allele:
        minor_allele = list(alleles)[0]
        major_allele = list(alleles)[1]
    #   Return the data structure
    return(minor_allele, major_allele)


file_data = []
#   Read the file in line-by-line
with open(sys.argv[1]) as f:
    for line in f:
        #   Skip the header lines - write them out without modification
        if line.startswith('#'):
            sys.stdout.write(line)
        else:
            file_data.append(line.strip().split('\t'))

# Remove header from matrix            
header = file_data.pop(0)

#   Transpose so we can iterate over columns
file_data_t = zip(*file_data)
#   We can get sample names as the first element
sample_names = file_data_t[0]

# 
geno_output = []

#   And iterate over columns
for snp in file_data_t[1:]:
    snp = map(str.upper,(snp))
    minor, major = GetMajorMinor(snp)
    #   Unfortunately, python doesn't have an elegant way to
    #   replace all values in a list like R does
    #   Start with minor
    #       For diagnostic purposes, print the SNP before and its minor and major states
    # print snp, minor, major
    snp = ['0' if (x == minor*2) else x for x in snp]
    #   Then major
    snp = ['2' if (x == major*2) else x for x in snp]
    #   And finally the hets
    snp = ['1' if (x == major + minor or x == minor + major) else x for x in snp]
    geno_output.append(snp)

geno_output_t = zip(*geno_output)
#for geno in geno_output_t:
#	print geno
print (*geno_output_t, sep='\n')	



