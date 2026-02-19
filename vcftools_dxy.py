#!/usr/bin/env python3

import pandas as pd
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Calculate Dxy from π values and pairwise diversity.")
parser.add_argument("--pop1", required=True, help="Path to the π values file for population 1 (pop1.sites.pi).")
parser.add_argument("--pop2", required=True, help="Path to the π values file for population 2 (pop2.sites.pi).")
parser.add_argument("--pairwise", required=True, help="Path to the pairwise diversity file (pairwise.weir.fst).")
parser.add_argument("--output", required=True, help="Path to the output file (e.g., dxy_results.txt).")
args = parser.parse_args()

# Load π values for each population
pop1_pi = pd.read_csv(args.pop1, sep="\t", names=["CHROM", "POS", "PI_POP1"], skiprows=1)
pop2_pi = pd.read_csv(args.pop2, sep="\t", names=["CHROM", "POS", "PI_POP2"], skiprows=1)

# Load pairwise diversity (Weir and Cockerham FST output)
pairwise_div = pd.read_csv(args.pairwise, sep="\t", names=["CHROM", "POS", "PAIRWISE_DIV"], skiprows=1)

# Merge the dataframes on CHROM and POS
merged = pd.merge(pop1_pi, pop2_pi, on=["CHROM", "POS"])
merged = pd.merge(merged, pairwise_div, on=["CHROM", "POS"])

# Calculate Dxy
merged["DXY"] = merged["PI_POP1"] + merged["PI_POP2"] + merged["PAIRWISE_DIV"]

# Save the results
merged.to_csv(args.output, sep="\t", index=False)