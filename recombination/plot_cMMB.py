#!/usr/bin/env python
"""
Read smoothed cM/Mb recombination rates from a bed file.
Expected format: Chromosome, LeftBP, RightBP, cMMb_smoothed
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Default file path - can be overridden by command line argument
bed_file = sys.argv[1] if len(sys.argv) > 1 else "merged_genetic_physical_map.bed"

# Read the bed file
print(f"Reading recombination data from {bed_file}...")
df = pd.read_csv(bed_file, sep='\t', header=0)

# Keep only the columns we need
df = df[['Chromosome', 'LeftBP', 'RightBP', 'cMMb_lowess']].copy()
df.columns = ['Chromosome', 'LeftBP', 'RightBP', 'cMMb_smoothed']

# Convert numeric columns to proper types
df['LeftBP'] = pd.to_numeric(df['LeftBP'], errors='coerce')
df['RightBP'] = pd.to_numeric(df['RightBP'], errors='coerce')
df['cMMb_smoothed'] = pd.to_numeric(df['cMMb_smoothed'], errors='coerce')

# Filter out contig segments (keep only numbered chromosomes like Vu01, Vu02, etc.)
df = df[~df['Chromosome'].str.contains('contig', case=False, na=False)].copy()

# Calculate midpoint position
df['Position'] = (df['LeftBP'] + df['RightBP']) / 2
df['Size'] = df['RightBP'] - df['LeftBP']

print(f"Loaded {len(df)} recombination data points (after filtering contigs)")
print(f"\nChromosomes: {sorted(df['Chromosome'].unique())}")

# Summary statistics
print(f"\nRecombination rate statistics:")
print(f"  Min: {df['cMMb_smoothed'].min():.3f}")
print(f"  Max: {df['cMMb_smoothed'].max():.3f}")
print(f"  Mean: {df['cMMb_smoothed'].mean():.3f}")
print(f"  Median: {df['cMMb_smoothed'].median():.3f}")
print(f"  Std Dev: {df['cMMb_smoothed'].std():.3f}")

# Per-chromosome statistics
print(f"\nPer-chromosome statistics:")
for chrom in sorted(df['Chromosome'].unique()):
    chrom_data = df[df['Chromosome'] == chrom]
    print(f"  {chrom}: {len(chrom_data)} points, Mean cM/Mb = {chrom_data['cMMb_smoothed'].mean():.3f}")

# Visualize recombination rates by chromosome
print("\nGenerating visualization...")
chromosomes = sorted(df['Chromosome'].unique())
n_cols = min(3, len(chromosomes))
n_rows = (len(chromosomes) + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
if n_rows * n_cols == 1:
    axes = [axes]
else:
    axes = axes.flatten()

# Determine global y-axis scale
max_cmmb = 0.12

for i, chrom in enumerate(chromosomes):
    ax = axes[i]
    chrom_data = df[df['Chromosome'] == chrom].sort_values('Position')
    
    # Plot recombination rate
    ax.plot(chrom_data['Position']/1e6, chrom_data['cMMb_smoothed'], 
            color='#E63946', linewidth=2, marker='o', markersize=4, alpha=0.7)
    
    # Fill area under curve
    ax.fill_between(chrom_data['Position']/1e6, chrom_data['cMMb_smoothed'], 
                    alpha=0.2, color='#E63946')
    
    ax.set_title(f"{chrom}", fontsize=12)
    ax.set_xlabel('Position (Mb)', fontsize=10)
    ax.set_ylabel('cM/Mb', fontsize=10)
    ax.set_ylim(0, max_cmmb)
    ax.grid(True, alpha=0.3)

# Hide unused subplots
for i in range(len(chromosomes), len(axes)):
    axes[i].set_visible(False)

plt.suptitle('Smoothed Recombination Rates by Chromosome', fontsize=14)
plt.tight_layout()
plt.savefig('recombination_rates.png', dpi=300, bbox_inches='tight')
print("Saved plot to recombination_rates.png")

# Save summary data
output_file = 'recombination_summary.tsv'
df[['Chromosome', 'LeftBP', 'RightBP', 'Position', 'cMMb_smoothed']].to_csv(
    output_file, sep='\t', index=False
)
print(f"Saved summary data to {output_file}")

print("\nDone!")
