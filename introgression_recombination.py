import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import pearsonr, spearmanr
import os

# File paths
recomb_file = "/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Analyses/recombination/bopa_9k_cM_100kb.txt"
introgression_file = "/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Analyses/recombination/interval_genetic_distances.tsv"
output_dir = os.path.dirname(recomb_file)

# Read recombination data
print("Reading recombination data...")
df_recomb = pd.read_csv(recomb_file, sep='\t')

# Calculate midpoint position for each genomic region
df_recomb['Position'] = (df_recomb['LeftBP'] + df_recomb['RightBP']) / 2
df_recomb['Size'] = df_recomb['RightBP'] - df_recomb['LeftBP']

# Set negative values to zero as requested
df_recomb['cMMb_nonneg'] = df_recomb['cMMb'].clip(lower=0)

# Filter extreme outliers using IQR method
Q1 = df_recomb['cMMb'].quantile(0.25)
Q3 = df_recomb['cMMb'].quantile(0.75)
IQR = Q3 - Q1
lower_bound = Q1 - 5*IQR
upper_bound = Q3 + 5*IQR
df_filtered = df_recomb[(df_recomb['cMMb'] >= lower_bound) & (df_recomb['cMMb'] <= upper_bound)].copy()
print(f"Filtered out {len(df_recomb) - len(df_filtered)} extreme values")

# Apply LOWESS smoothing by chromosome to get smoothed recombination rates
smoothed_data = {}
for chrom in df_filtered['Chromosome'].unique():
    chrom_data = df_filtered[df_filtered['Chromosome'] == chrom].sort_values('Position')
    
    if len(chrom_data) > 3:
        # Apply LOWESS smoothing
        smoothed = lowess(
            chrom_data['cMMb_nonneg'], 
            chrom_data['Position'],
            frac=0.15,  # Smoothing parameter
            it=3,       # Number of iterations
            return_sorted=True
        )
        
        smoothed_df = pd.DataFrame({
            'Chromosome': chrom,
            'Position': smoothed[:, 0],
            'cMMb_smoothed': smoothed[:, 1]
        })
        smoothed_data[chrom] = smoothed_df

# Combine all smoothed data
df_smoothed = pd.concat(smoothed_data.values())
print(f"Generated smoothed recombination rates for {len(df_smoothed)} positions")

# Save smoothed recombination data
smoothed_output_path = os.path.join(output_dir, "recombination_smoothed_data.tsv")
df_smoothed.to_csv(smoothed_output_path, sep='\t', index=False)
print(f"Saved smoothed recombination data to: {smoothed_output_path}")

# Read introgression data
print("Reading introgression data...")
df_introgression = pd.read_csv(introgression_file, sep='\t')
print(f"Read {len(df_introgression)} introgression segments")

# Get total number of unique samples
total_samples = df_introgression['Sample'].nunique()
print(f"Found {total_samples} unique samples")

# Create genomic bins for each chromosome (using recombination position data as anchors)
print("Creating genomic bins for summarizing introgression...")
bin_size = 1_000_000  # 1 Mb bins
chromosome_bins = {}

for chrom in df_smoothed['Chromosome'].unique():
    # Get chromosome length from max position in introgression data
    chrom_introgression = df_introgression[df_introgression['Chromosome'] == chrom]
    if len(chrom_introgression) == 0:
        continue
    
    max_pos = chrom_introgression['End'].max()
    bins = np.arange(0, max_pos + bin_size, bin_size)
    chromosome_bins[chrom] = bins

# Count introgression events per bin
print("Calculating introgression frequency...")
introgression_summary = []

for chrom, bins in chromosome_bins.items():
    # Initialize counters for each bin
    bin_counts = np.zeros(len(bins) - 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    
    # Get introgression segments for this chromosome
    chrom_introgression = df_introgression[df_introgression['Chromosome'] == chrom]
    
    # Count unique samples with introgression in each bin
    for i in range(len(bins) - 1):
        bin_start = bins[i]
        bin_end = bins[i+1]
        
        # Find samples with introgression overlapping this bin
        samples_in_bin = set()
        for _, row in chrom_introgression.iterrows():
            if row['End'] >= bin_start and row['Start'] <= bin_end:
                samples_in_bin.add(row['Sample'])
        
        # Count unique samples and normalize
        bin_counts[i] = len(samples_in_bin) / total_samples
    
    # Create summary dataframe
    chrom_summary = pd.DataFrame({
        'Chromosome': chrom,
        'Position': bin_centers,
        'IntrogressionRate': bin_counts
    })
    
    introgression_summary.append(chrom_summary)

# Combine all chromosome summaries
df_introgression_rate = pd.concat(introgression_summary)
print(f"Generated introgression rates for {len(df_introgression_rate)} genomic bins")

# Save introgression summary
intro_output_path = os.path.join(output_dir, "introgression_summary.tsv")
df_introgression_rate.to_csv(intro_output_path, sep='\t', index=False)
print(f"Saved introgression summary to: {intro_output_path}")

# Align recombination and introgression data
print("Aligning recombination and introgression data...")
merged_data = []

for chrom in df_smoothed['Chromosome'].unique():
    recomb_chrom = df_smoothed[df_smoothed['Chromosome'] == chrom]
    introgression_chrom = df_introgression_rate[df_introgression_rate['Chromosome'] == chrom]
    
    if len(recomb_chrom) == 0 or len(introgression_chrom) == 0:
        continue
    
    # For each recombination position, find nearest introgression bin
    for _, recomb_row in recomb_chrom.iterrows():
        recomb_pos = recomb_row['Position']
        
        # Find closest introgression bin
        distances = abs(introgression_chrom['Position'] - recomb_pos)
        nearest_idx = distances.idxmin()
        nearest_bin = introgression_chrom.loc[nearest_idx]
        
        # Only include if positions are within reasonable distance (2Mb)
        if abs(nearest_bin['Position'] - recomb_pos) <= 2_000_000:
            merged_data.append({
                'Chromosome': chrom,
                'Position': recomb_pos,
                'cMMb_smoothed': recomb_row['cMMb_smoothed'],
                'IntrogressionRate': nearest_bin['IntrogressionRate']
            })

df_merged = pd.DataFrame(merged_data)
print(f"Created {len(df_merged)} matched data points for correlation")

# Save merged data
merged_output_path = os.path.join(output_dir, "recombination_introgression_corr.tsv")
df_merged.to_csv(merged_output_path, sep='\t', index=False)
print(f"Saved correlation data to: {merged_output_path}")

# Calculate correlations
print("Calculating correlations...")
pearson_r, pearson_p = pearsonr(df_merged['cMMb_smoothed'], df_merged['IntrogressionRate'])
spearman_r, spearman_p = spearmanr(df_merged['cMMb_smoothed'], df_merged['IntrogressionRate'])

print(f"Pearson correlation: r={pearson_r:.3f}, p={pearson_p:.6f}")
print(f"Spearman correlation: r={spearman_r:.3f}, p={spearman_p:.6f}")

# Replace the current plotting section with this code (after line ~181)

# Create chromosome-wise profile plots
print("Creating chromosome profile plots...")

# Get ordered list of chromosomes
chromosomes = sorted(df_smoothed['Chromosome'].unique())

# Set up figure with subplots - one per chromosome
plt.figure(figsize=(18, 14))
n_cols = min(2, len(chromosomes))
n_rows = (len(chromosomes) + n_cols - 1) // n_cols

# Define colors for better distinction
recomb_color = "#E63946"  # red for recombination
introgression_color = "#457B9D"  # blue for introgression

# Store chromosome-specific correlations for summary
correlations = []

# Create a subplot for each chromosome
for i, chrom in enumerate(chromosomes):
    ax = plt.subplot(n_rows, n_cols, i+1)
    
    # Get data for this chromosome
    recomb_data = df_smoothed[df_smoothed['Chromosome'] == chrom].sort_values('Position')
    introgression_data = df_introgression_rate[df_introgression_rate['Chromosome'] == chrom].sort_values('Position')
    
    if len(recomb_data) == 0 or len(introgression_data) == 0:
        plt.text(0.5, 0.5, f"No data for {chrom}", ha='center', va='center')
        continue
    
    # Create primary y-axis for introgression rate (left axis)
    l1 = ax.plot(
        introgression_data['Position']/1e6, 
        introgression_data['IntrogressionRate'],
        color=introgression_color, 
        linewidth=2.5, 
        label='Introgression Rate'
    )
    
    # Fill area under the introgression curve
    ax.fill_between(
        introgression_data['Position']/1e6,
        introgression_data['IntrogressionRate'], 
        alpha=0.2, 
        color=introgression_color
    )
    
    # Set up secondary y-axis for recombination rate (right axis)
    ax2 = ax.twinx()
    
    # Plot recombination rate on secondary axis
    l2 = ax2.plot(
        recomb_data['Position']/1e6, 
        recomb_data['cMMb_smoothed'],
        color=recomb_color, 
        linewidth=2.5, 
        label='Recombination Rate (cM/Mb)'
    )
    
    # Fill area under the recombination curve
    ax2.fill_between(
        recomb_data['Position']/1e6,
        recomb_data['cMMb_smoothed'], 
        alpha=0.15, 
        color=recomb_color
    )
    
    # Set labels and title
    ax.set_title(f"{chrom}", fontsize=14)
    ax.set_xlabel('Position (Mb)', fontsize=11)
    ax.set_ylabel('Introgression Rate', color=introgression_color, fontsize=11)
    ax2.set_ylabel('Recombination Rate (cM/Mb)', color=recomb_color, fontsize=11)
    
    # Set tick colors to match data lines
    ax.tick_params(axis='y', colors=introgression_color)
    ax2.tick_params(axis='y', colors=recomb_color)
    
    # Set reasonable y-axis limits for recombination rate
    max_cmmb = min(15, np.percentile(recomb_data['cMMb_smoothed'], 99) * 1.2)
    ax2.set_ylim(0, max_cmmb)
    
    # Introgression rate is already 0-1
    ax.set_ylim(0, 1.05)
    
    # Calculate correlation between metrics at matched positions
    from scipy.interpolate import interp1d
    
    # Only calculate if sufficient data
    if len(introgression_data) > 3:
        # Create interpolation function for introgression data
        introgression_interp = interp1d(
            introgression_data['Position'],
            introgression_data['IntrogressionRate'],
            bounds_error=False, 
            fill_value='extrapolate'
        )
        
        # Get introgression values at recombination positions
        introgression_at_recomb = introgression_interp(recomb_data['Position'])
        
        # Calculate correlation
        if len(recomb_data) > 5:
            r, p = pearsonr(recomb_data['cMMb_smoothed'], introgression_at_recomb)
            rho, p_rho = spearmanr(recomb_data['cMMb_smoothed'], introgression_at_recomb)
            
            # Add correlation to plot
            ax.text(
                0.03, 0.95, 
                f"r = {r:.2f} (p = {p:.3f})\nρ = {rho:.2f} (p = {p_rho:.3f})", 
                transform=ax.transAxes, 
                fontsize=9,
                bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3')
            )
            
            correlations.append({
                'Chromosome': chrom,
                'Pearson_r': r,
                'Pearson_p': p,
                'Spearman_rho': rho,
                'Spearman_p': p_rho
            })
    
    # Add grid for readability
    ax.grid(True, linestyle='--', alpha=0.3)
    
    # Create combined legend
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    if i == 0:  # Only add legend to first plot
        ax.legend(
            lines, labels, 
            loc='upper right', 
            frameon=True, 
            framealpha=0.9, 
            fontsize=9
        )

# Add overall title
plt.suptitle('Introgression Rate and Recombination Rate Profiles by Chromosome', fontsize=16, y=0.98)

# Adjust layout
plt.tight_layout(rect=[0, 0, 1, 0.97])

# Save the figure
profile_output_path = os.path.join(output_dir, 'chromosome_profiles.png')
plt.savefig(profile_output_path, dpi=300, bbox_inches='tight')
print(f"Saved chromosome profiles to: {profile_output_path}")

# Save correlation results
if correlations:
    corr_df = pd.DataFrame(correlations)
    corr_path = os.path.join(output_dir, "chromosome_correlations.tsv")
    corr_df.to_csv(corr_path, sep='\t', index=False)
    print(f"Saved correlation results to: {corr_path}")

print("Analysis complete!")