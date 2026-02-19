#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
fst_data = pd.read_csv('/Users/pmorrell/Desktop/WBDC_Introgression/Fst/cult_intr.weir.fst', 
                        sep='\t', comment='//')

# Set the style
sns.set_theme(style="whitegrid")

# Create the figure and axes
plt.figure(figsize=(10, 6))

# Create the histogram
ax = sns.histplot(data=fst_data, x="WEIR_AND_COCKERHAM_FST", 
                 bins=50, kde=True, color="steelblue")

# Add vertical line at Fst = 0
plt.axvline(x=0, color='red', linestyle='--', alpha=0.7)

# Customize the plot
plt.title('Distribution of Weir & Cockerham Fst Values', fontsize=16)
plt.xlabel('Fst', fontsize=14)
plt.ylabel('Count', fontsize=14)

# Improve x-axis range - focus on the main distribution but show some extreme values
plt.xlim(-0.05, 1.0)

# Add text annotation showing mean and median
mean_fst = fst_data['WEIR_AND_COCKERHAM_FST'].mean()
median_fst = fst_data['WEIR_AND_COCKERHAM_FST'].median()
plt.text(0.7, 0.85, f'Mean: {mean_fst:.3f}\nMedian: {median_fst:.3f}', 
         transform=plt.gca().transAxes, 
         bbox=dict(facecolor='white', alpha=0.5))

# Save the figure
plt.tight_layout()
plt.savefig('/Users/pmorrell/Desktop/WBDC_Introgression/Fst/fst_histogram.png', dpi=300)

# Show the plot
plt.show()