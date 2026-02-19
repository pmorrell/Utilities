#!/usr/bin/env python
"""Calculate the cM/Mb rate for the BOPA SNPs."""

import sys

try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
except ImportError:
    lowess = None

physical_vcf = sys.argv[1]
genetic_map = sys.argv[2]
lowess_frac = float(sys.argv[3]) if len(sys.argv) > 3 else 0.2

if lowess is None:
    raise ImportError("statsmodels is required for LOWESS smoothing. Install with 'pip install statsmodels'.")


gen = {}
with open(genetic_map, 'r') as f:
    for index, line in enumerate(f):
        if index == 0 and line.startswith('#'):
            continue
        else:
            # BED format: chrom, 0-based pos, 1-based pos, snp name, genetic map pos
            t = line.strip().split('\t')
            if len(t) < 5:
                continue
            gen[t[3]] = (t[0], t[4])

phys = {}
ordered_snps = {}
with open(physical_vcf, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            t = line.strip().split()
            #   Only get the physical positions of the SNPs that are on the
            #   genetic map.
            if t[2] not in gen:
                continue
            else:
                phys[t[2]] = (t[0], t[1])
                if t[0] not in ordered_snps:
                    ordered_snps[t[0]] = []
                ordered_snps[t[0]].append(t[2])

print('\t'.join(['Chromosome', 'LeftBP', 'RightBP', 'cMMb', 'cMMb_lowess']))
for chrom, markers in ordered_snps.items():
    markers_sorted = sorted(markers, key=lambda snp: int(phys[snp][1]))
    cmmb_values = []
    midpoints = []
    rows = []
    for i in range(0, len(markers_sorted)-1):
        j = i + 1
        # Skip if genetic map values are missing ('--')
        if gen[markers_sorted[j]][1] == '--' or gen[markers_sorted[i]][1] == '--':
            continue
        delta_phys = int(phys[markers_sorted[j]][1]) - int(phys[markers_sorted[i]][1])
        delta_gen = float(gen[markers_sorted[j]][1]) - float(gen[markers_sorted[i]][1])
        if delta_phys > 0:
            #cmmb = (delta_gen/delta_phys) * 1000000
            cmmb = (delta_gen/delta_phys) * 100000 #run 100 kb to see how the results look!!!
            # Ensure cM/Mb is never negative
            cmmb = max(0, cmmb)
        else:
            cmmb = 'NA'
        left_bp = phys[markers_sorted[i]][1]
        right_bp = phys[markers_sorted[j]][1]
        rows.append([chrom, left_bp, right_bp, cmmb])
        if cmmb != 'NA':
            cmmb_values.append(float(cmmb))
            midpoints.append((int(left_bp) + int(right_bp)) / 2.0)

    smoothed = {}
    if cmmb_values and len(cmmb_values) > 3:
        # Filter to only finite values for LOWESS
        import numpy as np
        valid_idx = [i for i, (v, m) in enumerate(zip(cmmb_values, midpoints)) if np.isfinite(v)]
        if len(valid_idx) > 3:
            valid_cmmb = [cmmb_values[i] for i in valid_idx]
            valid_midpoints = [midpoints[i] for i in valid_idx]
            lowess_result = lowess(valid_cmmb, valid_midpoints, frac=lowess_frac, return_sorted=True)
            for x, y in lowess_result:
                smoothed[x] = y

    for chrom, left_bp, right_bp, cmmb in rows:
        if cmmb == 'NA':
            cmmb_lowess = 'NA'
        else:
            midpoint = (int(left_bp) + int(right_bp)) / 2.0
            cmmb_lowess = smoothed.get(midpoint, 'NA')
        toprint = '\t'.join([chrom, left_bp, right_bp, str(cmmb), str(cmmb_lowess)])
        print(toprint)
