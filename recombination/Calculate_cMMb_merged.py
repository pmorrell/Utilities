#!/usr/bin/env python3
"""
Calculate cM/Mb recombination rates from a PLINK .map file.

This script computes recombination rates by comparing genetic (cM) and physical (bp) 
distances between consecutive markers. Optionally applies LOWESS smoothing to 
reduce noise in the resulting rate estimates.

Output format:
Tab-separated with columns: Chromosome, LeftBP, RightBP, cMMb, cMMb_lowess (if smoothing enabled)

References:
Modified from: https://github.com/MorrellLAB/Env_Assoc/blob/master/script/recombination/Calculate_cMMb.py
"""

import sys
import argparse

try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
    import numpy as np
except ImportError:
    lowess = None
    np = None


def calculate_cm_mb(map_file, output_file, apply_smoothing=False, lowess_frac=0.2):
    """
    Calculate cM/Mb rates from a PLINK .map file.
    
    PLINK .map format:
    Col 1: Chromosome
    Col 2: SNP ID
    Col 3: Genetic Distance (cM)
    Col 4: Physical Position (bp)
    
    Args:
        map_file: Path to input PLINK .map file
        output_file: Path to output results file
        apply_smoothing: Whether to apply LOWESS smoothing
        lowess_frac: Fraction of data to use for local regression (default: 0.2)
    """
    
    print(f"Reading from: {map_file}")
    print(f"Writing to:   {output_file}")

    with open(map_file, 'r') as fin:
        # Parse PLINK .map file
        chrom_data = {}  # chrom -> list of (bp, cm)
        
        for line in fin:
            parts = line.strip().split()
            if not parts:
                continue
            
            try:
                chrom = parts[0]
                cm = float(parts[2])
                bp = int(parts[3])
            except (ValueError, IndexError):
                # Skip lines that don't look like data (e.g., headers)
                continue
            
            if chrom not in chrom_data:
                chrom_data[chrom] = []
            chrom_data[chrom].append((bp, cm))
        
        # Sort by physical position
        for chrom in chrom_data:
            chrom_data[chrom].sort()
    
    # Calculate rates and optionally apply smoothing
    with open(output_file, 'w') as fout:
        # Write header
        if apply_smoothing:
            fout.write("Chromosome\tLeftBP\tRightBP\tcMMb\tcMMb_lowess\n")
        else:
            fout.write("Chromosome\tLeftBP\tRightBP\tcMMb\n")
        
        total_intervals = 0
        
        for chrom, markers in chrom_data.items():
            if len(markers) < 2:
                continue
            
            # Calculate cM/Mb for each interval
            cmmb_values = []
            midpoints = []
            rows = []
            
            for i in range(len(markers) - 1):
                left_bp, left_cm = markers[i]
                right_bp, right_cm = markers[i + 1]
                
                delta_bp = right_bp - left_bp
                delta_cm = right_cm - left_cm
                
                if delta_bp > 0:
                    # Convert bp to Mb (1 Mb = 1,000,000 bp)
                    delta_mb = delta_bp / 1_000_000.0
                    cmmb = delta_cm / delta_mb
                    # Ensure cM/Mb is never negative
                    cmmb = max(0, cmmb)
                else:
                    cmmb = 'NA'
                
                rows.append([chrom, left_bp, right_bp, cmmb])
                total_intervals += 1
                
                if cmmb != 'NA':
                    cmmb_values.append(float(cmmb))
                    midpoints.append((left_bp + right_bp) / 2.0)
            
            # Apply LOWESS smoothing if requested
            smoothed = {}
            if apply_smoothing and cmmb_values and len(cmmb_values) > 3:
                if lowess is None or np is None:
                    raise ImportError("statsmodels and numpy are required for LOWESS smoothing. "
                                    "Install with 'pip install statsmodels numpy'")
                
                # Filter to only finite values
                valid_idx = [i for i, v in enumerate(cmmb_values) if np.isfinite(v)]
                
                if len(valid_idx) > 3:
                    valid_cmmb = [cmmb_values[i] for i in valid_idx]
                    valid_midpoints = [midpoints[i] for i in valid_idx]
                    
                    lowess_result = lowess(valid_cmmb, valid_midpoints, frac=lowess_frac, return_sorted=True)
                    for x, y in lowess_result:
                        smoothed[x] = y
            
            # Write results
            for chrom, left_bp, right_bp, cmmb in rows:
                if apply_smoothing:
                    if cmmb == 'NA':
                        cmmb_lowess = 'NA'
                    else:
                        midpoint = (left_bp + right_bp) / 2.0
                        cmmb_lowess = smoothed.get(midpoint, 'NA')
                    fout.write(f"{chrom}\t{left_bp}\t{right_bp}\t{cmmb}\t{cmmb_lowess}\n")
                else:
                    fout.write(f"{chrom}\t{left_bp}\t{right_bp}\t{cmmb}\n")
    
    print(f"Done. Processed {total_intervals} intervals.")


def main():
    parser = argparse.ArgumentParser(
        description='Calculate cM/Mb recombination rates from a PLINK .map file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage without smoothing
  python Calculate_cMMb.py input.map output.txt
  
  # With LOWESS smoothing
  python Calculate_cMMb.py input.map output.txt --smooth
  
  # With custom LOWESS fraction
  python Calculate_cMMb.py input.map output.txt --smooth --lowess-frac 0.3
        """)
    
    parser.add_argument('map_file', help='Input PLINK .map file')
    parser.add_argument('output_file', help='Output file')
    parser.add_argument('--smooth', action='store_true', 
                        help='Apply LOWESS smoothing to cM/Mb values')
    parser.add_argument('--lowess-frac', type=float, default=0.2,
                        help='Fraction of data for LOWESS local regression (default: 0.2)')
    
    args = parser.parse_args()
    
    calculate_cm_mb(args.map_file, args.output_file, 
                    apply_smoothing=args.smooth, 
                    lowess_frac=args.lowess_frac)


if __name__ == "__main__":
    main()
