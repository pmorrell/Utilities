import sys
import tskit
import tsinfer
import numpy as np
from Bio import AlignIO

def analyze_dhn7_trees(fasta_file):
    # 1. Dynamically figure out sample size and length
    try:
        alignment = AlignIO.read(fasta_file, "fasta")
    except Exception as e:
        return f"Error reading file: {e}"

    num_samples = len(alignment)
    alignment_len = alignment.get_alignment_length()
    
    print(f"File: {fasta_file}")
    print(f"Detected Samples (n): {num_samples}")
    print(f"Alignment Length (L): {alignment_len} bp")
    print("-" * 30)

    # 2. Prepare SampleData for tsinfer
    # We create a sample data object to hold our SNP data
    with tsinfer.SampleData(sequence_length=alignment_len) as sample_data:
        for pos in range(alignment_len):
            column = alignment[:, pos]
            unique_alleles = [a for a in set(column) if a != "-"]
            
            # Only process if the site is polymorphic (segregating)
            if len(unique_alleles) > 1:
                # Map nucleotides to integer IDs (0, 1, 2...)
                allele_map = {allele: i for i, allele in enumerate(unique_alleles)}
                # Default to -1 for gaps/missing data
                genotypes = [allele_map.get(base, -1) for base in column]
                
                # add_site requires genotypes for all samples
                if -1 not in genotypes:
                    sample_data.add_site(pos, genotypes, unique_alleles)

    # 3. Infer the Tree Sequence (the ARG)
    ts = tsinfer.infer(sample_data).simplify()

    # 4. Calculate Parameters
    # k+1 segments
    num_segments = ts.num_trees 
    
    # theta (diversity) based on branch lengths/pairwise diffs
    theta_locus = ts.diversity()
    theta_per_bp = theta_locus / alignment_len
    
    # rho (parametric estimate based on tree count)
    # rho = (K-1) / sum(1/i for i from 1 to n-1)
    harmonic_n = sum(1.0 / i for i in range(1, num_samples))
    rho_locus = (num_segments - 1) / harmonic_n
    rho_per_bp = rho_locus / alignment_len

    return {
        "segments": num_segments,
        "theta": theta_locus,
        "theta_bp": theta_per_bp,
        "rho": rho_locus,
        "rho_bp": rho_per_bp
    }

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py your_alignment.fasta")
    else:
        results = analyze_dhn7_trees(sys.argv[1])
        print(f"Tree Sequence Segments: {results['segments']}")
        print(f"Locus Theta (θ): {results['theta']:.4f} (per bp: {results['theta_bp']:.6f})")
        print(f"Locus Rho (ρ):   {results['rho']:.4f} (per bp: {results['rho_bp']:.6f})")
