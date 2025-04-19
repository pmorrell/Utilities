#!/usr/bin/env python3

import argparse
import gzip

def is_gzipped(filepath):
    """Check if a file is gzipped by looking at its first bytes"""
    with open(filepath, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def extract_variants_from_vcf(vcf_file):
    """Extract all variants from VCF file"""
    variants = []
    
    # Check if file is gzipped
    gzipped = vcf_file.endswith('.gz') or vcf_file.endswith('.bgz') or is_gzipped(vcf_file)
    
    # Open file with appropriate method
    opener = gzip.open if gzipped else open
    mode = 'rt' if gzipped else 'r'
    
    with opener(vcf_file, mode) as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                variants.append({
                    'chrom': parts[0],
                    'pos': parts[1],
                    'id': parts[2] if parts[2] != '.' else f"{parts[0]}_{parts[1]}",
                    'ref': parts[3],
                    'alt': parts[4].split(',')[0]  # Just take first alt allele for simplicity
                })
    
    return variants

def parse_fasta(fasta_file):
    """Simple FASTA parser that returns a dictionary of {sequence_id: sequence}"""
    sequences = {}
    current_id = None
    current_seq = []
    
    # Check if FASTA is gzipped
    gzipped = fasta_file.endswith('.gz') or is_gzipped(fasta_file)
    opener = gzip.open if gzipped else open
    mode = 'rt' if gzipped else 'r'
    
    with opener(fasta_file, mode) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Save the previous sequence if there was one
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                
                # Start a new sequence
                current_id = line[1:].split()[0]  # Extract ID
                current_seq = []
            else:
                current_seq.append(line)
    
    # Save the last sequence
    if current_id:
        sequences[current_id] = ''.join(current_seq)
        
    return sequences

def create_illumina_format(fasta_file, vcf_file, position=61, output_file="illumina_probes.txt", flank_length=60):
    """Create Illumina-style SNP assay design format"""
    # Get all variants from VCF
    variants = extract_variants_from_vcf(vcf_file)
    if not variants:
        print("Error: No variants found in VCF file")
        return
    
    print(f"Found {len(variants)} variants in VCF file")
    print(f"Reading FASTA file: {fasta_file}")
    
    # Process FASTA sequences
    sequences = parse_fasta(fasta_file)
    print(f"Found {len(sequences)} sequences in FASTA file")
    
    # Open output file
    with open(output_file, 'w') as out:
        sequence_count = 0
        
        # Iterate through variants and find matching sequences
        for i, variant in enumerate(variants):
            # For each variant, find FASTA sequences to match with it
            if i < len(sequences):
                # Get the i-th sequence (assumes ordered matching)
                seq_id = list(sequences.keys())[i]
                seq = sequences[seq_id]
                
                pos_index = position - 1  # Convert to 0-based
                
                if len(seq) <= pos_index or pos_index < flank_length or len(seq) - pos_index <= flank_length:
                    print(f"Warning: Sequence {seq_id} is too short for flanking extraction")
                    continue
                
                # Extract flanking sequences
                upstream = seq[pos_index-flank_length:pos_index]
                downstream = seq[pos_index+1:pos_index+flank_length+1]
                
                # Use the variant's ref/alt alleles
                ref_allele = variant['ref']
                alt_allele = variant['alt']
                
                # SNP name from VCF
                snp_name = f"{variant['chrom']}_{variant['pos']}"
                
                # Create the Illumina format string with actual variant from VCF
                illumina_string = f"{snp_name}\t{upstream}[{ref_allele}/{alt_allele}]{downstream}"
                out.write(illumina_string + "\n")
                
                sequence_count += 1
            
        print(f"Processed {sequence_count} sequences")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create Illumina SNP assay design format")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file")
    parser.add_argument("-p", "--position", type=int, default=61, help="Position to modify (1-based)")
    parser.add_argument("-o", "--output", default="illumina_probes.txt", help="Output file")
    parser.add_argument("--flank", type=int, default=60, help="Length of flanking sequence (default: 60)")
    
    args = parser.parse_args()
    create_illumina_format(args.fasta, args.vcf, args.position, args.output, args.flank)
    print(f"Illumina format file created: {args.output}")
