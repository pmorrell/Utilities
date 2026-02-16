 #!/usr/bin/env python3		

import argparse
import gzip
from multiprocessing import Pool

def is_gzipped(filepath):
    """Check if a file is gzipped."""
    with open(filepath, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def extract_variants_from_vcf(vcf_file):
    """Yield variants from VCF file without loading everything into memory."""
    gzipped = is_gzipped(vcf_file)
    opener = gzip.open if gzipped else open
    mode = 'rt' if gzipped else 'r'

    with opener(vcf_file, mode) as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue

            yield {
                'chrom': parts[0],
                'pos': parts[1],
                'id': parts[2] if parts[2] != '.' else f"{parts[0]}_{parts[1]}",
                'ref': parts[3],
                'alt': parts[4].split(',')[0]  # Just take first alt allele
            }

def process_fasta_entry(entry):
    """Process each FASTA entry for parallel parsing."""
    lines = entry.strip().split("\n")
    header = lines[0]
    seq = "".join(lines[1:])
    # Ensure we get the full header ID without accidentally dropping characters
    seq_id = header.split()[0]
    return seq_id, seq

def parse_fasta_parallel(fasta_file):
    """Use multiprocessing for faster FASTA parsing."""
    gzipped = is_gzipped(fasta_file)
    opener = gzip.open if gzipped else open
    mode = 'rt' if gzipped else 'r'

    with opener(fasta_file, mode) as f:
        content = f.read()
        entries = [e for e in content.split(">") if e.strip()]  # Split on FASTA headers

    with Pool() as pool:
        sequences = dict(pool.map(process_fasta_entry, entries))

    # Debug output to verify IDs
    sample_ids = list(sequences.keys())[:5]
    print(f"Sample sequence IDs: {sample_ids}")

    return sequences

def create_illumina_format(fasta_file, vcf_file, position=61, output_file="illumina_probes.txt", flank_length=60):
    """Generate Illumina-style SNP assay design file."""
    # Extract and build variant dictionary
    variants = list(extract_variants_from_vcf(vcf_file))
    if not variants:
        print("Error: No variants found in VCF file.")
        return

    # Show sample variant formats
    print("Sample variants:")
    for v in variants[:3]:
        print(f"  Chrom: '{v['chrom']}', Pos: {v['pos']}")

    # Build variant dictionary
    variant_dict = {f"{v['chrom']}_{v['pos']}": v for v in variants}
    print(f"Found {len(variant_dict)} variants.")

    # Process FASTA sequences in parallel
    sequences = parse_fasta_parallel(fasta_file)
    print(f"Found {len(sequences)} sequences in FASTA file.")

    # Open output file
    warnings = []
    sequence_count = 0
    processed_variants = set()

    with open(output_file, 'w') as out:
        for seq_id, seq in sequences.items():
            # Debug: check what the sequence ID looks like
            if sequence_count == 0:
                print(f"First sequence ID being processed: '{seq_id}'")

            variant_key = seq_id
            if variant_key in variant_dict:
                variant = variant_dict[variant_key]
                pos_index = position - 1  # Convert to 0-based index

                # Ensure sequence length is sufficient
                if len(seq) <= pos_index or pos_index < flank_length or len(seq) - pos_index <= flank_length:
                    warnings.append(f"Warning: Sequence {seq_id} is too short for flanking extraction.")
                    continue

                # Extract flanking regions
                upstream = seq[pos_index-flank_length:pos_index]
                downstream = seq[pos_index+1:pos_index+flank_length+1]

                illumina_string = f"{seq_id}\t{upstream}[{variant['ref']}/{variant['alt']}]{downstream}"
                out.write(illumina_string + "\n")

                sequence_count += 1
                processed_variants.add(variant_key)
            else:
                # Check if we have a chromosome name mismatch
                parts = seq_id.split('_')
                if len(parts) == 2:
                    chrom, pos = parts
                    # Try different variations of the chromosome name
                    alternatives = []
                    if chrom.startswith('hr') and len(chrom) > 2:
                        alternatives.append(f"c{chrom}_{pos}")
                    if not chrom.startswith('chr'):
                        alternatives.append(f"chr{chrom}_{pos}")

                    # Try each alternative
                    for alt_key in alternatives:
                        if alt_key in variant_dict:
                            variant = variant_dict[alt_key]
                            # Process this variant with the fixed key
                            pos_index = position - 1
                            if len(seq) <= pos_index or pos_index < flank_length or len(seq) - pos_index <= flank_length:
                                warnings.append(f"Warning: Sequence {seq_id} is too short for flanking extraction.")
                                continue

                            upstream = seq[pos_index-flank_length:pos_index]
                            downstream = seq[pos_index+1:pos_index+flank_length+1]

                            illumina_string = f"{seq_id}\t{upstream}[{variant['ref']}/{variant['alt']}]{downstream}"
                            out.write(illumina_string + "\n")

                            sequence_count += 1
                            processed_variants.add(alt_key)
                            print(f"Fixed ID: {seq_id} → {alt_key}")
                            break
                    else:
                        warnings.append(f"Warning: No matching variant found for sequence {seq_id}")
                else:
                    warnings.append(f"Warning: No matching variant found for sequence {seq_id}")

    print(f"Processed {sequence_count}/{len(variant_dict)} variants.")

    # Report missing variants
    missing_variants = set(variant_dict.keys()) - processed_variants
    if missing_variants:
        warnings.append(f"{len(missing_variants)} variants did not have matching sequences.")
        # Show a few examples
        examples = sorted(list(missing_variants))[:5]
        warnings.append(f"Examples of missing variants: {', '.join(examples)}")

    # Print first few warnings and count
    if warnings:
        for warning in warnings[:20]:
            print(warning)
        if len(warnings) > 20:
            print(f"...and {len(warnings)-20} more warnings.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create Illumina SNP assay design format.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file")
    parser.add_argument("-p", "--position", type=int, default=61, help="Position to modify (1-based)")
    parser.add_argument("-o", "--output", default="illumina_probes.txt", help="Output file")
    parser.add_argument("--flank", type=int, default=60, help="Length of flanking sequence")

    args = parser.parse_args()
    create_illumina_format(args.fasta, args.vcf, args.position, args.output, args.flank)
    print(f"Illumina format file created: {args.output}")

