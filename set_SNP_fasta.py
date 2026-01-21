#!/usr/bin/env python3
import sys
from collections import defaultdict

def create_utils_lookup(utils_file):
    lookup = {}
    with open(utils_file) as f:
        for line_num, line in enumerate(f, 1):
            try:
                # Parse the line
                coords, sequence = line.strip().split('\t')
                chrom_range, seq = coords, sequence

                # Extract chromosome and position
                chrom, pos_range = chrom_range.split(':')
                start, end = map(int, pos_range.split('-'))
                mid_pos = (start + end) // 2

                # Find the variant position and alleles
                left_bracket = sequence.find('[')
                right_bracket = sequence.find(']')
                if left_bracket != -1 and right_bracket != -1:
                    alleles = sequence[left_bracket+1:right_bracket]
                    ref, alt = alleles.split('/')
                    key = f"{chrom}_{mid_pos}"
                    lookup[key] = (ref, alt)

            except Exception as e:
                print(f"Error on line {line_num}: {line.strip()}", file=sys.stderr)
                print(f"Error details: {str(e)}", file=sys.stderr)

    print(f"Loaded {len(lookup)} positions from utils file", file=sys.stderr)
    return lookup

def process_fasta(fasta_file, utils_lookup):
    missing_pos = defaultdict(int)
    processed = 0
    replaced = 0
    with open(fasta_file) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()

            # Process header to get middle position
            header = header[1:]  # Remove '>'
            chrom, range_part = header.split(':')
            start, end = map(int, range_part.split('-'))
            mid_pos = (start + end) // 2

            # Create lookup key
            lookup_key = f"{chrom}_{mid_pos}"

            # Get reference nucleotide from utils lookup
            n_pos = seq.find('N')
            if n_pos != -1:
                if lookup_key in utils_lookup:
                    ref_nuc = utils_lookup[lookup_key][0]  # Use reference nucleotide
                    seq = seq[:n_pos] + ref_nuc + seq[n_pos + 1:]
                    replaced += 1
                else:
                    missing_pos[chrom] += 1

            print(f">{chrom}_{mid_pos}")
            print(seq)
            processed += 1

    print(f"\nProcessed {processed} sequences", file=sys.stderr)
    print(f"Successfully replaced {replaced} N's with reference nucleotides", file=sys.stderr)
    print("\nMissing positions by chromosome:", file=sys.stderr)
    for chrom, count in sorted(missing_pos.items()):
        print(f"{chrom}: {count} positions not found", file=sys.stderr)

if __name__ == '__main__':
    utils_lookup = create_utils_lookup('SSW_snps_morex_v2_SNP_Utils.txt')
    process_fasta('SSW_snps_final_morex_v2_renamed.fasta', utils_lookup)
