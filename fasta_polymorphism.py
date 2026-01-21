#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from cyvcf2 import VCF

def modify_sequences(fasta_file, vcf_file, position=61, output_ref="reference_sequences.fa", output_alt="alternate_sequences.fa"):
    """
    Modify the specified position in each FASTA sequence with VCF alleles
    """
    # Adjust position for 0-based indexing
    pos_index = position - 1

    # Parse VCF to get reference and alternate alleles
    vcf_reader = VCF(vcf_file)
    alleles = {}

    for record in vcf_reader:
        # Store reference and alternate alleles
        alleles[record.POS] = {
            'ref': record.REF,
            'alt': record.ALT
        }

    # Find variant closest to target position
    target_variant = min(alleles.keys(), key=lambda x: abs(x - position))
    ref_allele = alleles[target_variant]['ref']
    alt_alleles = alleles[target_variant]['alt']

    # Process FASTA sequences
    with open(output_ref, 'w') as ref_out, open(output_alt, 'w') as alt_out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq)

            if len(seq) <= pos_index:
                print(f"Warning: Sequence {record.id} is too short ({len(seq)}), skipping.")
                continue

            # Create reference version
            ref_seq = seq[:pos_index] + ref_allele + seq[pos_index+1:]
            ref_record = record
            ref_record.seq = Seq(ref_seq)
            ref_record.id = f"{record.id}_ref"

            # Write reference sequence
            SeqIO.write(ref_record, ref_out, "fasta")

            # Create alternate versions
            for i, alt in enumerate(alt_alleles):
                alt_seq = seq[:pos_index] + alt + seq[pos_index+1:]
                alt_record = record
                alt_record.seq = Seq(alt_seq)
                alt_record.id = f"{record.id}_alt{i+1}"

                # Write alternate sequence
                SeqIO.write(alt_record, alt_out, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Replace nucleotide at specific position with VCF variants.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file")
    parser.add_argument("-p", "--position", type=int, default=61, help="Position to modify (1-based, default=61)")
    parser.add_argument("--ref-out", default="reference_sequences.fa", help="Output file for reference sequences")
    parser.add_argument("--alt-out", default="alternate_sequences.fa", help="Output file for alternate sequences")

    args = parser.parse_args()
    modify_sequences(args.fasta, args.vcf, args.position, args.ref_out, args.alt_out)
    print(f"Created reference sequences in {args.ref_out}")
    print(f"Created alternate sequences in {args.alt_out}")
