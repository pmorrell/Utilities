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
    
    # Collect output lines
    output_lines = []
    
    # Iterate through variants and find matching sequences
    for i, variant in enumerate(variants):
        # Use the variant's chromosome directly without modification
        variant_chrom = variant['chrom']
        seq_id = variant_chrom
        
        # Check if sequence ID matches the variant chromosome
        if seq_id not in sequences:
            print(f"Warning: No matching sequence found for {variant_chrom} in FASTA")
            continue  # Skip if no matching sequence in FASTA file

        seq = sequences[seq_id]
        pos_index = position - 1  # Convert to 0-based
        
        if len(seq) <= pos_index or pos_index < flank_length or len(seq) - pos_index <= flank_length:
            print(f"Warning: Sequence {seq_id} is too short for flanking extraction")
            continue
        
        # Extract flanking sequences
        upstream = seq[pos_index-flank_length:pos_index]
        downstream = seq[pos_index+1:pos_index+flank_length+1]
        
        # SNP name from VCF
        snp_name = f"{variant['chrom']}_{variant['pos']}"
        
        # Create the Illumina format string with actual variant from VCF
        illumina_string = f"{snp_name}\t{upstream}[{variant['ref']}/{variant['alt']}] {downstream}"
        output_lines.append(illumina_string)
    
    # Write the output all at once
    with open(output_file, 'w') as out:
        out.write("\n".join(output_lines) + "\n")
    
    print(f"Illumina format file created: {output_file}")
    print(f"Processed {len(output_lines)} sequences")

