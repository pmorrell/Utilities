#!/usr/bin/env python3
#    Peter L. Morrell - Falcon Heights, MN - 20 February 2022
#    Find the coordinates of specific lists of genes in a GFF3 file.
#    Write out BED file that includes the coordinates of SNP intervals to include



import gzip
import gffutils

def read_list(s):
    gene_list = []
    with open(s, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                gene_list.append(line.strip())
    return gene_list


def read_to_database(d):
    gff3 = gffutils.example_filename(d)
    db = gffutils.create_db(gff3, dbfn='gff3.db', force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)
    return db

def main(genes, annotation):
    """Main function."""
    genes_retain = read_list(genes)
    gff3_db = read_to_database(annotation)
    for i in genes_retain:
        locus = gff3_db[i]
    print(locus.chrom, locus.start, locus.end, sep='\t')



def main(samples, ancestral):
    """Main function."""
    # Then iterate through the derived VCF and print out the relevant fields
    retain = read_list(samples)
    df = pd.read_csv(ancestral, compression='gzip', header=0, sep='\t')
    df_out = df[["Chromosome", "Pos", "SNPID", "Ancestral", "Derived", "Reference"] + retain]
    sample, variant_class = file_name(samples, ancestral)
    variant_class = file_name(ancestral)
    sample_file = sample + '_' + variant_class + 'anc.txt.gz'
    df_out.to_csv(sample_file, sep='\t', index=False, compression='gzip')


if len(sys.argv) <= 2:
    print("""Take two input files, a list of genes to extract and gzipped GFF3 file.
    1) List of genes to cut from GFF3
    2) A GFF3 """)
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
