# Utilities
=========

2014-10-22

This repository contains various scripts and utilities that don't belong to any particular project.

Some of the code is quite useful *as is* while other code requires additional information to be useful. I'll try to update all pieces of code and example data files so that the code can be tested.


2026-02-16

| Script                   | Purpose                                                                 | Category                          |
|--------------------------|-------------------------------------------------------------------------|-----------------------------------|
| Ab1tofastq.py             | Use BioPython to convert Sanger sequence files to FASTQ reads           | Reuse of barley Sanger data       |
| Alchemy_Calls.sh          | Bash script to drive [Alchemy](https://alchemy.sourceforge.net) genotype calls | Barley genotyping datasets        |
| FASTQ_dumper.sh           | An SRA download script                                                  | Reuse of next generation sequence |
| file_rename.py            | Bulk renames files programmatically based on a pattern or mapping list. | general_utilities                 |
| list_unique.py            | Parses a file or column and outputs only the unique, non-duplicated entries. | general_utilities                 |
| getCol2.pl                | A quick Perl utility to extract the second column from tabular data.    | general_utilities                 |
| splitter.py | Splits large files to smaller pieces for S3 storage. |  general_utilities  |
| create_Illumina_lookup.py | Generates the lookup dictionaries needed for Illumina SNP microarrays.  | SNP_location                      |
| Illumina_lookup.py        | Maps or formats Illumina manifest and array probe data.                 | SNP_location                      |
| Illumina_lookup.sh        | Maps or formats Illumina manifest and array probe data.                 | SNP_location                      |
| GLnexus.sh                | Joint genotype calling across multiple gVCFs.                           | variant_calling                   |
| Alchemy_Calls.sh          | Drives genotype calling using the Alchemy software (useful for legacy array data). | variant_calling                   |
| paftools_call.sh          | Calls variants from PAF alignment formats using Minimap2's paftools.     | variant_calling                   |
| smoove.sh                 | Runs Smoove (which wraps Lumpy) to call and filter structural variants. | variant_calling                   |
| sniffles.sh | Identifies structural variants from long-read sequence alignments using Sniffles. | variant_calling |
| svim-asm.sh | Calls structural variants from whole-genome assembly alignments using SVIM-asm. |  variant_calling |

SNP_Utils_BLAST.sh - Uses BLAST to map or annotate the flanking sequences of SNPs.

SNP_Utils_Morex_v2.sh - Utility for handling SNPs specific to the barley Morex v2 assembly.

SNPMeta.sh - Assigns context or metadata annotations to SNP sets. | annotation

