#!/usr/bin/env python3

"""Convert a clean FASTA multiple sequence alignment into a VCF.

The first alignment sequence is treated as the aligned reference backbone.
If an external reference FASTA is supplied, its header may be in bedtools style,
for example:

    >chr2H:1000-2000
    >name::chr2H:1000-2000
    >name::chr2H:1000-2000(+)

This provides contig plus 1-based inclusive start/end coordinates.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO


@dataclass
class VariantRecord:
    chrom: str
    pos: int
    ref: str
    alts: list[str]
    genotypes: list[str]
    is_indel: bool = False


_REVCOMP_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def fasta_header_token(record) -> str:
    """Return the primary FASTA header token exactly as written (no spaces)."""
    return record.description.strip().split()[0]


def reverse_complement_sequence(seq: str) -> str:
    return seq.translate(_REVCOMP_TABLE)[::-1]


def reverse_complement_gapped(seq: str) -> str:
    # Gaps are preserved in-place after reversing because '-' is unchanged.
    return seq.translate(_REVCOMP_TABLE)[::-1]


def parse_header_coords(header: str) -> tuple[str, int, int] | None:
    """Parse coordinates from FASTA header.

    Supported examples:
    - chr2H:1000-2000
    - chr2H:1000_2000
    - name::chr2H:1000-2000
    - name::chr2H:1000-2000(+)
    - name_chr2H:1000-2000
    """
    token = header.split()[0]
    # bedtools may append strand annotation like (+) / (-)
    token = re.sub(r"\([+-]\)$", "", token)
    if "::" in token:
        token = token.rsplit("::", 1)[1]

    m = re.match(r"^(.+):(\d+)[-_](\d+)$", token)
    if not m:
        return None

    chrom = m.group(1)
    # Handle names where a sample/gene prefix is prepended before chromosome,
    # e.g. Vrn1_Morex_chr5H:528148015-528149277.
    if "_chr" in chrom.lower():
        chrom = chrom.rsplit("_", 1)[1]

    start = int(m.group(2))
    end = int(m.group(3))
    if start < 1 or end < start:
        raise ValueError("Invalid coordinates in reference header.")
    return chrom, start, end


def parse_alignment_fasta(path: Path) -> tuple[list[str], list[str], str]:
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError("Alignment FASTA is empty.")
    names = [fasta_header_token(r) for r in records]
    seqs = [str(r.seq).upper() for r in records]
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError("Alignment is not clean: sequence lengths differ.")
    return names, seqs, records[0].description


def parse_reference_fasta(path: Path) -> tuple[str, int, int, str]:
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError("Reference FASTA is empty.")
    if len(records) != 1:
        raise ValueError("Reference FASTA must contain exactly one sequence.")
    parsed = parse_header_coords(records[0].description)
    if parsed is None:
        raise ValueError(
            "Could not parse coordinates from reference FASTA header. "
            "Expected bedtools-style chr:start-end, optionally name:: prefix."
        )
    chrom, start, end = parsed
    seq = str(records[0].seq).upper().replace("-", "")
    return chrom, start, end, seq


def fallback_coords_from_name(header: str) -> tuple[str, int]:
    """Fallback coordinate parser for alignment reference name.

    Supports chr:start_end and legacy chr:start-end style.
    """
    parsed = parse_header_coords(header)
    if parsed is not None:
        chrom, start, _ = parsed
        return chrom, start
    token = header.split()[0]
    if "::" in token:
        token = token.split("::", 1)[0]
    return token, 1


def build_pos_map(ref_aln: str, start_1based: int) -> list[int | None]:
    out: list[int | None] = []
    offset = 0
    for base in ref_aln:
        if base == "-":
            out.append(None)
        else:
            out.append(start_1based + offset)
            offset += 1
    return out


def find_left_anchor(ref_aln: str, idx: int) -> int | None:
    j = idx - 1
    while j >= 0:
        if ref_aln[j] != "-":
            return j
        j -= 1
    return None


def is_variant_column(ref_base: str, col: list[str]) -> bool:
    for b in col:
        if b in {"N", "n"}:
            continue
        if b != ref_base:
            return True
    return False


def has_indel(ref_base: str, col: list[str]) -> bool:
    if ref_base == "-":
        return True
    return any(b == "-" for b in col)


def build_snp_record(chrom: str, i: int, ref_base: str, col: list[str], pos_map: list[int | None], sep: str = "|") -> VariantRecord | None:
    if ref_base == "-":
        return None

    ref = ref_base.upper()
    alleles = [ref]
    gts: list[str] = []
    for b in col:
        b = b.upper()
        if b in {"N", "-"}:
            gts.append(f".{sep}.")
            continue
        if b not in alleles:
            alleles.append(b)
        idx = alleles.index(b)
        gts.append(f"{idx}{sep}{idx}")

    if len(alleles) == 1:
        return None
    pos = pos_map[i]
    if pos is None:
        return None
    return VariantRecord(chrom=chrom, pos=pos, ref=ref, alts=alleles[1:], genotypes=gts)


def build_indel_record(
    chrom: str,
    block_start: int,
    block_end: int,
    ref_aln: str,
    aligned: list[str],
    pos_map: list[int | None],
    sep: str = "|",
) -> VariantRecord | None:
    anchor_idx = find_left_anchor(ref_aln, block_start)
    if anchor_idx is None:
        return None

    anchor = ref_aln[anchor_idx].upper()
    ref_mid = "".join(b.upper() for b in ref_aln[block_start : block_end + 1] if b != "-")
    ref_allele = anchor + ref_mid

    alleles = [ref_allele]
    gts: list[str] = []
    for seq in aligned:
        left = seq[anchor_idx].upper()
        mid_raw = seq[block_start : block_end + 1].upper()
        if left in {"N", "-"} or "N" in mid_raw:
            gts.append(f".{sep}.")
            continue

        mid = "".join(b for b in mid_raw if b != "-")
        allele = left + mid
        if allele not in alleles:
            alleles.append(allele)
        idx = alleles.index(allele)
        gts.append(f"{idx}{sep}{idx}")

    if len(alleles) == 1:
        return None
    pos = pos_map[anchor_idx]
    if pos is None:
        return None

    return VariantRecord(
        chrom=chrom,
        pos=pos,
        ref=ref_allele,
        alts=alleles[1:],
        genotypes=gts,
        is_indel=True,
    )


def build_records(chrom: str, start_1based: int, aligned: list[str], sep: str = "|") -> list[VariantRecord]:
    ref_aln = aligned[0]
    pos_map = build_pos_map(ref_aln, start_1based)

    out: list[VariantRecord] = []
    i = 0
    n = len(ref_aln)
    while i < n:
        col = [seq[i] for seq in aligned]
        ref_base = ref_aln[i]

        if not is_variant_column(ref_base, col):
            i += 1
            continue

        if has_indel(ref_base, col):
            j = i
            while j + 1 < n:
                nxt_col = [seq[j + 1] for seq in aligned]
                if has_indel(ref_aln[j + 1], nxt_col):
                    j += 1
                else:
                    break

            rec = build_indel_record(chrom, i, j, ref_aln, aligned, pos_map, sep)
            if rec is not None:
                out.append(rec)
            i = j + 1
            continue

        rec = build_snp_record(chrom, i, ref_base, col, pos_map, sep)
        if rec is not None:
            out.append(rec)
        i += 1

    return out


def write_vcf(
    path: Path,
    samples: list[str],
    records: list[VariantRecord],
    input_fastas: list[str],
    reverse_complemented_alignment: bool,
    unphased: bool = False,
) -> None:
    def build_info_field(rec: VariantRecord) -> str:
        called = [g for g in rec.genotypes if g not in {".|.", "./."}]
        ns = len(called)

        counts = [0] * (1 + len(rec.alts))
        for gt in called:
            sep = "/" if "/" in gt else "|"
            for allele in gt.split(sep):
                if allele == ".":
                    continue
                idx = int(allele)
                if 0 <= idx < len(counts):
                    counts[idx] += 1

        nonzero = [c for c in counts if c > 0]
        total_alleles = sum(counts)
        if len(nonzero) >= 2:
            mac = min(nonzero)
            maf = mac / float(total_alleles) if total_alleles else 0.0
        else:
            mac = 0
            maf = 0.0

        tag = "INDEL" if rec.is_indel else "SNP"
        return f"MAC={mac};MAF={maf:.12g};NS={ns};{tag}"

    with path.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.3\n")
        handle.write("##source=VCF_From_FASTA.py\n")
        for fasta_name in input_fastas:
            handle.write(f"##input_fasta={fasta_name}\n")
        handle.write(f"##reverse_complemented_alignment={'true' if reverse_complemented_alignment else 'false'}\n")
        handle.write("##INFO=<ID=MAC,Number=1,Type=Integer,Description=\"Minor allele count\">\n")
        handle.write("##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n")
        handle.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n")
        handle.write("##INFO=<ID=SNP,Number=0,Type=Flag,Description=\"Variant is a SNP\">\n")
        handle.write("##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Variant is an INDEL\">\n")
        phased_desc = "phased" if not unphased else "unphased"
        handle.write(
            f"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Diploid {phased_desc} genotype (homozygous for inbred lines)\">\n"
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for s in samples:
            handle.write(f"\t{s}")
        handle.write("\n")

        for rec in sorted(records, key=lambda r: (r.chrom, r.pos, r.ref, ",".join(r.alts))):
            row = [
                rec.chrom,
                str(rec.pos),
                ".",
                rec.ref,
                ",".join(rec.alts),
                ".",
                "PASS",
                build_info_field(rec),
                "GT",
            ]
            handle.write("\t".join(row + rec.genotypes) + "\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Convert aligned FASTA to VCF with SNPs and indels.")
    p.add_argument("alignment", help="Input aligned FASTA.")
    p.add_argument("-o", "--output", help="Output VCF path (default: <alignment stem>.vcf)")
    p.add_argument(
        "--reference-fasta",
        help=(
            "External reference FASTA using bedtools-style header, e.g. "
            "chr2H:1000-2000 or name::chr2H:1000-2000."
        ),
    )
    p.add_argument(
        "--reverse-complement-alignment",
        action="store_true",
        help="Reverse-complement all aligned sequences before variant calling.",
    )
    p.add_argument(
        "--unphased",
        action="store_true",
        help="Emit unphased genotypes (0/0, 1/1) instead of phased (0|0, 1|1).",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    aln_path = Path(args.alignment)

    sample_names, aligned, alignment_ref_header = parse_alignment_fasta(aln_path)
    ref_backbone = aligned[0].replace("-", "")
    input_fastas = [str(aln_path)]
    reverse_complemented_alignment = args.reverse_complement_alignment

    if args.reference_fasta:
        ref_path = Path(args.reference_fasta)
        chrom, start, end, ext_ref = parse_reference_fasta(ref_path)
        input_fastas.append(str(ref_path))
        expected_len = end - start + 1
        if len(ext_ref) != expected_len:
            raise ValueError(
                f"Reference FASTA length ({len(ext_ref)}) does not match header span ({expected_len})."
            )
        if ref_backbone != ext_ref:
            if reverse_complement_sequence(ref_backbone) == ext_ref:
                reverse_complemented_alignment = True
            else:
                raise ValueError("Ungapped alignment reference does not match provided reference FASTA sequence.")
    else:
        chrom, start = fallback_coords_from_name(alignment_ref_header)

    if reverse_complemented_alignment:
        aligned = [reverse_complement_gapped(seq) for seq in aligned]

    sep = "/" if args.unphased else "|"
    records = build_records(chrom=chrom, start_1based=start, aligned=aligned, sep=sep)

    out_path = Path(args.output) if args.output else aln_path.with_suffix("").with_name(f"{aln_path.stem}.vcf")
    write_vcf(out_path, sample_names, records, input_fastas, reverse_complemented_alignment, args.unphased)
    print(f"Wrote {len(records)} variants to {out_path}")


if __name__ == "__main__":
    main()
