#!/usr/bin/env python3
"""
Writes the working copy of the panel intervals BED used by every downstream module:
  - Harmonizes the chromosome naming with the reference fasta ("chr1" vs "1")
    when a mismatch was detected (Snakefile.smk).
  - Drops intervals on contigs absent from the reference (logged).
  - Drops track/browser/comment header lines (not accepted by every tool).
  - Guarantees a 4th (name) column, required by VarDict (-g 4). Missing names
    are filled with "chrom:start-end".
"""
import argparse
import sys
from collections import Counter


def main(bed, fai, output, needs_harmonization, ref_has_chr, log=sys.stderr):
    # Contigs actually present in the reference
    with open(fai) as f:
        valid_contigs = {line.split("\t")[0] for line in f}

    kept_lines = []
    dropped_contigs = Counter()
    n_headers = 0
    n_named = 0

    with open(bed) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith(("track", "browser", "#")):
                n_headers += 1
                continue

            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]

            # Only rewrite the naming if a mismatch was detected upfront
            if needs_harmonization:
                if ref_has_chr and not chrom.startswith("chr"):
                    chrom = f"chr{chrom}"
                elif not ref_has_chr and chrom.startswith("chr"):
                    chrom = chrom[3:]

            if chrom not in valid_contigs:
                dropped_contigs[chrom] += 1
                continue

            fields[0] = chrom
            if len(fields) < 4 or not fields[3].strip():
                fields = fields[:3] + [f"{chrom}:{fields[1]}-{fields[2]}"] + fields[4:]
                n_named += 1
            kept_lines.append("\t".join(fields) + "\n")

    if not kept_lines:
        raise ValueError(f"No intervals left in {bed} after harmonization: check the reference build.")

    with open(output, "w") as out:
        out.writelines(kept_lines)

    log.write(f"Chromosome naming harmonization applied: {needs_harmonization}\n")
    log.write(f"Intervals written: {len(kept_lines)}\n")
    log.write(f"Header lines removed: {n_headers}\n")
    log.write(f"Intervals without a name column (filled with coordinates): {n_named}\n")
    if dropped_contigs:
        total = sum(dropped_contigs.values())
        log.write(
            f"Dropped {total} interval(s) on {len(dropped_contigs)} contig(s) "
            f"absent from the reference fasta (likely alt-haplotype/random/"
            f"unplaced scaffolds not included in this assembly):\n"
        )
        for contig, n in sorted(dropped_contigs.items()):
            log.write(f"  {contig}: {n} interval(s)\n")
        log.write(
            "\nReview this list: if any of these contigs were expected to be "
            "covered by the panel, the reference fasta may be the wrong build.\n"
        )
    else:
        log.write("No intervals dropped; all contigs are present in the reference.\n")


if __name__ == "__main__":
    # Parameters from Snakemake
    if "snakemake" in locals() or "snakemake" in globals():
        with open(snakemake.log[0], "w") as log_handle:
            main(
                bed=snakemake.input.bed,
                fai=snakemake.input.fai,
                output=snakemake.output.bed,
                needs_harmonization=bool(snakemake.params.needs_harmonization),
                ref_has_chr=bool(snakemake.params.ref_has_chr),
                log=log_handle,
            )
    # Manual parameter parsing
    else:
        parser = argparse.ArgumentParser(description="Harmonize a panel BED with the reference fasta")
        parser.add_argument("--bed", required=True)
        parser.add_argument("--fai", required=True)
        parser.add_argument("--output", required=True)
        parser.add_argument("--ref_has_chr", action="store_true", help="Reference contigs use the 'chr' prefix")
        args = parser.parse_args()

        # Detect the BED naming from its first data line
        with open(args.bed) as f:
            first = next(l for l in f if l.strip() and not l.startswith(("track", "browser", "#")))
        bed_has_chr = first.startswith("chr")
        main(args.bed, args.fai, args.output, bed_has_chr != args.ref_has_chr, args.ref_has_chr)
