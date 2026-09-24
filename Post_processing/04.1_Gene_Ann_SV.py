#!/usr/bin/env python3
"""
Gene annotation of structural variants (Manta / Lumpy VCFs).

For every SV that passes the support (PE + SR) and FILTER criteria:
  - Both breakpoints are annotated with the protein-coding genes they hit:
    BND records use the mate position from the ALT field; DEL/DUP/INV/INS
    use POS and END.
  - If the two breakpoints hit different genes, the SV is reported as a
    fusion candidate "GeneA::GeneB". This includes intrachromosomal events
    such as the interstitial deletions STIL::TAL1 or FIP1L1::PDGFRA and
    inversions like CBFB::MYH11 when they are called as DEL/INV/DUP.
  - Otherwise, one row per affected gene is reported: breakpoint genes plus,
    for DEL/DUP, the genes fully or partially contained in the event
    (e.g. focal IKZF1 or CDKN2A/B deletions). Events spanning more than
    max_spanned_genes genes are collapsed into a single row.

Note: the order of the genes in "GeneA::GeneB" follows the VCF coordinates,
not the 5'->3' orientation of the fusion transcript.
"""
import argparse
import itertools
import re
import sys

import pandas as pd
import pyranges as pr
import pysam

INTERGENIC = "Intergenic"
DOSAGE_TYPES = {"DEL", "DUP", "CNV", "DUP:TANDEM"}
OUTPUT_COLUMNS = ["Gene", "SVType", "Support", "Coordinates", "SVLen",
                  "Fusion_candidate", "Genes_spanned", "Filter", "VAR_ID"]


def parse_bnd_alt(alt_string):
    """Extract mate chrom and pos from a BND ALT string (e.g. N]9:20441690])."""
    match = re.search(r'[\[\]]([^\[\]:]+):(\d+)[\[\]]', alt_string)
    if match:
        return match.group(1), int(match.group(2))
    return None, None


def strip_chr(chrom):
    return re.sub(r"^chr", "", str(chrom), flags=re.IGNORECASE)


def load_genes(gtf_file):
    """Protein-coding genes from the GTF (all genes if no biotype column)."""
    gtf = pr.read_gtf(gtf_file).df
    genes = gtf[gtf["Feature"] == "gene"]
    for biotype_col in ("gene_biotype", "gene_type"):
        if biotype_col in genes.columns:
            genes = genes[genes[biotype_col] == "protein_coding"]
            break
    genes = genes[["Chromosome", "Start", "End", "gene_name"]].dropna().copy()
    genes["Chromosome"] = genes["Chromosome"].astype(str).map(strip_chr)
    return pr.PyRanges(genes)


def overlapping_genes(intervals, genes):
    """intervals: list of dicts with Chromosome/Start/End/key.
    Returns {key: sorted list of overlapping gene names}."""
    if not intervals:
        return {}
    query = pr.PyRanges(pd.DataFrame(intervals))
    joined = query.join(genes).df
    if joined.empty:
        return {}
    joined["Chromosome"] = joined["Chromosome"].astype(str)
    return joined.groupby("key")["gene_name"].apply(lambda s: sorted(set(s))).to_dict()


def get_support(rec):
    """PE + SR alt support for Manta (FORMAT PR/SR) or Lumpy (INFO PE/SR)."""
    pe = sr = 0
    sample = rec.samples[0] if len(rec.samples) else {}
    # Manta: FORMAT PR/SR = (ref, alt)
    if "PR" in sample and sample["PR"] is not None:
        vals = sample["PR"]
        pe = vals[1] if isinstance(vals, (tuple, list)) and vals[1] is not None else 0
    if "SR" in sample and sample["SR"] is not None:
        vals = sample["SR"]
        sr = vals[1] if isinstance(vals, (tuple, list)) else vals
        sr = sr or 0
    # Lumpy: INFO PE/SR
    if "PE" in rec.info:
        val = rec.info["PE"]
        pe = sum(val) if isinstance(val, (tuple, list)) else val
    if "SR" in rec.info:
        val = rec.info["SR"]
        sr = sum(val) if isinstance(val, (tuple, list)) else val
    return int(pe or 0) + int(sr or 0)


def main(vcf_file, gtf_file, min_support, output, pass_only=True, max_spanned_genes=50, log=sys.stderr):
    print("Loading GTF...", file=log)
    genes = load_genes(gtf_file)

    # -------------------------------
    # 1. Read and filter the VCF
    # -------------------------------
    print("Reading VCF...", file=log)
    svs = []
    seen_mate_ids = set()
    n_total = n_filter = n_support = 0

    for rec in pysam.VariantFile(vcf_file):
        # For BND pairs keep only ONE representative record
        if rec.id in seen_mate_ids:
            continue
        n_total += 1
        svtype = str(rec.info.get("SVTYPE", "NA"))
        if svtype == "BND":
            mate_id = rec.info.get("MATEID")
            if mate_id:
                mate_id = mate_id[0] if isinstance(mate_id, (tuple, list)) else mate_id
                seen_mate_ids.add(mate_id)
            elif rec.id and rec.id.endswith("_1"):
                # Fallback for callers without MATEID (e.g. Lumpy "_1"/"_2" ids)
                seen_mate_ids.add(rec.id[:-1] + "2")

        filters = list(rec.filter.keys())
        filter_str = ";".join(filters) if filters else "."
        if pass_only and filters and filters != ["PASS"]:
            n_filter += 1
            continue

        support = get_support(rec)
        if support < min_support:
            n_support += 1
            continue

        chrom_a, pos_a = strip_chr(rec.chrom), rec.pos
        if svtype == "BND":
            chrom_b, pos_b = parse_bnd_alt(str(rec.alts[0]))
            if chrom_b is None:
                continue
            chrom_b = strip_chr(chrom_b)
            coord = f"{chrom_a}:{pos_a}-{chrom_b}:{pos_b}"
            svlen = ""
        else:
            chrom_b, pos_b = chrom_a, rec.stop
            coord = f"{chrom_a}:{pos_a}-{pos_b}"
            svlen = rec.info.get("SVLEN", "")
            svlen = svlen[0] if isinstance(svlen, (tuple, list)) else svlen

        svs.append({
            "VAR_ID": rec.id or f"sv{n_total}", "SVType": svtype, "Support": support,
            "Coordinates": coord, "SVLen": svlen, "Filter": filter_str,
            "chrom_a": chrom_a, "pos_a": pos_a, "chrom_b": chrom_b, "pos_b": pos_b,
        })

    print(f"SV records: {n_total} | removed by FILTER: {n_filter} | "
          f"removed by support < {min_support}: {n_support} | kept: {len(svs)}", file=log)

    # -------------------------------
    # 2. Annotate breakpoints and spanned regions
    # -------------------------------
    bp_intervals, span_intervals = [], []
    for i, sv in enumerate(svs):
        bp_intervals.append({"Chromosome": sv["chrom_a"], "Start": sv["pos_a"] - 1, "End": sv["pos_a"], "key": f"{i}_A"})
        bp_intervals.append({"Chromosome": sv["chrom_b"], "Start": sv["pos_b"] - 1, "End": sv["pos_b"], "key": f"{i}_B"})
        if sv["SVType"] in DOSAGE_TYPES and sv["pos_b"] > sv["pos_a"]:
            span_intervals.append({"Chromosome": sv["chrom_a"], "Start": sv["pos_a"] - 1, "End": sv["pos_b"], "key": str(i)})

    bp_genes = overlapping_genes(bp_intervals, genes)
    span_genes = overlapping_genes(span_intervals, genes)

    # -------------------------------
    # 3. Build the output rows
    # -------------------------------
    rows = []
    for i, sv in enumerate(svs):
        genes_a = bp_genes.get(f"{i}_A", [])
        genes_b = bp_genes.get(f"{i}_B", [])
        spanned = [g for g in span_genes.get(str(i), []) if g not in genes_a and g not in genes_b]
        base = {k: sv[k] for k in ("SVType", "Support", "Coordinates", "SVLen", "Filter", "VAR_ID")}
        base["Genes_spanned"] = ",".join(spanned) if len(spanned) <= max_spanned_genes \
            else f"{len(spanned)} genes ({','.join(spanned[:10])},...)"

        different_genes = genes_a and genes_b and set(genes_a).isdisjoint(genes_b)
        if sv["SVType"] == "BND" or different_genes:
            # Fusion candidate(s): one row per gene pair
            for gene_a, gene_b in itertools.product(genes_a or [INTERGENIC], genes_b or [INTERGENIC]):
                if gene_a == INTERGENIC and gene_b == INTERGENIC:
                    continue
                rows.append({**base, "Gene": f"{gene_a}::{gene_b}",
                             "Fusion_candidate": gene_a != INTERGENIC and gene_b != INTERGENIC})

        # Genes hit by a non-fusion event (same gene at both ends / one side only)
        if not (sv["SVType"] == "BND" or different_genes):
            for gene in sorted(set(genes_a) | set(genes_b)):
                rows.append({**base, "Gene": gene, "Fusion_candidate": False})

        # Genes contained in a deletion/duplication
        if spanned:
            if len(spanned) <= max_spanned_genes:
                for gene in spanned:
                    rows.append({**base, "Gene": gene, "Fusion_candidate": False})
            else:
                rows.append({**base, "Gene": "MULTIPLE_GENES", "Fusion_candidate": False})

    # -------------------------------
    # 4. Output (always written, even if empty)
    # -------------------------------
    df = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    df["SVLen"] = pd.to_numeric(df["SVLen"], errors="coerce").astype("Int64")
    if not df.empty:
        df = df.drop_duplicates().sort_values(by=["Fusion_candidate", "Support"], ascending=[False, False])
    df.to_csv(output, sep="\t", index=False)
    print(f"Table saved to {output}: {len(df)} rows from {df['VAR_ID'].nunique()} SVs.", file=log)


if __name__ == "__main__":
    # Parameters from Snakemake
    if "snakemake" in locals() or "snakemake" in globals():
        with open(snakemake.log[0], "w") as log_handle:
            main(
                vcf_file=snakemake.input.vcf,
                gtf_file=snakemake.input.gtf,
                min_support=int(snakemake.params.get("min_support", 3)),
                pass_only=bool(snakemake.params.get("pass_only", True)),
                max_spanned_genes=int(snakemake.params.get("max_spanned_genes", 50)),
                output=snakemake.output.tsv,
                log=log_handle,
            )
    # Manual parameter parsing
    else:
        parser = argparse.ArgumentParser(description="Annotate Manta/Lumpy SV VCFs with genes")
        parser.add_argument("--vcf_file", required=True)
        parser.add_argument("--gtf_file", required=True)
        parser.add_argument("--min_support", type=int, default=3)
        parser.add_argument("--keep_filtered", action="store_true", help="Keep non-PASS SV records")
        parser.add_argument("--max_spanned_genes", type=int, default=50)
        parser.add_argument("--output", required=True)
        args = parser.parse_args()

        main(args.vcf_file, args.gtf_file, args.min_support, args.output,
             pass_only=not args.keep_filtered, max_spanned_genes=args.max_spanned_genes)
