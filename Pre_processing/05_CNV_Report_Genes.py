#!/usr/bin/env python3
"""
Builds the list of genes shown in the gene-level sheet of the CNV report.

  - Default: protein-coding genes (GTF) overlapping the panel target intervals.
  - If a gene list is given (config cnv.report_genes, one symbol per line),
    only those genes are kept; symbols not found in the GTF are logged.

Output: tab-separated Chromosome, Start, End, Gene (chromosomes without "chr").
"""
import argparse
import sys

import pandas as pd
import pyranges as pr


def strip_chr(series):
    return series.astype(str).str.replace(r"^chr", "", case=False, regex=True)


def main(bed, gtf, output, genes_file=None, log=sys.stderr):
    genes = pr.read_gtf(gtf).df
    genes = genes[genes["Feature"] == "gene"]
    for col in ("gene_biotype", "gene_type"):
        if col in genes.columns:
            genes = genes[genes[col] == "protein_coding"]
            break
    genes = genes[["Chromosome", "Start", "End", "gene_name"]].dropna().rename(columns={"gene_name": "Gene"})
    genes["Chromosome"] = strip_chr(genes["Chromosome"])

    if genes_file:
        wanted = {l.strip() for l in open(genes_file) if l.strip() and not l.startswith("#")}
        selected = genes[genes["Gene"].isin(wanted)]
        missing = sorted(wanted - set(selected["Gene"]))
        if missing:
            log.write(f"Genes not found in the GTF (ignored): {', '.join(missing)}\n")
    else:
        targets = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2], names=["Chromosome", "Start", "End"])
        targets["Chromosome"] = strip_chr(targets["Chromosome"])
        hits = pr.PyRanges(genes).overlap(pr.PyRanges(targets)).df
        selected = genes[genes["Gene"].isin(set(hits["Gene"]))] if not hits.empty else genes.iloc[0:0]

    # One row per gene symbol (the longest locus if a symbol is duplicated)
    selected = (selected.assign(_len=selected["End"] - selected["Start"])
                .sort_values("_len", ascending=False)
                .drop_duplicates("Gene")
                .drop(columns="_len")
                .sort_values(["Chromosome", "Start"]))
    selected.to_csv(output, sep="\t", index=False)
    log.write(f"{len(selected)} report genes written to {output}\n")


if __name__ == "__main__":
    if "snakemake" in locals() or "snakemake" in globals():
        with open(snakemake.log[0], "w") as log_handle:
            genes_in = snakemake.input.get("genes")
            main(
                bed=snakemake.input.bed,
                gtf=snakemake.input.gtf,
                output=snakemake.output[0],
                genes_file=(genes_in[0] if isinstance(genes_in, list) and genes_in else genes_in) or None,
                log=log_handle,
            )
    else:
        parser = argparse.ArgumentParser(description="Genes for the gene-level CNV report")
        parser.add_argument("--bed", required=True)
        parser.add_argument("--gtf", required=True)
        parser.add_argument("--genes", help="Optional gene list (one symbol per line)")
        parser.add_argument("--output", required=True)
        args = parser.parse_args()
        main(args.bed, args.gtf, args.output, args.genes)
