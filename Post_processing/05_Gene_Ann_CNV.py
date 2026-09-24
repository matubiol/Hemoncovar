#!/usr/bin/env python3
"""
Annotates CNV segments (GATK CallCopyRatioSegments or CNVkit call) with genes
and cytobands, and harmonizes both outputs into one table:
    Sample, Chromosome, Start, End, Num_Probes, Log2, [CN], Call, Gene, Cytoband, Arm

  - GATK: genes from the GTF (protein-coding); Call from the CALL column (+/-/0).
  - CNVkit: genes from the CNVkit 'gene' column (antitarget bins removed);
    Call from the absolute copy number, taking the sample sex into account
    for chrX/chrY (a single X copy in a male is not a loss).
"""
import argparse
import sys

import numpy as np
import pandas as pd
import pyranges as pr

GATK_CALLS = {"+": "Gain", "-": "Loss", "0": "Neutral"}
CNVKIT_NON_GENES = {"-", "Antitarget", "Background", ""}


def strip_chr(series):
    return series.astype(str).str.replace(r"^chr", "", case=False, regex=True)


def chrom_key(chrom):
    chrom = str(chrom).upper()
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if chrom in special:
        return special[chrom]
    try:
        return int(chrom)
    except ValueError:
        return 99


def expected_copies(chrom, sex):
    """Normal copy number of a chromosome given the sample sex."""
    chrom = str(chrom).upper()
    if sex == "male" and chrom in ("X", "Y"):
        return 1
    if sex == "female" and chrom == "Y":
        return 0
    return 2


def load_gatk(call, sample_name):
    df = pd.read_csv(call, sep="\t", comment="@")
    df = df.rename(columns={
        "CONTIG": "Chromosome", "START": "Start", "END": "End",
        "NUM_POINTS_COPY_RATIO": "Num_Probes", "MEAN_LOG2_COPY_RATIO": "Log2", "CALL": "Call",
    })
    df["Call"] = df["Call"].astype(str).map(GATK_CALLS).fillna(df["Call"].astype(str))
    df["Sample"] = sample_name or "Unknown"
    return df


def load_cnvkit(call, sample_name, sex):
    df = pd.read_csv(call, sep="\t")
    df = df.rename(columns={
        "chromosome": "Chromosome", "start": "Start", "end": "End",
        "probes": "Num_Probes", "log2": "Log2", "gene": "Gene", "cn": "CN",
    })
    df["Sample"] = sample_name or "Unknown"
    chroms = strip_chr(df["Chromosome"])
    expected = chroms.map(lambda c: expected_copies(c, sex))
    df["Call"] = np.where(df["CN"] > expected, "Gain", np.where(df["CN"] < expected, "Loss", "Neutral"))
    df["Gene"] = df["Gene"].apply(
        lambda x: ",".join(sorted({g for g in str(x).split(",") if g not in CNVKIT_NON_GENES}))
        if pd.notnull(x) else ""
    )
    return df


def genes_from_gtf(df, gtf):
    """Protein-coding genes overlapping each segment (by segment _id)."""
    genes = pr.read_gtf(gtf).df
    genes = genes[genes["Feature"] == "gene"]
    for col in ("gene_biotype", "gene_type"):
        if col in genes.columns:
            genes = genes[genes[col] == "protein_coding"]
            break
    genes = genes[["Chromosome", "Start", "End", "gene_name"]].dropna().copy()
    genes["Chromosome"] = strip_chr(genes["Chromosome"])
    query = df[["Chromosome", "Start", "End", "_id"]]
    joined = pr.PyRanges(query).join(pr.PyRanges(genes)).df
    if joined.empty:
        return pd.Series("", index=df["_id"])
    return joined.groupby("_id")["gene_name"].apply(lambda s: ",".join(sorted(set(s))))


def cytobands(df, cyto):
    cyto_df = pd.read_csv(cyto, sep="\t", header=None, names=["Chromosome", "Start", "End", "Band", "Stain"])
    cyto_df["Chromosome"] = strip_chr(cyto_df["Chromosome"])
    joined = pr.PyRanges(df[["Chromosome", "Start", "End", "_id"]]).join(pr.PyRanges(cyto_df)).df
    if joined.empty:
        return pd.DataFrame(columns=["_id", "Cytoband", "Arm"])
    joined = joined.sort_values(["_id", "Start_b"])

    def band_range(g):
        chrom = g["Chromosome"].iloc[0]
        first, last = g["Band"].iloc[0], g["Band"].iloc[-1]
        return f"{chrom}{first}" if first == last else f"{chrom}{first}-{last}"

    def arms(g):
        return "".join(sorted({str(b)[0] for b in g["Band"] if pd.notnull(b) and str(b)}))

    grouped = joined.groupby("_id")
    return pd.DataFrame({
        "Cytoband": grouped.apply(band_range, include_groups=False) if _pandas_has_include_groups() else grouped.apply(band_range),
        "Arm": grouped.apply(arms, include_groups=False) if _pandas_has_include_groups() else grouped.apply(arms),
    }).reset_index()


def _pandas_has_include_groups():
    major, minor = (int(x) for x in pd.__version__.split(".")[:2])
    return (major, minor) >= (2, 2)


def main(call, cyto, output, min_probes=0, format_type="gatk", sample_name=None, gtf=None, sex="", log=sys.stderr):
    # 1. Load and standardize input
    if format_type == "gatk":
        df = load_gatk(call, sample_name)
    elif format_type == "cnvkit":
        df = load_cnvkit(call, sample_name, sex)
    else:
        raise ValueError(f"Unknown format: {format_type}")

    df["Chromosome"] = strip_chr(df["Chromosome"])
    df = df[df["Num_Probes"] >= min_probes].reset_index(drop=True)
    df["_id"] = np.arange(len(df))

    if df.empty:
        cols = ["Sample", "Chromosome", "Start", "End", "Num_Probes", "Log2", "CN", "Call", "Gene", "Cytoband", "Arm"]
        pd.DataFrame(columns=[c for c in cols if c != "CN" or format_type == "cnvkit"]).to_csv(output, sep="\t", index=False)
        print("No segments left after filtering.", file=log)
        return

    # 2. Gene annotation
    if format_type == "gatk":
        if gtf:
            df["Gene"] = df["_id"].map(genes_from_gtf(df, gtf)).fillna("")
        else:
            df["Gene"] = ""

    # 3. Cytoband annotation
    df = df.merge(cytobands(df, cyto), on="_id", how="left")
    df[["Cytoband", "Arm"]] = df[["Cytoband", "Arm"]].fillna("")

    # 4. Sort and save
    df["_key"] = df["Chromosome"].map(chrom_key)
    df = df.sort_values(["Sample", "_key", "Start"]).drop(columns=["_key", "_id"])
    lead = [c for c in ["Sample", "Chromosome", "Start", "End", "Num_Probes", "Log2", "CN", "Call", "Gene", "Cytoband", "Arm"] if c in df.columns]
    df = df[lead + [c for c in df.columns if c not in lead]]
    df.to_csv(output, sep="\t", index=False)
    print(f"{len(df)} segments written to: {output} "
          f"(gains: {(df['Call'] == 'Gain').sum()}, losses: {(df['Call'] == 'Loss').sum()})", file=log)


if __name__ == "__main__":
    # Parameters from Snakemake
    if "snakemake" in locals() or "snakemake" in globals():
        with open(snakemake.log[0], "w") as log_handle:
            main(
                call=snakemake.input.call,
                cyto=snakemake.input.cyto,
                gtf=snakemake.input.get("gtf"),
                output=snakemake.output.tsv,
                min_probes=snakemake.params.get("min_probes", 0),
                format_type=snakemake.params.format,
                sample_name=snakemake.params.sample,
                sex=snakemake.params.get("sex", ""),
                log=log_handle,
            )
    # Manual parameter parsing
    else:
        parser = argparse.ArgumentParser(description="Annotate CNV calls (GATK or CNVkit) with genes and cytobands.")
        parser.add_argument("--call", required=True, help="Input CNV file")
        parser.add_argument("--format", choices=["gatk", "cnvkit"], required=True, help="Input format")
        parser.add_argument("--gtf", help="GTF for gene annotation (GATK format)")
        parser.add_argument("--cyto", required=True, help="Cytoband file")
        parser.add_argument("--min_probes", type=int, default=0)
        parser.add_argument("--sample", help="Sample name")
        parser.add_argument("--sex", default="", choices=["", "male", "female"])
        parser.add_argument("--output", required=True)
        args = parser.parse_args()

        main(args.call, args.cyto, args.output, args.min_probes, args.format, args.sample, args.gtf, args.sex)
