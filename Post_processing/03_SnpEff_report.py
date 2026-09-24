#!/usr/bin/env python3
"""
Builds the per-sample small variant report from the SnpSift extractFields TSV.

Filters (in this order):
  1. Caller filter: only PASS (or unfiltered '.') calls, unless keep_filtered.
  2. Technical: alt-supporting reads >= min_depth and VAF >= min_vaf.
  3. Functional: drop synonymous and MODIFIER-impact variants unless they are
     splice-related.
  4. Classification: drop variants classified as (likely) benign.

The classification is a PRELIMINARY triage, not an ACMG/AMP call:
  - ClinVar significance when available (ClinVar is mostly germline-oriented).
  - Otherwise, the predicted SnpEff effect/impact.
It is meant to prioritise variants for expert review, not to replace it.
"""
import argparse
import re
import sys

import pandas as pd

# -----------------------------------------------------------------------------
# Classification
# -----------------------------------------------------------------------------
# ClinVar CLNSIG terms (lower case) -> report label. Order = priority when a
# record carries several terms (most informative first).
CLINVAR_LABELS = [
    ("pathogenic/likely_pathogenic", "Pathogenic/Likely pathogenic"),
    ("pathogenic", "Pathogenic"),
    ("likely_pathogenic", "Likely pathogenic"),
    ("uncertain_significance", "Uncertain significance"),
    ("benign/likely_benign", "Benign/Likely benign"),
    ("likely_benign", "Likely benign"),
    ("benign", "Benign"),
]
BENIGN_LABELS = {"Benign", "Likely benign", "Benign/Likely benign", "Likely benign (predicted)"}

# SnpEff effects considered splice-related (kept even if synonymous/MODIFIER)
SPLICE_RE = re.compile(r"splice", re.IGNORECASE)
# Inframe indels / duplications: key for FLT3-ITD, CALR, NPM1, CEBPA...
INFRAME_RE = re.compile(r"inframe_insertion|inframe_deletion|duplication", re.IGNORECASE)


def clinvar_label(clnsig):
    """Map a ClinVar CLNSIG value to a report label, or None if uninformative."""
    if clnsig is None or pd.isna(clnsig):
        return None
    value = str(clnsig).strip().lower()
    if value in ("", ".", "nan"):
        return None
    # Several values may be joined by ',' (alleles), '|' (secondary terms) or '&'
    terms = {t.strip() for t in re.split(r"[,|&;]", value) if t.strip()}
    # Conflicting interpretations/classifications -> uncertain
    if any(t.startswith("conflicting") for t in terms):
        return "Uncertain significance"
    for term, label in CLINVAR_LABELS:
        if term in terms:
            return label
    # drug_response, risk_factor, not_provided, other... carry no classification
    return None


def predicted_label(effect, impact):
    """Fallback label from the SnpEff prediction."""
    effect = str(effect).lower()
    impact = str(impact).upper()
    # HIGH: frameshift, stop_gained, splice_acceptor/donor, start/stop_lost, exon_loss...
    if impact == "HIGH":
        return "Predicted high impact"
    if SPLICE_RE.search(effect) or INFRAME_RE.search(effect):
        return "Uncertain significance"
    if impact == "MODERATE" or "missense" in effect:
        return "Uncertain significance"
    return "Likely benign (predicted)"


def classify(row):
    label = clinvar_label(row["CLNSIG"])
    if label is not None:
        return pd.Series({"Class": label, "Class_source": "ClinVar"})
    return pd.Series({
        "Class": predicted_label(row["ANN[0].EFFECT"], row["ANN[0].IMPACT"]),
        "Class_source": "Predicted effect",
    })


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
RENAME = {
    "CHROM": "Chrom",
    "POS": "Pos",
    "REF": "Ref",
    "ALT": "Alt",
    "ANN[0].GENE": "Gene",
    "ANN[0].FEATUREID": "Transcript",
    "ANN[0].RANK": "Exon",
    "ANN[0].HGVS_C": "HGVS_c",
    "ANN[0].HGVS_P": "HGVS_p",
    "ANN[0].EFFECT": "Effect",
    "ANN[0].IMPACT": "Impact",
    "GEN[0].DP": "Depth",
    "GEN[0].AD[1]": "Alt_depth",
    "GEN[0].AF": "VAF",
    "FILTER": "Filter",
    "CLNSIG": "ClinVar",
    "VARTYPE": "Var_type",
}
OUTPUT_COLUMNS = [
    "Chrom", "Pos", "Ref", "Alt", "Gene", "Transcript", "Exon", "HGVS_c", "HGVS_p",
    "Effect", "Impact", "Var_type", "Depth", "Alt_depth", "VAF", "Filter",
    "ClinVar", "Class", "Class_source",
]


def main(input, min_depth, output, min_vaf=0.0, keep_filtered=False, log=sys.stderr):
    df = pd.read_csv(input, sep="\t", dtype=str)
    n_in = len(df)

    # --- 1. Caller filter ---
    if not keep_filtered:
        df = df[df["FILTER"].fillna(".").isin(["PASS", "."])]
    n_pass = len(df)

    # --- 2. Technical filters ---
    for col in ("GEN[0].DP", "GEN[0].AD[1]", "GEN[0].AF"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[(df["GEN[0].AD[1]"] >= min_depth) & (df["GEN[0].AF"] >= min_vaf)]
    n_tech = len(df)

    # --- 3. Functional filters ---
    effect = df["ANN[0].EFFECT"].fillna("")
    impact = df["ANN[0].IMPACT"].fillna("").str.upper()
    is_splice = effect.str.contains(SPLICE_RE)
    is_synonymous = effect.str.contains("synonymous", case=False)
    df = df[~(is_synonymous & ~is_splice) & ((impact != "MODIFIER") | is_splice)]
    n_func = len(df)

    # --- 4. Classification ---
    if df.empty:
        df = df.assign(Class=pd.Series(dtype=str), Class_source=pd.Series(dtype=str))
    else:
        df = pd.concat([df, df.apply(classify, axis=1)], axis=1)
    df = df[~df["Class"].isin(BENIGN_LABELS)]

    # --- Formatting ---
    df["ANN[0].RANK"] = df["ANN[0].RANK"].astype(str).str.split("/").str[0].replace("nan", "")
    df["GEN[0].AF"] = df["GEN[0].AF"].round(4)
    df = df.rename(columns=RENAME)
    df = df[[c for c in OUTPUT_COLUMNS if c in df.columns]]
    df = df.sort_values(["Gene", "Pos"], kind="stable")
    df.to_csv(output, sep="\t", index=False)

    print(
        f"Input records: {n_in} | caller filter: {n_pass} | technical (alt reads >= {min_depth}, "
        f"VAF >= {min_vaf}): {n_tech} | functional: {n_func} | reported (non-benign): {len(df)}\n"
        f"Report written to: {output}",
        file=log,
    )


if __name__ == "__main__":
    # Parameters from Snakemake
    if "snakemake" in locals() or "snakemake" in globals():
        with open(snakemake.log[0], "w") as log_handle:
            main(
                input=snakemake.input.vars,
                min_depth=snakemake.params.get("min_depth", 8),
                min_vaf=float(snakemake.params.get("min_vaf", 0.0)),
                keep_filtered=bool(snakemake.params.get("keep_filtered", False)),
                output=snakemake.output.report,
                log=log_handle,
            )
    # Manual parameter parsing
    else:
        parser = argparse.ArgumentParser(description="Build the small variant report from a SnpSift TSV")
        parser.add_argument("--vars", required=True)
        parser.add_argument("--min_depth", type=int, default=8, help="Minimum alt-supporting reads")
        parser.add_argument("--min_vaf", type=float, default=0.0)
        parser.add_argument("--keep_filtered", action="store_true", help="Also report non-PASS calls")
        parser.add_argument("--report", required=True)
        args = parser.parse_args()

        main(args.vars, args.min_depth, args.report, args.min_vaf, args.keep_filtered)
