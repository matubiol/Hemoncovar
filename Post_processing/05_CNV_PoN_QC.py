#!/usr/bin/env python3
"""
Quality control of the CNV panel of normals (leave-one-out).

Every normal has been normalized, segmented and called against a PoN built
without it. A good normal should look flat: low noise and no (or very few)
CNV calls. For each normal and caller this script reports:
  - noise: median absolute deviation (MAD) of the bin log2 ratios,
  - altered_segments / altered_Mb: gain/loss segments with >= min_probes bins
    (false positives by definition, since the sample is a healthy normal),
  - CNVkit: inferred sex and chrX/chrY log2 ratios (cnvkit.py sex).

Flags (column 'Flags'; 'Recommendation' = review if any flag):
  - sex_mismatch: declared sex (sample sheet) differs from the inferred sex.
  - high_noise: noise above median + 3 robust SD of all normals (per caller).
  - cnv_calls: >= 3 altered segments with >= min_probes bins (per caller).
  - possible_mosaic_LOY: male whose chrY log2 ratio is >= 0.5 below the
    median of the other males (mosaic loss of chrY is common in older men
    and distorts chrY in the reference).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def mad(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    return float(np.median(np.abs(values - np.median(values))))


def expected_cn(chrom, sex):
    chrom = str(chrom).upper().replace("CHR", "")
    if sex == "male" and chrom in ("X", "Y"):
        return 1
    if sex == "female" and chrom == "Y":
        return 0
    return 2


def cnvkit_metrics(cnr, call, sex, min_probes):
    bins = pd.read_csv(cnr, sep="\t")
    autosomes = ~bins["chromosome"].astype(str).str.upper().str.replace("CHR", "").isin(["X", "Y"])
    segs = pd.read_csv(call, sep="\t")
    if "probes" in segs.columns:
        segs = segs[segs["probes"] >= min_probes]
    exp = segs["chromosome"].map(lambda c: expected_cn(c, sex))
    if not sex:
        # Unknown sex: sex chromosomes cannot be evaluated
        exp = exp.where(~segs["chromosome"].astype(str).str.upper().str.replace("CHR", "").isin(["X", "Y"]))
    altered = segs[exp.notna() & (segs["cn"] != exp)]
    return {
        "cnvkit_noise_MAD": round(mad(bins.loc[autosomes, "log2"]), 4),
        "cnvkit_altered_segments": int(len(altered)),
        "cnvkit_altered_Mb": round(float((altered["end"] - altered["start"]).sum()) / 1e6, 2),
    }


def gatk_metrics(dcr, called, min_probes):
    bins = pd.read_csv(dcr, sep="\t", comment="@")
    segs = pd.read_csv(called, sep="\t", comment="@")
    segs = segs[segs["NUM_POINTS_COPY_RATIO"] >= min_probes]
    altered = segs[segs["CALL"].astype(str).isin(["+", "-"])]
    return {
        "gatk_noise_MAD": round(mad(bins["LOG2_COPY_RATIO"]), 4),
        "gatk_altered_segments": int(len(altered)),
        "gatk_altered_Mb": round(float((altered["END"] - altered["START"]).sum()) / 1e6, 2),
    }


def match_file(files, sample, suffix):
    for f in files:
        if Path(f).name == f"{sample}{suffix}":
            return f
    return None


def robust_outliers(series):
    values = series.astype(float)
    med = values.median()
    robust_sd = max(1.4826 * mad(values), 0.01)
    return values > med + 3 * robust_sd


def main(normals, sexes, output, min_probes=10, cnvkit_cnr=(), cnvkit_call=(), cnvkit_sex=None,
         gatk_dcr=(), gatk_called=(), log=sys.stderr):
    rows = []
    sex_table = None
    if cnvkit_sex:
        sex_table = pd.read_csv(cnvkit_sex, sep="\t")
        sex_table["sample"] = sex_table["sample"].map(lambda s: Path(str(s)).name.replace(".cnr", ""))

    for normal, sex in zip(normals, sexes):
        sex = (sex or "").lower()
        row = {"Normal": normal, "Declared_sex": sex or "unknown"}
        if sex_table is not None:
            s = sex_table[sex_table["sample"] == normal]
            if not s.empty:
                row["Inferred_sex"] = str(s["sex"].iloc[0]).lower()
                row["X_log2ratio"] = pd.to_numeric(s["X_logratio"].iloc[0], errors="coerce")
                row["Y_log2ratio"] = pd.to_numeric(s["Y_logratio"].iloc[0], errors="coerce")
        # Sex used to evaluate chrX/chrY: declared, otherwise inferred by CNVkit
        eval_sex = sex or (row.get("Inferred_sex") if row.get("Inferred_sex") in ("male", "female") else "")
        if cnvkit_cnr:
            row.update(cnvkit_metrics(match_file(cnvkit_cnr, normal, ".cnr"),
                                      match_file(cnvkit_call, normal, ".call.cns"), eval_sex, min_probes))
        if gatk_dcr:
            row.update(gatk_metrics(match_file(gatk_dcr, normal, ".denoisedCR.tsv"),
                                    match_file(gatk_called, normal, ".called.seg"), min_probes))
        rows.append(row)

    df = pd.DataFrame(rows)
    flags = {n: [] for n in df["Normal"]}

    if "Inferred_sex" in df.columns:
        for _, r in df.iterrows():
            if r["Declared_sex"] in ("male", "female") and r["Inferred_sex"] in ("male", "female") \
                    and r["Declared_sex"] != r["Inferred_sex"]:
                flags[r["Normal"]].append("sex_mismatch")
        males = df[(df["Inferred_sex"] == "male") & df["Y_log2ratio"].notna()]
        if len(males) >= 3:
            for _, r in males.iterrows():
                others = males.loc[males["Normal"] != r["Normal"], "Y_log2ratio"]
                if r["Y_log2ratio"] < others.median() - 0.5:
                    flags[r["Normal"]].append("possible_mosaic_LOY")

    for caller in ("cnvkit", "gatk"):
        noise_col, seg_col = f"{caller}_noise_MAD", f"{caller}_altered_segments"
        if noise_col in df.columns:
            for n in df.loc[robust_outliers(df[noise_col]), "Normal"]:
                flags[n].append(f"{caller}_high_noise")
            for n in df.loc[df[seg_col] >= 3, "Normal"]:
                flags[n].append(f"{caller}_cnv_calls")

    df["Flags"] = df["Normal"].map(lambda n: ",".join(flags[n]))
    df["Recommendation"] = np.where(df["Flags"] == "", "keep", "review")
    df.to_csv(output, sep="\t", index=False)
    log.write(f"PoN QC of {len(df)} normals written to {output}; "
              f"{(df['Recommendation'] == 'review').sum()} flagged for review.\n")


if __name__ == "__main__":
    if "snakemake" in locals() or "snakemake" in globals():
        with open(snakemake.log[0], "w") as log_handle:
            sex_in = list(snakemake.input.cnvkit_sex)
            main(
                normals=list(snakemake.params.normals),
                sexes=list(snakemake.params.sexes),
                output=snakemake.output[0],
                min_probes=int(snakemake.params.min_probes),
                cnvkit_cnr=list(snakemake.input.cnvkit_cnr),
                cnvkit_call=list(snakemake.input.cnvkit_call),
                cnvkit_sex=sex_in[0] if sex_in else None,
                gatk_dcr=list(snakemake.input.gatk_dcr),
                gatk_called=list(snakemake.input.gatk_called),
                log=log_handle,
            )
    else:
        parser = argparse.ArgumentParser(description="Leave-one-out QC of the CNV panel of normals")
        parser.add_argument("--normals", nargs="+", required=True)
        parser.add_argument("--sexes", nargs="+", required=True, help="male/female/unknown, same order as --normals")
        parser.add_argument("--cnvkit_cnr", nargs="*", default=[])
        parser.add_argument("--cnvkit_call", nargs="*", default=[])
        parser.add_argument("--cnvkit_sex")
        parser.add_argument("--gatk_dcr", nargs="*", default=[])
        parser.add_argument("--gatk_called", nargs="*", default=[])
        parser.add_argument("--min_probes", type=int, default=10)
        parser.add_argument("--output", required=True)
        a = parser.parse_args()
        sexes = ["" if s == "unknown" else s for s in a.sexes]
        main(a.normals, sexes, a.output, a.min_probes, a.cnvkit_cnr, a.cnvkit_call, a.cnvkit_sex,
             a.gatk_dcr, a.gatk_called)
