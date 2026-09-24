#!/usr/bin/env python3
"""
Per-sample CNV report (Excel) combining CNVkit and GATK.

Sheets:
  Summary            Cytogenetic-style summary per caller (e.g. "+8, -7, 5q-, 17p-")
                     and events found by both callers.
  Arm_level          Every chromosome arm: fraction of its segmented length that
                     is gained / lost, arm call per caller and concordance.
                     An arm is called when >= arm_min_fraction of its segmented
                     length is altered and it is covered by >= min_probes bins.
  Segments           Altered segments (>= min_probes bins) of each caller, with
                     report genes, cytoband and support by the other caller
                     (>= 50% reciprocal overlap, same direction).
  Genes              Report genes (panel genes by default) with the overlapping
                     segment of each caller.
  Allelic_imbalance  GATK ModelSegments segments with minor allele fraction
                     <= loh_max_maf (>= loh_min_hets het SNPs): copy-neutral LOH,
                     LOH with loss or imbalanced gain; CNVkit segments with
                     minor allele copy number 0 (cn2 == 0).
  Parameters         Settings and inputs used.

Notes:
  - GATK excludes chrX/chrY by default (mixed-sex PoN); sex chromosomes are
    then reported by CNVkit only (sex-aware when the sample sex is known).
  - Tumour-only allelic imbalance: in high-purity samples, heterozygous SNPs
    inside LOH regions look homozygous and may be dropped, so large LOH regions
    can appear with few supporting SNPs. Review with the plots.
"""
import argparse
import sys

import numpy as np
import pandas as pd

CALLERS = ("CNVkit", "GATK")
ALTERED = ("Gain", "Loss")


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def strip_chr(series):
    return series.astype(str).str.replace(r"^chr", "", case=False, regex=True)


def chrom_key(chrom):
    chrom = str(chrom).upper()
    return {"X": 23, "Y": 24}.get(chrom, int(chrom) if chrom.isdigit() else 99)


def load_segments(path):
    """Annotated segments from 05_Gene_Ann_CNV.py (all segments, any call)."""
    if not path:
        return None
    df = pd.read_csv(path, sep="\t", dtype={"Chromosome": str, "Gene": str})
    df["Chromosome"] = strip_chr(df["Chromosome"])
    df["Gene"] = df["Gene"].fillna("")
    return df


def load_arms(cyto):
    """Chromosome arms from the cytoband table (centromere = 'acen' bands)."""
    bands = pd.read_csv(cyto, sep="\t", header=None, names=["Chromosome", "Start", "End", "Band", "Stain"])
    bands["Chromosome"] = strip_chr(bands["Chromosome"])
    arms = []
    for chrom, b in bands.groupby("Chromosome"):
        if chrom_key(chrom) == 99:
            continue
        acen = b[b["Stain"] == "acen"]
        end = int(b["End"].max())
        if acen.empty:
            arms.append((chrom, "p", 0, end))
            continue
        arms.append((chrom, "p", 0, int(acen["Start"].min())))
        arms.append((chrom, "q", int(acen["End"].max()), end))
    arms = pd.DataFrame(arms, columns=["Chromosome", "Arm", "Start", "End"])
    arms["_key"] = arms["Chromosome"].map(chrom_key)
    return arms.sort_values(["_key", "Arm"]).drop(columns="_key").reset_index(drop=True), bands


def band_range(bands, chrom, start, end):
    b = bands[(bands["Chromosome"] == chrom) & (bands["End"] > start) & (bands["Start"] < end)].sort_values("Start")
    if b.empty:
        return ""
    first, last = b["Band"].iloc[0], b["Band"].iloc[-1]
    return f"{chrom}{first}" if first == last else f"{chrom}{first}-{last}"


def overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def reportable(segs, min_probes):
    return segs[segs["Call"].isin(ALTERED) & (segs["Num_Probes"] >= min_probes)].copy()


# -----------------------------------------------------------------------------
# Arm-level events
# -----------------------------------------------------------------------------
def arm_calls(segs, arms, min_probes, min_fraction, sex):
    rows = []
    for _, arm in arms.iterrows():
        chrom = arm["Chromosome"]
        if chrom == "Y" and sex == "female":
            continue
        s = segs[(segs["Chromosome"] == chrom) & (segs["End"] > arm["Start"]) & (segs["Start"] < arm["End"])]
        if s.empty:
            rows.append({"Chromosome": chrom, "Arm": arm["Arm"], "Call": "Not covered"})
            continue
        ov = s.apply(lambda r: overlap(r["Start"], r["End"], arm["Start"], arm["End"]), axis=1)
        total = ov.sum()
        # Bins inside the arm (approximated by the probes of overlapping segments, weighted by overlap)
        probes = float((s["Num_Probes"] * ov / (s["End"] - s["Start"]).clip(lower=1)).sum())
        gain = ov[(s["Call"] == "Gain") & (s["Num_Probes"] >= min_probes)].sum() / total if total else 0
        loss = ov[(s["Call"] == "Loss") & (s["Num_Probes"] >= min_probes)].sum() / total if total else 0
        if probes < min_probes:
            call = "Not assessable"
        elif gain >= min_fraction and gain >= loss:
            call = "Gain"
        elif loss >= min_fraction:
            call = "Loss"
        else:
            call = "Neutral"
        rows.append({
            "Chromosome": chrom, "Arm": arm["Arm"], "Call": call,
            "Gain_fraction": round(gain, 3), "Loss_fraction": round(loss, 3),
            "Mean_log2": round(float(np.average(s["Log2"], weights=ov)) if total else np.nan, 3),
            "Bins": int(round(probes)),
        })
    return pd.DataFrame(rows)


def consensus_arms(arm_df):
    """Arm calls kept only when both callers agree; arms assessable by a single
    caller (e.g. chrX/chrY, excluded from GATK) are left out of the consensus."""
    rows = []
    for _, r in arm_df.iterrows():
        a, b = r.get("CNVkit_Call"), r.get("GATK_Call")
        if a not in ("Gain", "Loss", "Neutral") or b not in ("Gain", "Loss", "Neutral"):
            call = "Not assessable"
        else:
            call = a if a == b else "Neutral"
        rows.append({"Chromosome": r["Chromosome"], "Arm": r["Arm"], "Call": call})
    return pd.DataFrame(rows)


def cytogenetic_summary(arm_df):
    """'+8, -7, 5q-, 17p-' style summary from arm calls."""
    if arm_df is None or arm_df.empty:
        return ""
    events = []
    for chrom, g in arm_df.groupby("Chromosome", sort=False):
        assessable = g[g["Call"].isin(("Gain", "Loss", "Neutral"))]
        if assessable.empty:
            continue
        calls = set(assessable["Call"])
        # Whole chromosome: every assessable arm altered in the same direction
        if len(calls) == 1 and calls <= set(ALTERED) and (len(assessable) == len(g[g["Call"] != "Not covered"])):
            events.append(("+" if "Gain" in calls else "-") + chrom)
            continue
        for _, r in assessable.iterrows():
            if r["Call"] in ALTERED:
                events.append(f"{chrom}{r['Arm']}{'+' if r['Call'] == 'Gain' else '-'}")
    return ", ".join(events) if events else "No arm-level events"


# -----------------------------------------------------------------------------
# Segment concordance and genes
# -----------------------------------------------------------------------------
def supported(seg, other):
    """>= 50% reciprocal overlap with a segment of the other caller in the same direction."""
    if other is None:
        return "NA"
    o = other[(other["Chromosome"] == seg["Chromosome"]) & (other["Call"] == seg["Call"])]
    for _, r in o.iterrows():
        ov = overlap(seg["Start"], seg["End"], r["Start"], r["End"])
        if ov >= 0.5 * (seg["End"] - seg["Start"]) and ov >= 0.5 * (r["End"] - r["Start"]):
            return "Yes"
    return "No"


def segment_table(segs, other, caller, report_genes, min_probes):
    alt = reportable(segs, min_probes)
    if alt.empty:
        return pd.DataFrame()
    gene_set = set(report_genes["Gene"])
    alt["Caller"] = caller
    alt["Size_Mb"] = ((alt["End"] - alt["Start"]) / 1e6).round(3)
    alt["Report_genes"] = alt["Gene"].map(
        lambda g: ",".join(x for x in str(g).split(",") if x in gene_set))
    alt["N_genes"] = alt["Gene"].map(lambda g: len([x for x in str(g).split(",") if x]))
    alt["Supported_by_other_caller"] = alt.apply(lambda r: supported(r, other), axis=1)
    cols = ["Caller", "Chromosome", "Start", "End", "Size_Mb", "Num_Probes", "Log2", "CN", "Call",
            "Cytoband", "Arm", "Report_genes", "N_genes", "Supported_by_other_caller"]
    return alt[[c for c in cols if c in alt.columns]]


def gene_table(report_genes, segs_by_caller):
    rows = []
    for _, g in report_genes.iterrows():
        row = {"Gene": g["Gene"], "Chromosome": g["Chromosome"], "Start": g["Start"], "End": g["End"]}
        for caller, segs in segs_by_caller.items():
            s = segs[(segs["Chromosome"] == g["Chromosome"]) & (segs["End"] > g["Start"]) & (segs["Start"] < g["End"])]
            if s.empty:
                row[f"{caller}_call"] = "Not covered"
                continue
            best = s.loc[s.apply(lambda r: overlap(r["Start"], r["End"], g["Start"], g["End"]), axis=1).idxmax()]
            row[f"{caller}_call"] = best["Call"]
            row[f"{caller}_log2"] = round(float(best["Log2"]), 3)
            if "CN" in best.index and pd.notna(best["CN"]):
                row[f"{caller}_CN"] = best["CN"]
            row[f"{caller}_segment_bins"] = int(best["Num_Probes"])
            row[f"{caller}_segment"] = f"{best['Chromosome']}:{int(best['Start'])}-{int(best['End'])}"
        rows.append(row)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    call_cols = [c for c in df.columns if c.endswith("_call")]
    df["Altered"] = df[call_cols].isin(ALTERED).any(axis=1)
    if len(call_cols) == 2:
        df["Concordance"] = df.apply(lambda r: concordance(r[call_cols[0]], r[call_cols[1]]), axis=1)
    df["_key"] = df["Chromosome"].map(chrom_key)
    df = df.sort_values(["Altered", "_key", "Start"], ascending=[False, True, True]).drop(columns=["_key"])
    return df


def concordance(a, b):
    a_alt, b_alt = a in ALTERED, b in ALTERED
    if a_alt and b_alt:
        return "Both" if a == b else "Discordant"
    if a_alt:
        return "CNVkit only"
    if b_alt:
        return "GATK only"
    return "-"


# -----------------------------------------------------------------------------
# Allelic imbalance / LOH
# -----------------------------------------------------------------------------
def allelic_imbalance(gatk_model, cnvkit, bands, report_genes, max_maf, min_hets, neutral, min_probes):
    rows = []
    if gatk_model:
        m = pd.read_csv(gatk_model, sep="\t", comment="@")
        m["CONTIG"] = strip_chr(m["CONTIG"])
        m = m[(m["NUM_POINTS_ALLELE_FRACTION"] >= min_hets) & (m["MINOR_ALLELE_FRACTION_POSTERIOR_50"] <= max_maf)]
        for _, r in m.iterrows():
            ratio = 2 ** r["LOG2_COPY_RATIO_POSTERIOR_50"]
            state = ("Copy-neutral LOH" if neutral[0] <= ratio <= neutral[1]
                     else "LOH with loss" if ratio < neutral[0] else "Imbalanced gain")
            rows.append({
                "Caller": "GATK", "Chromosome": r["CONTIG"], "Start": int(r["START"]), "End": int(r["END"]),
                "Size_Mb": round((r["END"] - r["START"]) / 1e6, 3), "Het_SNPs": int(r["NUM_POINTS_ALLELE_FRACTION"]),
                "Minor_allele_fraction": round(float(r["MINOR_ALLELE_FRACTION_POSTERIOR_50"]), 3),
                "Log2": round(float(r["LOG2_COPY_RATIO_POSTERIOR_50"]), 3), "State": state,
            })
    if cnvkit is not None and {"cn1", "cn2", "CN"} <= set(cnvkit.columns):
        c = cnvkit[(cnvkit["cn2"] == 0) & (cnvkit["CN"] >= 1) & (cnvkit["Num_Probes"] >= min_probes)
                   & ~cnvkit["Chromosome"].isin(["X", "Y"])]
        for _, r in c.iterrows():
            state = "Copy-neutral LOH" if r["CN"] == 2 else "LOH with loss" if r["CN"] < 2 else "Imbalanced gain"
            rows.append({
                "Caller": "CNVkit", "Chromosome": r["Chromosome"], "Start": int(r["Start"]), "End": int(r["End"]),
                "Size_Mb": round((r["End"] - r["Start"]) / 1e6, 3), "Het_SNPs": np.nan,
                "Minor_allele_fraction": round(float(r["baf"]), 3) if "baf" in r and pd.notna(r["baf"]) else np.nan,
                "Log2": round(float(r["Log2"]), 3), "State": state,
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["Cytoband"] = df.apply(lambda r: band_range(bands, r["Chromosome"], r["Start"], r["End"]), axis=1)
    df["Report_genes"] = df.apply(lambda r: ",".join(report_genes.loc[
        (report_genes["Chromosome"] == r["Chromosome"]) & (report_genes["End"] > r["Start"])
        & (report_genes["Start"] < r["End"]), "Gene"]), axis=1)
    df["_key"] = df["Chromosome"].map(chrom_key)
    return df.sort_values(["_key", "Start", "Caller"]).drop(columns="_key")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main(sample, output, cyto, genes, cnvkit=None, gatk=None, gatk_model=None, sex="", purity="",
         min_probes=10, arm_min_fraction=0.5, loh_max_maf=0.35, loh_min_hets=10,
         neutral_bounds=(0.9, 1.1), pon_mode="", pon_dir="", version="", log=sys.stderr):
    sex = (sex or "").lower()
    arms, bands = load_arms(cyto)
    report_genes = pd.read_csv(genes, sep="\t", dtype={"Chromosome": str})
    segs = {name: df for name, df in (("CNVkit", load_segments(cnvkit)), ("GATK", load_segments(gatk)))
            if df is not None}

    # Arm level
    arm_tables, summary_rows = {}, []
    for caller, s in segs.items():
        a = arm_calls(s, arms, min_probes, arm_min_fraction, sex)
        arm_tables[caller] = a
        summary_rows.append({"Item": f"{caller} arm-level events", "Value": cytogenetic_summary(a)})
    arm_df = None
    for caller, a in arm_tables.items():
        a = a.rename(columns={c: f"{caller}_{c}" for c in a.columns if c not in ("Chromosome", "Arm")})
        arm_df = a if arm_df is None else arm_df.merge(a, on=["Chromosome", "Arm"], how="outer")
    if arm_df is not None and len(arm_tables) == 2:
        arm_df["Concordance"] = arm_df.apply(
            lambda r: concordance(r.get("CNVkit_Call"), r.get("GATK_Call")), axis=1)
        summary_rows.append({"Item": "Arm-level events found by both callers",
                             "Value": cytogenetic_summary(consensus_arms(arm_df))})
    if arm_df is not None:
        arm_df["_key"] = arm_df["Chromosome"].map(chrom_key)
        arm_df = arm_df.sort_values(["_key", "Arm"]).drop(columns="_key")

    # Segments
    seg_tables = []
    for caller, s in segs.items():
        other = [d for c, d in segs.items() if c != caller]
        seg_tables.append(segment_table(s, other[0] if other else None, caller, report_genes, min_probes))
    seg_df = pd.concat([t for t in seg_tables if not t.empty], ignore_index=True) if any(
        not t.empty for t in seg_tables) else pd.DataFrame(columns=["Caller", "Chromosome", "Start", "End", "Call"])
    for caller in segs:
        n = int((seg_df["Caller"] == caller).sum()) if not seg_df.empty else 0
        summary_rows.append({"Item": f"{caller} altered segments (>= {min_probes} bins)", "Value": n})

    # Genes and allelic imbalance
    gene_df = gene_table(report_genes, segs)
    if not gene_df.empty:
        altered = gene_df.loc[gene_df["Altered"], "Gene"].tolist()
        summary_rows.append({"Item": "Report genes in altered segments",
                             "Value": ", ".join(altered) if altered else "None"})
    ai_df = allelic_imbalance(gatk_model, segs.get("CNVkit"), bands, report_genes,
                              loh_max_maf, loh_min_hets, neutral_bounds, min_probes)
    cnloh = ai_df[ai_df["State"] == "Copy-neutral LOH"] if not ai_df.empty else ai_df
    summary_rows.append({"Item": "Copy-neutral LOH segments",
                         "Value": ", ".join(f"{r['Caller']}: {r['Cytoband'] or r['Chromosome']}"
                                            for _, r in cnloh.iterrows()) if not cnloh.empty else "None"})

    summary = pd.DataFrame(
        [{"Item": "Sample", "Value": sample},
         {"Item": "Sex (sample sheet)", "Value": sex or "unknown"},
         {"Item": "Purity (sample sheet)", "Value": purity or "not provided (no purity correction)"},
         {"Item": "Callers", "Value": ", ".join(segs)}] + summary_rows)

    params = pd.DataFrame([
        ("min_probes", min_probes), ("arm_min_fraction", arm_min_fraction),
        ("loh_max_maf", loh_max_maf), ("loh_min_hets", loh_min_hets),
        ("neutral_bounds (GATK)", f"{neutral_bounds[0]}-{neutral_bounds[1]}"),
        ("pon_mode", pon_mode), ("pon_dir", pon_dir), ("pipeline_version", version),
        ("note", "GATK excludes chrX/chrY by default; sex chromosomes are reported by CNVkit only."),
        ("note", "Arm-level calls are based on the segmented part of each arm covered by the panel."),
    ], columns=["Parameter", "Value"])

    with pd.ExcelWriter(output, engine="openpyxl") as xw:
        summary.to_excel(xw, sheet_name="Summary", index=False)
        (arm_df if arm_df is not None else pd.DataFrame()).to_excel(xw, sheet_name="Arm_level", index=False)
        seg_df.to_excel(xw, sheet_name="Segments", index=False)
        gene_df.to_excel(xw, sheet_name="Genes", index=False)
        ai_df.to_excel(xw, sheet_name="Allelic_imbalance", index=False)
        params.to_excel(xw, sheet_name="Parameters", index=False)
        for ws in xw.book.worksheets:
            for col in ws.columns:
                width = max(len(str(c.value)) if c.value is not None else 0 for c in col)
                ws.column_dimensions[col[0].column_letter].width = min(max(10, width + 2), 60)

    log.write(f"CNV report written to {output}\n")
    for _, r in summary.iterrows():
        log.write(f"  {r['Item']}: {r['Value']}\n")


if __name__ == "__main__":
    if "snakemake" in locals() or "snakemake" in globals():
        with open(snakemake.log[0], "w") as log_handle:
            p = snakemake.params
            main(
                sample=p.sample, output=snakemake.output[0],
                cyto=snakemake.input.cyto, genes=snakemake.input.genes,
                cnvkit=snakemake.input.get("cnvkit"), gatk=snakemake.input.get("gatk"),
                gatk_model=snakemake.input.get("gatk_model"),
                sex=p.sex, purity=p.purity, min_probes=int(p.min_probes),
                arm_min_fraction=float(p.arm_min_fraction), loh_max_maf=float(p.loh_max_maf),
                loh_min_hets=int(p.loh_min_hets), neutral_bounds=[float(x) for x in p.neutral_bounds],
                pon_mode=p.pon_mode, pon_dir=p.pon_dir, version=p.version, log=log_handle,
            )
    else:
        parser = argparse.ArgumentParser(description="Per-sample CNV report (CNVkit + GATK)")
        parser.add_argument("--sample", required=True)
        parser.add_argument("--cyto", required=True)
        parser.add_argument("--genes", required=True, help="Report genes TSV (05_CNV_Report_Genes.py)")
        parser.add_argument("--cnvkit", help="CNVkit annotated segments")
        parser.add_argument("--gatk", help="GATK annotated segments")
        parser.add_argument("--gatk_model", help="GATK modelFinal.seg")
        parser.add_argument("--sex", default="")
        parser.add_argument("--purity", default="")
        parser.add_argument("--min_probes", type=int, default=10)
        parser.add_argument("--arm_min_fraction", type=float, default=0.5)
        parser.add_argument("--loh_max_maf", type=float, default=0.35)
        parser.add_argument("--loh_min_hets", type=int, default=10)
        parser.add_argument("--output", required=True)
        a = parser.parse_args()
        main(a.sample, a.output, a.cyto, a.genes, a.cnvkit, a.gatk, a.gatk_model, a.sex, a.purity,
             a.min_probes, a.arm_min_fraction, a.loh_max_maf, a.loh_min_hets)
