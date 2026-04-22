"""
Repository usage summary
Description: Generate additional manuscript-oriented figures and tables from processed result tables.
Usage: Run after the main processed text tables have been generated; pass --data-dir and optionally --outdir.
Outputs: Additional Excel tables and PNG summary figures for manuscript support files.
Note: Update local input/output paths or pass explicit command-line arguments before running.

Public repository version: file paths are intentionally left blank or configurable.
This script was lightly sanitized for sharing and may require project-specific inputs.
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Additional figures and tables for EOAD vs LOAD vs HC manuscript.
This script is designed to be simple, editable, and not path-hardcoded.
It does NOT run automatically here; the user can download and run locally.

Expected files in data_dir:
- roi_level_coupling_stats_age_sex_fdr.txt
- coupling_ancova_table_age_sex.txt
- coupling_ancova_coefficients_age_sex.txt
- gretna_bonferroni_master_summary_autoTF.txt
- gretna_bonferroni_detail_summary_autoTF.txt
- roi_level_sdi_stats_age_sex.txt            (optional)
- roi_level_fa_md_stats_age_sex.txt          (optional)
- run_coupling_sensitivity_and_clinical_minimal_fix2.txt (optional; only for parsing hints)

Outputs in outdir:
- Table_S1_full_roi_coupling_results.xlsx
- Table_S2_mean_coupling_model.xlsx
- Table_S3_full_graph_bonferroni_results.xlsx
- Figure_S1_frontal_roi_effect_sizes.png
- Figure_S2_graph_metric_counts.png
- Figure_S3_negative_results_summary.png   (if optional files exist)
- Figure_S4_top_graph_regions.png
"""

from __future__ import annotations
import argparse
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_text(path: Path) -> str:
    return path.read_text(encoding='utf-8', errors='ignore')


def smart_read_table(path: Path) -> pd.DataFrame:
    for sep in ["\t", ",", None]:
        try:
            if sep is None:
                return pd.read_csv(path, sep=None, engine="python")
            return pd.read_csv(path, sep=sep)
        except Exception:
            continue
    raise ValueError(f"Could not read table: {path}")


ROI_PATTERNS = [
    r'Frontal[_\s-].*', r'.*Frontal.*', r'.*Superior_Frontal.*', r'.*Middle_Frontal.*',
    r'.*Inferior_Frontal.*', r'.*Orbital.*', r'.*Precentral.*', r'.*Supplementary_Motor.*'
]


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    return df


def find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    for c in df.columns:
        cl = c.lower()
        if any(cand.lower() in cl for cand in candidates):
            return c
    return None


def load_roi_coupling(path: Path) -> pd.DataFrame:
    df = normalize_columns(smart_read_table(path))
    roi_col = find_col(df, ["ROI", "Region", "roi_name", "AAL", "Label"])
    q_col = find_col(df, ["Q", "q", "FDR_q", "q_value", "Q_raw", "Q_group_2", "Q_group_3"])
    p_col = find_col(df, ["P", "p", "p_value", "P_group_2", "P_group_3"])
    beta_col = find_col(df, ["Beta", "beta", "Estimate", "Coef", "Coefficient", "group_2_beta", "group_3_beta"])

    # Keep original full table for export.
    out = df.copy()

    # Try to identify EOAD vs LOAD comparison.
    contrast_col = find_col(df, ["contrast", "comparison", "term"])
    if contrast_col is not None:
        mask = df[contrast_col].astype(str).str.contains("LOAD|EOAD|group_2|group2", case=False, na=False)
        if mask.any():
            df = df.loc[mask].copy()

    if roi_col is None:
        # Fallback: assume first column is ROI-like.
        roi_col = df.columns[0]
    df["ROI_plot"] = df[roi_col].astype(str)

    if q_col is not None:
        q = pd.to_numeric(df[q_col], errors="coerce")
    else:
        q = pd.Series(np.nan, index=df.index)
    if p_col is not None:
        p = pd.to_numeric(df[p_col], errors="coerce")
    else:
        p = pd.Series(np.nan, index=df.index)
    if beta_col is not None:
        beta = pd.to_numeric(df[beta_col], errors="coerce")
    else:
        beta = pd.Series(np.nan, index=df.index)

    df["Q_plot"] = q
    df["P_plot"] = p
    df["Beta_plot"] = beta
    return out, df[["ROI_plot", "Q_plot", "P_plot", "Beta_plot"]].copy()


def select_frontal_hits(df: pd.DataFrame, q_threshold: float = 0.05) -> pd.DataFrame:
    mask_q = df["Q_plot"].notna() & (df["Q_plot"] < q_threshold)
    mask_frontal = pd.Series(False, index=df.index)
    for pat in ROI_PATTERNS:
        mask_frontal = mask_frontal | df["ROI_plot"].str.contains(pat, case=False, regex=True, na=False)
    hits = df.loc[mask_q & mask_frontal].copy()
    if hits.empty:
        hits = df.loc[mask_frontal].copy()
    hits = hits.sort_values(["Q_plot", "P_plot"], na_position="last").head(8)
    return hits


def plot_frontal_effect_sizes(df: pd.DataFrame, outpath: Path) -> None:
    if df.empty:
        return
    plot_df = df.copy().sort_values("Beta_plot")
    fig, ax = plt.subplots(figsize=(8, max(4, 0.5 * len(plot_df))))
    ax.barh(plot_df["ROI_plot"], plot_df["Beta_plot"])
    ax.set_xlabel("Effect size (LOAD vs EOAD)")
    ax.set_ylabel("Region")
    ax.set_title("Frontal regional structural-functional coupling")
    for i, (_, row) in enumerate(plot_df.iterrows()):
        qtxt = "NA" if pd.isna(row["Q_plot"]) else f"q={row['Q_plot']:.3f}"
        ax.text(row["Beta_plot"], i, f"  {qtxt}", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_ancova_table(path: Path) -> pd.DataFrame:
    txt = read_text(path)
    lines = [ln.strip() for ln in txt.splitlines() if ln.strip()]
    rows = []
    for ln in lines:
        if re.search(r'group|age|sex|model|error|residual', ln, re.I):
            parts = re.split(r'\s{2,}|\t', ln)
            if len(parts) >= 2:
                rows.append(parts)
    if not rows:
        return pd.DataFrame({"raw_line": lines})
    maxlen = max(len(r) for r in rows)
    rows = [r + [""] * (maxlen - len(r)) for r in rows]
    cols = [f"col{i+1}" for i in range(maxlen)]
    return pd.DataFrame(rows, columns=cols)


def parse_coeff_table(path: Path) -> pd.DataFrame:
    txt = read_text(path)
    lines = [ln.strip() for ln in txt.splitlines() if ln.strip()]
    rows = []
    for ln in lines:
        if re.search(r'intercept|group|age|sex', ln, re.I):
            parts = re.split(r'\s{2,}|\t', ln)
            if len(parts) >= 2:
                rows.append(parts)
    if not rows:
        return pd.DataFrame({"raw_line": lines})
    maxlen = max(len(r) for r in rows)
    rows = [r + [""] * (maxlen - len(r)) for r in rows]
    cols = [f"col{i+1}" for i in range(maxlen)]
    return pd.DataFrame(rows, columns=cols)


def load_graph_detail(path: Path) -> pd.DataFrame:
    txt = read_text(path)
    # Flexible parser: metric header followed by region lines.
    blocks = []
    current_metric = None
    for ln in txt.splitlines():
        s = ln.strip()
        if not s:
            continue
        if re.match(r'^(metric|measure)\s*[:=]', s, re.I):
            current_metric = re.split(r'[:=]', s, maxsplit=1)[1].strip()
            continue
        if re.search(r'betweenness|degree|efficiency|clustering', s, re.I) and len(s.split()) < 8:
            current_metric = s
            continue
        if current_metric is None:
            continue
        parts = re.split(r'\s{2,}|\t', s)
        if len(parts) >= 1:
            blocks.append([current_metric] + parts)
    if not blocks:
        # Fallback to generic table read.
        df = normalize_columns(smart_read_table(path))
        metric_col = find_col(df, ["metric", "measure"])
        if metric_col is None:
            df["metric"] = "graph_metric"
        else:
            df.rename(columns={metric_col: "metric"}, inplace=True)
        return df
    maxlen = max(len(r) for r in blocks)
    blocks = [r + [""] * (maxlen - len(r)) for r in blocks]
    cols = ["metric"] + [f"col{i}" for i in range(1, maxlen)]
    df = pd.DataFrame(blocks, columns=cols)
    return df


def summarize_graph_counts(df: pd.DataFrame) -> pd.DataFrame:
    if "metric" not in df.columns:
        df["metric"] = "graph_metric"
    counts = df.groupby("metric").size().reset_index(name="n_regions")
    return counts.sort_values("n_regions", ascending=False)


def plot_graph_counts(df: pd.DataFrame, outpath: Path) -> None:
    if df.empty:
        return
    plot_df = df.sort_values("n_regions")
    fig, ax = plt.subplots(figsize=(8, max(4, 0.6 * len(plot_df))))
    ax.barh(plot_df["metric"], plot_df["n_regions"])
    ax.set_xlabel("Number of significant regions")
    ax.set_ylabel("Graph metric")
    ax.set_title("Bonferroni-significant nodal abnormalities")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_top_graph_regions(df: pd.DataFrame, outpath: Path) -> None:
    if df.empty:
        return
    data = df.copy()
    # Use first non-metric column as region label.
    region_col = next((c for c in data.columns if c != "metric"), None)
    if region_col is None:
        return
    data["region"] = data[region_col].astype(str)
    top = data.groupby("region").size().reset_index(name="count").sort_values("count", ascending=False).head(15)
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(top["region"].iloc[::-1], top["count"].iloc[::-1])
    ax.set_xlabel("Number of significant graph metrics")
    ax.set_ylabel("Region")
    ax.set_title("Top recurrent graph-theoretical regions")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def summarize_negative_results(sdi_path: Path | None, dti_path: Path | None) -> pd.DataFrame:
    rows = []
    for name, path in [("SDI", sdi_path), ("DTI", dti_path)]:
        if path is None or not path.exists():
            continue
        txt = read_text(path)
        sig = len(re.findall(r'\b(q|Q)\s*[<=>]\s*0?\.0?5', txt))
        rows.append({"domain": name, "reported_corrected_significant_hits": sig})
    return pd.DataFrame(rows)


def plot_negative_summary(df: pd.DataFrame, outpath: Path) -> None:
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(df["domain"], df["reported_corrected_significant_hits"])
    ax.set_ylabel("Reported corrected significant hits")
    ax.set_title("Summary of secondary negative-result domains")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default=".", help="Directory containing input files")
    parser.add_argument("--outdir", default="paper_outputs_additional", help="Output directory")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    roi_file = data_dir / "roi_level_coupling_stats_age_sex_fdr.txt"
    ancova_file = data_dir / "coupling_ancova_table_age_sex.txt"
    coef_file = data_dir / "coupling_ancova_coefficients_age_sex.txt"
    graph_file = data_dir / "gretna_bonferroni_detail_summary_autoTF.txt"
    sdi_file = data_dir / "roi_level_sdi_stats_age_sex.txt"
    dti_file = data_dir / "run_roi_level_fa_md_stats_age_sex.txt"
    if not dti_file.exists():
        dti_file = data_dir / "run_roi_level_dti_stats_minimal.txt"

    # Table S1 + Figure S1
    if roi_file.exists():
        roi_full, roi_plot = load_roi_coupling(roi_file)
        roi_full.to_excel(outdir / "Table_S1_full_roi_coupling_results.xlsx", index=False)
        frontal = select_frontal_hits(roi_plot)
        if not frontal.empty:
            frontal.to_excel(outdir / "Table_S1b_frontal_roi_hits.xlsx", index=False)
            plot_frontal_effect_sizes(frontal, outdir / "Figure_S1_frontal_roi_effect_sizes.png")

    # Table S2 mean coupling model
    ancova_df = parse_ancova_table(ancova_file) if ancova_file.exists() else pd.DataFrame()
    coef_df = parse_coeff_table(coef_file) if coef_file.exists() else pd.DataFrame()
    with pd.ExcelWriter(outdir / "Table_S2_mean_coupling_model.xlsx") as writer:
        ancova_df.to_excel(writer, index=False, sheet_name="ancova")
        coef_df.to_excel(writer, index=False, sheet_name="coefficients")

    # Table S3 + Figure S2/S4
    if graph_file.exists():
        graph_df = load_graph_detail(graph_file)
        graph_df.to_excel(outdir / "Table_S3_full_graph_bonferroni_results.xlsx", index=False)
        counts = summarize_graph_counts(graph_df)
        counts.to_excel(outdir / "Table_S3b_graph_metric_counts.xlsx", index=False)
        plot_graph_counts(counts, outdir / "Figure_S2_graph_metric_counts.png")
        plot_top_graph_regions(graph_df, outdir / "Figure_S4_top_graph_regions.png")

    # Negative summary if possible
    neg_df = summarize_negative_results(sdi_file if sdi_file.exists() else None,
                                       dti_file if dti_file.exists() else None)
    if not neg_df.empty:
        neg_df.to_excel(outdir / "Table_S4_negative_results_summary.xlsx", index=False)
        plot_negative_summary(neg_df, outdir / "Figure_S3_negative_results_summary.png")

    print(f"Done. Outputs saved to: {outdir.resolve()}")


if __name__ == "__main__":
    main()
