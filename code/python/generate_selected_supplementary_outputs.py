"""
Repository usage summary
Description: Generate selected supplementary tables and figures from processed summary files.
Usage: Run after the relevant processed summary files are available; pass --data-dir.
Outputs: Supplementary Excel tables and PNG figure files.
Note: Update local input/output paths or pass explicit command-line arguments before running.

Public repository version: file paths are intentionally left blank or configurable.
This script was lightly sanitized for sharing and may require project-specific inputs.
"""

# -*- coding: utf-8 -*-
"""
Generate selected supplementary figures/tables for the EOAD/LOAD/HC paper.

Requested outputs:
- Figure S2
- Figure S4
- Table S1
- Table S2
- Table S3
- Table S3b

Default data directory:
    

Required files in :
    roi_level_coupling_stats_age_sex_fdr.txt
    coupling_ancova_table_age_sex.txt
    coupling_ancova_coefficients_age_sex.txt
    gretna_bonferroni_detail_summary_autoTF.txt
"""

from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

DEFAULT_DATA_DIR = ""
OUT_DIR_NAME = "paper_outputs_sup"


def ensure_outdir(base: Path) -> Path:
    outdir = base / OUT_DIR_NAME
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


def read_txt_table(path: Path) -> pd.DataFrame:
    for sep in ["\t", ",", r"\s{2,}", r"\s+"]:
        try:
            df = pd.read_csv(path, sep=sep, engine="python")
            if df.shape[1] >= 2:
                return df
        except Exception:
            continue
    raise ValueError(f"Cannot parse table: {path}")


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    return df


def find_col(df: pd.DataFrame, candidates):
    for cand in candidates:
        for c in df.columns:
            if c.lower() == cand.lower():
                return c
    for cand in candidates:
        for c in df.columns:
            if cand.lower() in c.lower():
                return c
    return None


def load_roi_coupling_table(data_dir: Path) -> pd.DataFrame:
    path = data_dir / "roi_level_coupling_stats_age_sex_fdr.txt"
    df = normalize_columns(read_txt_table(path))

    roi_col = find_col(df, ["ROI", "Region", "roi", "region"])
    hemi_col = find_col(df, ["Hemisphere", "hemi", "Side"])
    contrast_col = find_col(df, ["Contrast", "comparison"])
    beta_col = find_col(df, ["Beta", "Estimate", "Coef", "Coefficient"])
    p_col = find_col(df, ["P", "p", "PValue", "p_value", "raw_p"])
    q_col = find_col(df, ["Q", "q", "FDR", "q_value", "adj_p"])

    if roi_col is None:
        roi_col = df.columns[0]

    out = pd.DataFrame()
    out["ROI"] = df[roi_col]
    out["Hemisphere"] = df[hemi_col] if hemi_col else ""
    out["Contrast"] = df[contrast_col] if contrast_col else "LOAD vs EOAD"
    out["Beta"] = pd.to_numeric(df[beta_col], errors="coerce") if beta_col else np.nan
    out["P"] = pd.to_numeric(df[p_col], errors="coerce") if p_col else np.nan
    out["Q"] = pd.to_numeric(df[q_col], errors="coerce") if q_col else np.nan
    out["Direction"] = np.where(out["Beta"] > 0, "Higher in first-named group", "Higher in second-named group")
    return out


def save_table_s1(df: pd.DataFrame, outdir: Path):
    df.to_excel(outdir / "Table_S1_full_roi_coupling_results.xlsx", index=False)


def load_mean_coupling_tables(data_dir: Path):
    table_path = data_dir / "coupling_ancova_table_age_sex.txt"
    coef_path = data_dir / "coupling_ancova_coefficients_age_sex.txt"
    table_df = normalize_columns(read_txt_table(table_path))
    coef_df = normalize_columns(read_txt_table(coef_path))
    return table_df, coef_df


def save_table_s2(table_df: pd.DataFrame, coef_df: pd.DataFrame, outdir: Path):
    with pd.ExcelWriter(outdir / "Table_S2_mean_coupling_model.xlsx") as writer:
        table_df.to_excel(writer, sheet_name="ANCOVA_table", index=False)
        coef_df.to_excel(writer, sheet_name="Coefficients", index=False)


def load_gretna_detail(data_dir: Path) -> pd.DataFrame:
    path = data_dir / "gretna_bonferroni_detail_summary_autoTF.txt"
    df = normalize_columns(read_txt_table(path))

    metric_col = find_col(df, ["Metric", "measure", "index"])
    roi_col = find_col(df, ["ROI", "Region", "Node"])
    contrast_col = find_col(df, ["Contrast", "comparison"])
    p_col = find_col(df, ["P", "p", "PValue", "p_value"])
    q_col = find_col(df, ["Q", "q", "Bonferroni", "adj_p"])
    beta_col = find_col(df, ["Beta", "Estimate", "Coef", "Coefficient"])

    out = pd.DataFrame()
    out["Metric"] = df[metric_col] if metric_col else ""
    out["ROI"] = df[roi_col] if roi_col else df.iloc[:, 0]
    out["Contrast"] = df[contrast_col] if contrast_col else ""
    out["P"] = pd.to_numeric(df[p_col], errors="coerce") if p_col else np.nan
    out["AdjustedP"] = pd.to_numeric(df[q_col], errors="coerce") if q_col else np.nan
    out["Beta"] = pd.to_numeric(df[beta_col], errors="coerce") if beta_col else np.nan
    return out


def save_table_s3(df: pd.DataFrame, outdir: Path):
    df.to_excel(outdir / "Table_S3_full_graph_bonferroni_results.xlsx", index=False)


def save_table_s3b(df: pd.DataFrame, outdir: Path):
    counts = (
        df.groupby("Metric", dropna=False)
          .size()
          .reset_index(name="Significant_ROI_Count")
          .sort_values("Significant_ROI_Count", ascending=False)
    )
    counts.to_excel(outdir / "Table_S3b_graph_metric_counts.xlsx", index=False)
    return counts


def make_figure_s2(counts: pd.DataFrame, outdir: Path):
    plt.figure(figsize=(8, 5))
    plt.bar(counts["Metric"].astype(str), counts["Significant_ROI_Count"])
    plt.xticks(rotation=30, ha="right")
    plt.ylabel("Number of significant regions")
    plt.title("Figure S2. Significant graph findings by metric")
    plt.tight_layout()
    plt.savefig(outdir / "Figure_S2_graph_metric_counts.png", dpi=300)
    plt.close()


def make_figure_s4(df: pd.DataFrame, outdir: Path):
    top = (
        df.groupby("ROI", dropna=False)
          .size()
          .reset_index(name="Count")
          .sort_values("Count", ascending=False)
          .head(15)
    )
    plt.figure(figsize=(8, 6))
    plt.barh(top["ROI"].astype(str)[::-1], top["Count"][::-1])
    plt.xlabel("Number of significant appearances")
    plt.title("Figure S4. Top graph regions across significant metrics")
    plt.tight_layout()
    plt.savefig(outdir / "Figure_S4_top_graph_regions.png", dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default=".", help="Directory containing input files")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    outdir = ensure_outdir(data_dir)

    roi_df = load_roi_coupling_table(data_dir)
    save_table_s1(roi_df, outdir)

    ancova_df, coef_df = load_mean_coupling_tables(data_dir)
    save_table_s2(ancova_df, coef_df, outdir)

    gretna_df = load_gretna_detail(data_dir)
    save_table_s3(gretna_df, outdir)
    counts = save_table_s3b(gretna_df, outdir)
    make_figure_s2(counts, outdir)
    make_figure_s4(gretna_df, outdir)

    print("Done.")
    print(f"Outputs saved to: {outdir}")


if __name__ == "__main__":
    main()
