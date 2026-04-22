"""
Repository usage summary
Description: Generate clean paper-style tables and figures from selected processed cohort files.
Usage: Run after core processed tables are available; pass --data-dir and optionally --out-dir.
Outputs: Clean CSV/XLSX tables and PNG summary figures suitable for manuscript drafting.
Note: Update local input/output paths or pass explicit command-line arguments before running.

Public repository version: file paths are intentionally left blank or configurable.
This script was lightly sanitized for sharing and may require project-specific inputs.
"""

import argparse
from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROI_LABELS = {
    8: 'Frontal_Mid_R',
    10: 'Frontal_Mid_Orb_R',
    14: 'Frontal_Inf_Tri_R',
    23: 'Frontal_Sup_Medial_L',
}
GROUP_MAP = {1: 'EOAD', 2: 'LOAD', 3: 'HC', 'group1': 'EOAD', 'group2': 'LOAD', 'group3': 'HC'}
CLIN_COLS = ['TOTSCORE', 'CDRSB', 'FAQTOTAL', 'GDTOTAL', 'MMSCORE', 'MOCA']
CLIN_LABELS = {
    'TOTSCORE': 'ADAS-Cog total score',
    'CDRSB': 'CDR-SB',
    'FAQTOTAL': 'FAQ total score',
    'GDTOTAL': 'GDS total score',
    'MMSCORE': 'MMSE',
    'MOCA': 'MoCA',
}


def read_inputs(base: Path):
    files = {
        'biao': base / 'biao.xlsx',
        'coupling': base / 'sc_fc_coupling_summary.txt',
        'roi_fdr': base / 'roi_level_coupling_stats_age_sex_fdr.txt',
        'gretna': base / 'gretna_bonferroni_detail_summary_autoTF.txt',
        'matched': base / 'matched_sc_fc_subjects.txt',
    }
    missing = [str(p.name) for p in files.values() if not p.exists()]
    if missing:
        raise FileNotFoundError('Missing required files: ' + ', '.join(missing))

    biao = pd.read_excel(files['biao'])
    biao = biao.rename(columns={'组别': 'group', '新编号': 'sub', '年龄': 'age', '性别': 'sex'})
    biao['group_label'] = biao['group'].map(GROUP_MAP)
    matched = pd.read_csv(files['matched'], sep='\t')
    matched['group_num'] = matched['group'].str.extract(r'(\d+)').astype(int)
    matched['group_label'] = matched['group'].map(GROUP_MAP)
    coupling = pd.read_csv(files['coupling'], sep='\t')
    coupling['group_num'] = coupling['group'].str.extract(r'(\d+)').astype(int)
    coupling['group_label'] = coupling['group'].map(GROUP_MAP)
    roi_fdr = pd.read_csv(files['roi_fdr'], sep='\t')
    gretna = pd.read_csv(files['gretna'], sep='\t')
    return biao, matched, coupling, roi_fdr, gretna


def mean_sd(x):
    x = pd.to_numeric(x, errors='coerce').dropna()
    return f"{x.mean():.2f} ± {x.std(ddof=1):.2f}" if len(x) else ''


def build_table1_panel(df: pd.DataFrame):
    rows = []
    groups = [('EOAD', 1), ('LOAD', 2), ('HC', 3)]

    def subdf(g):
        return df[df['group'] == g]

    rows.append(['N'] + [len(subdf(g)) for _, g in groups])
    rows.append(['Age, years'] + [mean_sd(subdf(g)['age']) for _, g in groups])
    for sex_value, label in [(1, 'Male, n (%)'), (2, 'Female, n (%)')]:
        vals = []
        for _, g in groups:
            d = subdf(g)
            n = int((d['sex'] == sex_value).sum())
            vals.append(f"{n} ({n/len(d)*100:.1f})" if len(d) else '')
        rows.append([label] + vals)
    for col in CLIN_COLS:
        rows.append([CLIN_LABELS[col]] + [mean_sd(subdf(g)[col]) for _, g in groups])
    out = pd.DataFrame(rows, columns=['Variable', 'EOAD', 'LOAD', 'HC'])
    return out


def save_table1(biao, matched, outdir: Path):
    panel_a = build_table1_panel(biao.copy())
    matched_ids = set(matched['sub'])
    panel_b_df = biao[biao['sub'].isin(matched_ids)].copy()
    panel_b = build_table1_panel(panel_b_df)
    xlsx = outdir / 'Table1_demographics.xlsx'
    with pd.ExcelWriter(xlsx, engine='openpyxl') as writer:
        panel_a.to_excel(writer, index=False, sheet_name='Panel_A_functional')
        panel_b.to_excel(writer, index=False, sheet_name='Panel_B_multimodal')
    panel_a.to_csv(outdir / 'Table1_panel_A.csv', index=False, encoding='utf-8-sig')
    panel_b.to_csv(outdir / 'Table1_panel_B.csv', index=False, encoding='utf-8-sig')
    return panel_a, panel_b


def plot_mean_coupling(coupling, outdir: Path):
    data = [coupling.loc[coupling['group_label'] == g, 'mean_coupling'].dropna().values for g in ['EOAD', 'LOAD', 'HC']]
    fig, ax = plt.subplots(figsize=(6.2, 4.4), dpi=200)
    bp = ax.boxplot(data, labels=['EOAD', 'LOAD', 'HC'], patch_artist=True, widths=0.55)
    for patch in bp['boxes']:
        patch.set_facecolor('white')
    for i, vals in enumerate(data, start=1):
        x = np.random.normal(i, 0.05, size=len(vals))
        ax.scatter(x, vals, s=10, alpha=0.5)
    ax.set_ylabel('Mean structural-functional coupling')
    ax.set_xlabel('Group')
    ax.set_title('Mean structural-functional coupling by group')
    ax.axhline(0, lw=0.8)
    fig.tight_layout()
    fig.savefig(outdir / 'Figure_mean_coupling.png', bbox_inches='tight')
    plt.close(fig)


def save_coupling_table_and_plot(roi_fdr, outdir: Path):
    sig = roi_fdr[roi_fdr['fdr_pass_group2'] == 1].copy()
    sig['ROI label'] = sig['ROI'].map(ROI_LABELS).fillna(sig['ROI'].astype(str))
    sig['Contrast'] = 'LOAD vs EOAD'
    sig['Direction'] = np.where(sig['beta_group2'] > 0, 'LOAD > EOAD', 'EOAD > LOAD')
    table = sig[['ROI', 'ROI label', 'Contrast', 'beta_group2', 'p_group2', 'q_group2', 'Direction']].rename(columns={
        'beta_group2': 'Beta', 'p_group2': 'P', 'q_group2': 'Q'
    }).sort_values('Q')
    table.to_csv(outdir / 'Table2_regional_coupling_main_results.csv', index=False, encoding='utf-8-sig')
    table.to_excel(outdir / 'Table2_regional_coupling_main_results.xlsx', index=False)

    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=200)
    plot_df = table.sort_values('Beta')
    ax.barh(plot_df['ROI label'], plot_df['Beta'])
    ax.set_xlabel('Adjusted beta for LOAD vs EOAD')
    ax.set_ylabel('Region')
    ax.set_title('Significant regional structural-functional coupling differences')
    for i, (_, row) in enumerate(plot_df.iterrows()):
        ax.text(row['Beta'], i, f"  q={row['Q']:.3f}", va='center', fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / 'Figure_regional_coupling_main_results.png', bbox_inches='tight')
    plt.close(fig)
    return table


def save_graph_table_and_plot(gretna, outdir: Path):
    use = gretna[gretna['Contrast'] == 'group1_vs_group2'].copy()
    summary = use.groupby('Metric').agg(
        Significant_ROIs=('ROI', 'count'),
        Example_region=('Label', 'first'),
        Min_P=('P', 'min')
    ).reset_index().sort_values(['Significant_ROIs', 'Min_P'], ascending=[False, True])
    summary.to_csv(outdir / 'Table3_graph_nodal_summary.csv', index=False, encoding='utf-8-sig')
    summary.to_excel(outdir / 'Table3_graph_nodal_summary.xlsx', index=False)

    fig, ax = plt.subplots(figsize=(7, 4.2), dpi=200)
    plot_df = summary.sort_values('Significant_ROIs')
    ax.barh(plot_df['Metric'], plot_df['Significant_ROIs'])
    ax.set_xlabel('Number of significant ROIs')
    ax.set_ylabel('Graph metric')
    ax.set_title('Bonferroni-significant nodal graph findings\n(EOAD vs LOAD)')
    fig.tight_layout()
    fig.savefig(outdir / 'Figure_graph_nodal_summary.png', bbox_inches='tight')
    plt.close(fig)
    return summary


def write_readme(outdir: Path):
    text = (
        'Files generated by generate_clean_figures_tables.py\n\n'
        'Table1_demographics.xlsx: demographics and clinical characteristics with Panel A and Panel B.\n'
        'Figure_mean_coupling.png: boxplot of mean structural-functional coupling by group.\n'
        'Table2_regional_coupling_main_results.*: significant regional coupling results after BH-FDR.\n'
        'Figure_regional_coupling_main_results.png: bar plot of significant regional coupling regions.\n'
        'Table3_graph_nodal_summary.*: Bonferroni-significant graph metrics summary for EOAD vs LOAD.\n'
        'Figure_graph_nodal_summary.png: number of significant ROIs by graph metric.\n'
    )
    (outdir / 'README.txt').write_text(text, encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description='Generate clean paper-style tables and figures from current project files.')
    parser.add_argument('--data-dir', default='.', help='Folder containing input files. Default: current folder')
    parser.add_argument('--out-dir', default='paper_outputs_clean', help='Output folder name. Default: paper_outputs_clean')
    args = parser.parse_args()

    data_dir = Path(args.data_dir).resolve()
    outdir = Path(args.out_dir)
    if not outdir.is_absolute():
        outdir = data_dir / outdir
    outdir.mkdir(parents=True, exist_ok=True)

    biao, matched, coupling, roi_fdr, gretna = read_inputs(data_dir)
    save_table1(biao, matched, outdir)
    plot_mean_coupling(coupling, outdir)
    save_coupling_table_and_plot(roi_fdr, outdir)
    save_graph_table_and_plot(gretna, outdir)
    write_readme(outdir)
    print(f'Output folder: {outdir}')


if __name__ == '__main__':
    main()
