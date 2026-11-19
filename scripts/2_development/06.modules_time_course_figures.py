"""Figures of 2_development/03 and 04: A/B/C archetype gene modules (RNA) and chromatin region modules (ATAC),
M1-M5, over postnatal day in Yoo25 L2/3, NR vs DR, per cell group (all L2/3, A', B', C' cells).

One panel per (cell group, gene / region set, module): y = module score per sample (mean over the module's
genes / regions of log2(CPM+1) minus their P21 NR mean, all L2/3; dashed line at 0); lines through the
mean across samples per age x condition, dots = samples. Title: n = module genes / regions passing the
filter (genes: per-gene = fraction of them with per-gene DESeq2 padj < FDR, DR and slope); DR and slope =
OLS estimate (padj, BH over the 15 modules of the cell group). Marked (bracket right of P21): padj < FDR,
no effect-size cutoff. y is shared within a row. Per modality, the same panels in two layouts:
  *_by_cells    one page per cell group; rows = A/B/C sets, columns = M1-M5
  *_by_module   one page per set (A, B, C); rows = M1-M5, columns = cell groups

Reads:
  results/2_development/03.yoo25_l23_gene_modules_{samples,stats}.tsv
  results/2_development/04.yoo25_l23_atac_region_modules_{samples,stats}.tsv
Outputs:
  figures/2_development/06.yoo25_l23_gene_modules_time_course_by_{cells,module}.pdf
  figures/2_development/06.yoo25_l23_atac_region_modules_time_course_by_{cells,module}.pdf
"""

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from l23.constants import FDR
from l23.timecourse import stars, draw_nr_dr

matplotlib.rcParams['pdf.fonttype'] = 42

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES_DIR   = os.path.join(REPO_ROOT, 'results', '2_development')
FIG_DIR   = os.path.join(REPO_ROOT, 'figures', '2_development')

GROUPS     = ['all L2/3', "A' cells", "B' cells", "C' cells"]
ARCHETYPES = ['A', 'B', 'C']
MODULES    = ['M1', 'M2', 'M3', 'M4', 'M5']
PANEL_W    = 2.5   # inches
# modality -> inputs, outputs and wording
KINDS = {
    'genes': dict(
        samples='03.yoo25_l23_gene_modules_samples.tsv', stats='03.yoo25_l23_gene_modules_stats.tsv',
        out='06.yoo25_l23_gene_modules_time_course_by_{}.pdf', panel_h=3.2,
        n_title=lambda r: (f'n = {int(r["n_genes"])}; per-gene: DR {r["frac_genes_DR_14"]:.0%}, '
                           f'slope {r["frac_genes_slope_14"]:.0%}'),
        row='{} genes', page_cells='{}: A/B/C gene modules (clst 0-4 = M1-M5)',
        page_module='{} gene modules (clst 0-4 = M1-M5) by cell group',
        ylabel='mean log$_2$(CPM+1) − P21 NR mean',
        n_note=f'n = module genes passing the filter; per-gene = fraction of genes with DESeq2 padj < {FDR}'),
    'atac': dict(
        samples='04.yoo25_l23_atac_region_modules_samples.tsv', stats='04.yoo25_l23_atac_region_modules_stats.tsv',
        out='06.yoo25_l23_atac_region_modules_time_course_by_{}.pdf', panel_h=3.0,
        n_title=lambda r: f'n = {int(r["n_regions"])} regions',
        row='{} regions', page_cells='{}: A/B/C ATAC region modules (M1-M5)',
        page_module='{} ATAC region modules (M1-M5) by cell group',
        ylabel='mean log$_2$(CPM in peaks + 1) − P21 NR mean',
        n_note='n = module regions passing the filter (>= 10 counts in >= 2 samples)'),
}

os.makedirs(FIG_DIR, exist_ok=True)


def draw_panel(ax, kind, samp, st, days, grp, k, m):
    r = st.loc[(grp, k, m)]
    s = samp[(samp['cells'] == grp) & (samp['archetype'] == k) & (samp['module'] == m)]
    draw_nr_dr(ax, s, 'score', days, [f'{name} {stars(r[f"padj_{key}"])}'
                                      for name, key in [('DR', 'DR_d16'), ('slope', 'slope')] if stars(r[f'padj_{key}'])])
    ax.set_title(f'{kind["n_title"](r)}\n'
                 f'DR {r["DR_d16"]:.2f} ({r["padj_DR_d16"]:.2g})\n'
                 f'slope {r["slope_per_day"]:+.3f} ({r["padj_slope"]:.2g})', fontsize=6.5)


def draw_page(pdf, kind, samp, st, days, title, rows, cols, key):
    """One page: rows x cols panels; key(row, col) -> (cell group, set, module)."""
    ph = kind['panel_h']
    fig, axes = plt.subplots(len(rows), len(cols), figsize=(PANEL_W * len(cols), ph * len(rows)),
                             sharey='row', squeeze=False)
    for i, rl in enumerate(rows):
        for j, cl in enumerate(cols):
            draw_panel(axes[i, j], kind, samp, st, days, *key(rl, cl))
            if i == 0:
                axes[i, j].text(0.5, 1.3, cl, transform=axes[i, j].transAxes, ha='center', va='bottom', fontsize=10)
        axes[i, 0].set_ylabel(rl, fontsize=11)
    for ax in axes[-1]:
        ax.set_xlabel('postnatal day')
    axes[0, 0].legend(fontsize=8, frameon=False, loc='lower left')
    fig.supylabel(kind['ylabel'])
    fig.suptitle(f'{title}, NR vs DR; reference = P21 NR mean, all L2/3 (dashed line)\n'
                 'OLS of the module score over P12-P21, ~ age + DR + DR·(day - 16)\n'
                 'DR = average DR effect (at day 16); slope = change in DR effect per day. '
                 'Panel: estimate (padj, BH over 15 modules)\n'
                 f'marked (bracket): padj < {FDR}; * < 0.05, ** < 0.01, *** < 0.001\n'
                 f'{kind["n_note"]}', fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.9 / (ph * len(rows))))
    pdf.savefig(fig)
    plt.close(fig)


for kind in KINDS.values():
    samp = pd.read_csv(os.path.join(RES_DIR, kind['samples']), sep='\t')
    st = pd.read_csv(os.path.join(RES_DIR, kind['stats']), sep='\t').set_index(['cells', 'archetype', 'module'])
    absent = [x for x in pd.MultiIndex.from_product([GROUPS, ARCHETYPES, MODULES]) if x not in st.index]
    if absent:
        raise ValueError(f'{absent} not in {kind["stats"]}')
    days = sorted(samp['day'].unique())
    args = (kind, samp, st, days)

    out = os.path.join(FIG_DIR, kind['out'].format('cells'))
    with PdfPages(out) as pdf:
        for grp in GROUPS:
            draw_page(pdf, *args, kind['page_cells'].format(grp), [kind['row'].format(k) for k in ARCHETYPES], MODULES,
                      lambda rl, cl: (grp, rl.split()[0], cl))
    print(f'wrote {out}')

    out = os.path.join(FIG_DIR, kind['out'].format('module'))
    with PdfPages(out) as pdf:
        for k in ARCHETYPES:
            draw_page(pdf, *args, kind['page_module'].format(k), [f'{k}-{m}' for m in MODULES], GROUPS,
                      lambda rl, cl: (cl, k, rl.split('-')[1]))
    print(f'wrote {out}')
