"""Figures of 2_development/06: Yoo25 L2/3 cells in P21 PC space colored by archetype score, and the
variance of each Age along the A-C and B axes, NR vs DR.

  08.yoo25_l23_pc12_abc_scores.pdf      one column per Age (NR P6-P21, then DR P12-P21); cells in the
                                        P21 (PC1, PC2) with the P21 archetype triangle (dashed). Row 1
                                        colored by score C - A, row 2 by score B (2_development/01); color
                                        range = 5th-95th percentile over all cells.
  08.yoo25_l23_ac_b_variance_over_age.pdf  variance over postnatal day of P8 PC1, the A-C axis and the B axis;
                                        lines = per Age, markers = per Sample, dashed = shuffled-gene null.
                                        P8 PC1 and A-C are divided by the max A-C variance over Ages,
                                        B by the max B variance.

Reads:
  results/2_development/06.yoo25_l23_hvg_pcs_cells.tsv
  results/2_development/06.yoo25_l23_hvg_pcs_archetypes.tsv
  results/2_development/06.yoo25_l23_hvg_pcs_variance.tsv
  results/2_development/01.yoo25_l23_abc_scores.tsv
Outputs:
  figures/2_development/09.yoo25_l23_pc12_abc_scores.pdf
  figures/2_development/09.yoo25_l23_ac_b_variance_over_age.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.ticker import MaxNLocator

from l23.constants import AGES, COND_COLORS
from l23.utils import plot_archetype

matplotlib.rcParams['pdf.fonttype'] = 42

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_CELLS  = os.path.join(REPO_ROOT, 'results', '2_development', '06.yoo25_l23_hvg_pcs_cells.tsv')
INPUT_ARCH   = os.path.join(REPO_ROOT, 'results', '2_development', '06.yoo25_l23_hvg_pcs_archetypes.tsv')
INPUT_VAR    = os.path.join(REPO_ROOT, 'results', '2_development', '06.yoo25_l23_hvg_pcs_variance.tsv')
INPUT_SCORES = os.path.join(REPO_ROOT, 'results', '2_development', '01.yoo25_l23_abc_scores.tsv')
OUT_PC12     = os.path.join(REPO_ROOT, 'figures', '2_development', '09.yoo25_l23_pc12_abc_scores.pdf')
OUT_VAR      = os.path.join(REPO_ROOT, 'figures', '2_development', '09.yoo25_l23_ac_b_variance_over_age.pdf')

PCTILES   = (5, 95)      # color range of the scores
POINT_S   = 4
FS_TITLE  = 10
FS_LABEL  = 9
FS_TICK   = 8
NULL_COLORS = {'NR': 'gray', 'DR': 'lightgray'}
MARKERS     = {'NR': 'o', 'DR': 's'}


def mix(c1, c2, w=0.5):
    return tuple(w * np.array(to_rgb(c1)) + (1 - w) * np.array(to_rgb(c2)))


def ac_cmap():
    """A (C0) -> gray -> C (C2), for score C - A."""
    return LinearSegmentedColormap.from_list('ac', ['C0', mix('C0', 'lightgray'), 'lightgray',
                                                    mix('C2', 'lightgray'), 'C2'])


def b_cmap():
    """gray -> B (C1), for score B."""
    return LinearSegmentedColormap.from_list('b', ['lightgray', mix('C1', 'white', 0.2), 'C1'])


def score_row(fig, axs, cells, aa, col, cmap, label):
    """One row of PC1/PC2 panels, one per Age, colored by `col`, with a colorbar at the row's end."""
    vmin, vmax = np.percentile(cells[col], PCTILES)
    for ax, age in zip(axs, AGES):
        c = cells[cells['Age'] == age]
        sc = ax.scatter(c['pc_later1'], c['pc_later2'], c=c[col], cmap=cmap, vmin=vmin, vmax=vmax,
                        s=POINT_S, edgecolor='none', rasterized=True)
        plot_archetype(ax, aa, fmt='--', color='gray', lw=0.8, zorder=2)
        ax.set_aspect('equal')
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.tick_params(labelsize=FS_TICK)
        ax.spines[['top', 'right']].set_visible(False)
    cb = fig.colorbar(sc, ax=axs, fraction=0.01, pad=0.01, aspect=15)
    cb.set_label(label, fontsize=FS_LABEL)
    cb.ax.tick_params(labelsize=FS_TICK)


def variance_panel(ax, var, col, scale, null_col=None, legend=False):
    """Variance of `col` over postnatal day, NR and DR: lines per Age, open markers per Sample."""
    for cond in ['NR', 'DR']:
        a = var[(var['level'] == 'Age') & (var['condition'] == cond)]
        s = var[(var['level'] == 'Sample') & (var['condition'] == cond)]
        ax.plot(a['day'], a[col] / scale, color=COND_COLORS[cond], label=cond)
        ax.plot(s['day'], s[col] / scale, MARKERS[cond], ms=5, fillstyle='none', color=COND_COLORS[cond])
        if null_col:
            ax.plot(a['day'], a[null_col] / scale, '--', color=NULL_COLORS[cond], label=f'shuffled {cond}')
            ax.plot(s['day'], s[null_col] / scale, MARKERS[cond], ms=5, fillstyle='none', color=NULL_COLORS[cond])
    ax.set_xticks(sorted({d for d, _ in AGES.values()}))
    ax.set_xlabel('Postnatal day (P)', fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.spines[['top', 'right']].set_visible(False)
    if legend:
        ax.legend(fontsize=FS_TICK, frameon=False, bbox_to_anchor=(1, 1), loc='upper left')


cells = pd.read_csv(INPUT_CELLS, sep='\t', index_col='cell')
scores = pd.read_csv(INPUT_SCORES, sep='\t', index_col='cell')
if not scores.index.sort_values().equals(cells.index.sort_values()):
    raise ValueError(f'{INPUT_SCORES} cells differ from {INPUT_CELLS}')
scores = scores.reindex(cells.index)
cells['score_c_minus_a'] = scores['score_c'] - scores['score_a']
cells['score_b'] = scores['score_b']
arch = pd.read_csv(INPUT_ARCH, sep='\t').set_index('archetype').loc[['A', 'B', 'C']]
aa = arch[['pc_later1', 'pc_later2']].to_numpy().T
var = pd.read_csv(INPUT_VAR, sep='\t')
os.makedirs(os.path.dirname(OUT_PC12), exist_ok=True)

# cells in P21 PC space, colored by archetype score
n = len(AGES)
fig, axs = plt.subplots(2, n, figsize=(1.75 * n, 3.7), sharex=True, sharey=True, layout='constrained')
score_row(fig, axs[0], cells, aa, 'score_c_minus_a', ac_cmap(), 'score C − A')
score_row(fig, axs[1], cells, aa, 'score_b', b_cmap(), 'score B')
for ax, age in zip(axs[0], AGES):
    ax.set_title(age, fontsize=FS_TITLE)
axs[1, 0].set_xlabel('PC1 (P21)', fontsize=FS_LABEL)
for ax in axs[:, 0]:
    ax.set_ylabel('PC2 (P21)', fontsize=FS_LABEL)
fig.savefig(OUT_PC12, dpi=300)
print(f'wrote {OUT_PC12}')

# variance along the axes over age
scale_ac = var.loc[var['level'] == 'Age', 'ac'].max()
scale_b = var.loc[var['level'] == 'Age', 'b'].max()
fig, axs = plt.subplots(1, 3, figsize=(9.5, 3), sharex=True, sharey=True, layout='constrained')
variance_panel(axs[0], var, 'pc_early1', scale_ac)
variance_panel(axs[1], var, 'ac', scale_ac, null_col='ac_shuff')
variance_panel(axs[2], var, 'b', scale_b, null_col='b_shuff', legend=True)
for ax, title in zip(axs, ['P8 PC1 (early A-C axis)', 'A-C axis', 'B axis']):
    ax.set_title(title, fontsize=FS_TITLE)
axs[0].set_ylabel('norm. variance', fontsize=FS_LABEL)
axs[0].set_ylim(bottom=0)
fig.savefig(OUT_VAR)
print(f'wrote {OUT_VAR}')
