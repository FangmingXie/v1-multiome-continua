"""Figure of 02: archetype scores along the A'-C' axis (snRNA) and depth (MERFISH), next to the EASI-FISH profiles.

Panels left to right: Yoo25 P21 NR snRNA along the normalized A'-C' axis; P21 NR MERFISH along normalized
depth; EASI-FISH LM vs RL retrograde labels along normalized depth. Score panels: archetype score
z-scored within the dataset, mean +/- SEM per equal-count bin (dots on the snRNA panel); EASI-FISH:
fraction of cells per group, mean +/- SD across samples. y: 0 = A' end / superficial (top), 1 = C' end / deep.

Reads:
  results/5_easifish/02.l23_archetype_scores_and_easifish_proj_depth_{bins,easifish_bins,panels}.tsv
Outputs:
  figures/5_easifish/03.l23_archetype_scores_and_easifish_proj_depth.pdf
"""

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from l23.constants import ARCH_COLORS

matplotlib.rcParams['pdf.fonttype'] = 42

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

IN_PREFIX = os.path.join(REPO_ROOT, 'results', '5_easifish', '02.l23_archetype_scores_and_easifish_proj_depth')
OUT_FIG   = os.path.join(REPO_ROOT, 'figures', '5_easifish', '03.l23_archetype_scores_and_easifish_proj_depth.pdf')

PROJ_COLORS = {'LM': 'magenta', 'RL': 'green'}

os.makedirs(os.path.dirname(OUT_FIG), exist_ok=True)

bins = pd.read_csv(f'{IN_PREFIX}_bins.tsv', sep='\t')
easi = pd.read_csv(f'{IN_PREFIX}_easifish_bins.tsv', sep='\t')
pan = pd.read_csv(f'{IN_PREFIX}_panels.tsv', sep='\t').set_index('panel')
score_panels = list(dict.fromkeys(bins['panel']))   # snRNA first, then the MERFISH datasets
n_bins = bins.groupby(['panel', 'archetype']).size().iloc[0]
cond = pan.loc[score_panels[0], 'marker_condition']
y_pct = pan.loc[score_panels[0], ['y_pctile_lo', 'y_pctile_hi']].astype(int).tolist()


def genes_txt(p):
    return ', '.join(f'{k} {int(pan.loc[p, f"n_genes_{k}"])}' for k in ARCH_COLORS)


fig, axes = plt.subplots(1, 1 + len(score_panels), figsize=(3.5 * (1 + len(score_panels)), 4.7))
ax1, ax_easi = axes[0], axes[-1]
for ax, p in zip(axes, score_panels):
    for k in ARCH_COLORS:
        t = bins[(bins['panel'] == p) & (bins['archetype'] == k)]
        ax.fill_betweenx(t['y_norm'], t['z_mean'] - t['z_sem'], t['z_mean'] + t['z_sem'],
                         color=ARCH_COLORS[k], alpha=0.3, lw=0)
        ax.plot(t['z_mean'], t['y_norm'], color=ARCH_COLORS[k], lw=1, label=f'{k} score',
                marker='o' if ax is ax1 else None, ms=2)   # dots on panel 1 only
    ax.set_xlabel('archetype score (z, within dataset)')
    ax.spines[['top', 'right']].set_visible(False)
ax1.set_ylabel("C'  ←  normalized A'-C' axis  →  A'")   # reads bottom to top
r = pan.loc[score_panels[0]]
ax1.set_title(f'{r["label"]}, L2/3 (n = {int(r["n_cells"])})\n'
              + genes_txt(score_panels[0]) + f' genes; axis {r["y0"]:.2f} to {r["y1"]:.2f}', fontsize=8)
for ax, p in zip(axes[1:-1], score_panels[1:]):
    r = pan.loc[p]
    ax.set_title(f'{r["label"]}, L2/3 (n = {int(r["n_cells"])})\n'
                 + genes_txt(p) + f' panel genes; {r["y0"]:.0f}-{r["y1"]:.0f} µm', fontsize=8)
for g, color in PROJ_COLORS.items():
    t = easi[easi['group'] == g]
    ax_easi.fill_betweenx(t['y_norm'], t['mean'] - t['sd'], t['mean'] + t['sd'], color=color, alpha=0.3, lw=0)
    ax_easi.plot(t['mean'], t['y_norm'], color=color, lw=1,
                 label=f'{g}-projecting\n(n = {pan.loc[f"easifish_{g}", "n_cells_per_sample"]})')
ax_easi.set_xlabel('fraction of cells per group')
r = pan.loc['easifish_LM']
ax_easi.set_title(f'EASI-FISH retrograde labels, mean ± SD\nacross 3 samples; {r["y0"]:.0f}-{r["y1"]:.0f} µm', fontsize=8)
ax_easi.legend(fontsize=7, frameon=False, loc='upper left', bbox_to_anchor=(1.0, 1.0))
ax_easi.spines[['top', 'right']].set_visible(False)
for ax in axes[1:]:
    ax.set_ylabel('deep  ←  normalized depth  →  superficial')
for ax in axes:   # all panels: normalized 0 (A' end / superficial, top) to 1 (C' end / deep)
    ax.set_ylim(1, 0)
lo, hi = y_pct
fig.suptitle(f'{cond} archetype scores (panels 1-2: mean ± SEM in {n_bins} equal-count bins) '
             'and EASI-FISH projection depth profiles (panel 3)'
             f"\n(y: panels 1-2, {lo}st-{hi}th cell percentiles mapped to 0-1; "
             f"panel 3, {r['y0']:.0f}-{r['y1']:.0f} µm mapped to 0-1; "
             "panel 1: A'/C' scores partly define the axis)", fontsize=9)
fig.legend(*ax1.get_legend_handles_labels(), fontsize=8, frameon=False, ncol=3,
           loc='lower center')
fig.tight_layout(rect=(0, 0.05, 1, 1))
fig.savefig(OUT_FIG)
print(f'wrote {OUT_FIG}')
