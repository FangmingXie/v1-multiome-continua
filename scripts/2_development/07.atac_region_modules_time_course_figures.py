"""Figures of 2_development/06: A/B/C chromatin region modules (M1-M5) over postnatal day in Yoo25 L2/3 ATAC, NR vs DR,
per cell group (all L2/3, A', B', C' cells).

One panel per (cell group, region set, module): y = module score per sample (mean over the module's
regions of log2(CPM in peaks + 1) minus the region's P21 NR mean, all L2/3; dashed line at 0); lines
through the mean across samples per age x condition, dots = samples. Title: n = module regions passing
the filter; DR and slope = OLS estimate (padj, BH over the 15 modules of the cell group) of 2_development/06. Marked
(bracket right of P21): padj < FDR, no effect-size cutoff. y is shared within a row. The same panels in
two layouts:
  OUT_FIG_BY_CELLS   one page per cell group; rows = A/B/C region sets, columns = M1-M5
  OUT_FIG_BY_MODULE  one page per region set (A, B, C); rows = M1-M5, columns = cell groups

Reads:
  results/2_development/06.yoo25_l23_atac_region_modules_samples.tsv
  results/2_development/06.yoo25_l23_atac_region_modules_stats.tsv
Outputs:
  figures/2_development/07.yoo25_l23_atac_region_modules_time_course_by_cells.pdf
  figures/2_development/07.yoo25_l23_atac_region_modules_time_course_by_module.pdf
"""

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_SAMPLES     = os.path.join(REPO_ROOT, 'results', '2_development', '06.yoo25_l23_atac_region_modules_samples.tsv')
INPUT_STATS       = os.path.join(REPO_ROOT, 'results', '2_development', '06.yoo25_l23_atac_region_modules_stats.tsv')
OUT_FIG_BY_CELLS  = os.path.join(REPO_ROOT, 'figures', '2_development', '07.yoo25_l23_atac_region_modules_time_course_by_cells.pdf')
OUT_FIG_BY_MODULE = os.path.join(REPO_ROOT, 'figures', '2_development', '07.yoo25_l23_atac_region_modules_time_course_by_module.pdf')

GROUPS      = ['all L2/3', "A' cells", "B' cells", "C' cells"]
ARCHETYPES  = ['A', 'B', 'C']
MODULES     = ['M1', 'M2', 'M3', 'M4', 'M5']
COND_COLORS = {'NR': 'C1', 'DR': 'black'}
FDR         = 0.05
STARS       = [(0.001, '***'), (0.01, '**'), (0.05, '*')]
PANEL_W, PANEL_H = 2.5, 3.0   # inches

os.makedirs(os.path.dirname(OUT_FIG_BY_CELLS), exist_ok=True)

samp = pd.read_csv(INPUT_SAMPLES, sep='\t')
st = pd.read_csv(INPUT_STATS, sep='\t').set_index(['cells', 'archetype', 'module'])
expected = pd.MultiIndex.from_product([GROUPS, ARCHETYPES, MODULES])
absent = [x for x in expected if x not in st.index]
if absent:
    raise ValueError(f'{absent} not in {INPUT_STATS}')
days = sorted(samp['day'].unique())


def stars(padj):
    return next((s for cut, s in STARS if padj < cut), '')


def draw_panel(ax, grp, k, m):
    ax.axhline(0, color='gray', lw=0.8, ls='--')
    ends = {}
    for c, color in COND_COLORS.items():
        s = samp[(samp['cells'] == grp) & (samp['archetype'] == k) & (samp['module'] == m) & (samp['condition'] == c)]
        t = s.groupby('day')['score'].mean()
        ax.plot(t.index, t.values, color=color, lw=1.5, label=c)
        ax.scatter(s['day'], s['score'], color=color, s=10, zorder=3)
        ends[c] = t.loc[t.index.max()]
    r = st.loc[(grp, k, m)]
    lab = [f'{name} {stars(r[f"padj_{key}"])}' for name, key in [('DR', 'DR_d16'), ('slope', 'slope')]
           if stars(r[f'padj_{key}'])]
    if lab:
        # bracket: short ticks from each line end out to a vertical bar
        ax.plot([21.8, 22.1, 22.1, 21.8], [ends['NR'], ends['NR'], ends['DR'], ends['DR']], color='black', lw=1)
        ax.text(22.4, (ends['NR'] + ends['DR']) / 2, '\n'.join(lab), fontsize=7, va='center', ha='left')
    ax.set_xlim(5, 25.5)
    ax.set_xticks(days)
    ax.tick_params(axis='x', labelsize=6.5)
    ax.set_title(f'n = {int(r["n_regions"])} regions\n'
                 f'DR {r["DR_d16"]:.2f} ({r["padj_DR_d16"]:.2g})\n'
                 f'slope {r["slope_per_day"]:+.3f} ({r["padj_slope"]:.2g})', fontsize=6.5)
    ax.spines[['top', 'right']].set_visible(False)


def draw_page(pdf, title, rows, cols, key):
    """One page: rows x cols panels; key(row, col) -> (cell group, region set, module)."""
    fig, axes = plt.subplots(len(rows), len(cols), figsize=(PANEL_W * len(cols), PANEL_H * len(rows)),
                             sharey='row', squeeze=False)
    for i, rl in enumerate(rows):
        for j, cl in enumerate(cols):
            draw_panel(axes[i, j], *key(rl, cl))
            if i == 0:
                axes[i, j].text(0.5, 1.3, cl, transform=axes[i, j].transAxes, ha='center', va='bottom', fontsize=10)
        axes[i, 0].set_ylabel(rl, fontsize=11)
    for ax in axes[-1]:
        ax.set_xlabel('postnatal day')
    axes[0, 0].legend(fontsize=8, frameon=False, loc='lower left')
    fig.supylabel('mean log$_2$(CPM in peaks + 1) − P21 NR mean')
    fig.suptitle(f'{title}, NR vs DR; reference = P21 NR mean, all L2/3 (dashed line)\n'
                 'OLS of the module score over P12-P21, ~ age + DR + DR·(day - 16)\n'
                 'DR = average DR effect (at day 16); slope = change in DR effect per day. '
                 'Panel: estimate (padj, BH over 15 modules)\n'
                 f'marked (bracket): padj < {FDR}; * < 0.05, ** < 0.01, *** < 0.001\n'
                 'n = module regions passing the filter (>= 10 counts in >= 2 samples)',
                 fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.9 / (PANEL_H * len(rows))))
    pdf.savefig(fig)
    plt.close(fig)


with PdfPages(OUT_FIG_BY_CELLS) as pdf:
    for grp in GROUPS:
        draw_page(pdf, f'{grp}: A/B/C ATAC region modules (M1-M5)', [f'{k} regions' for k in ARCHETYPES], MODULES,
                  lambda rl, cl: (grp, rl.split()[0], cl))
print(f'wrote {OUT_FIG_BY_CELLS}')

with PdfPages(OUT_FIG_BY_MODULE) as pdf:
    for k in ARCHETYPES:
        draw_page(pdf, f'{k} ATAC region modules (M1-M5) by cell group', [f'{k}-{m}' for m in MODULES], GROUPS,
                  lambda rl, cl: (cl, k, rl.split('-')[1]))
print(f'wrote {OUT_FIG_BY_MODULE}')
