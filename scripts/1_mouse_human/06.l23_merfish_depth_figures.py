"""Figures of 04: conserved L2/3 A/C genes in mouse vs human MERFISH, on one relative L2/3 depth scale.

  OUT_FIG_PROFILE  per gene (rows = A, C side): mean log expression per relative-depth bin divided by the
                   section's peak bin, mean +- SEM over sections, mouse vs human.
  OUT_FIG_ALLSEC   depth-alignment check: the in situ genes (columns) in every section (rows; mouse, then
                   human); all panels the same size, y = relative depth over ALLSEC_YLIM, x = position
                   stretched to the panel; dashed lines at the top and bottom of L2/3 (0, 1); colour = log
                   expression / the gene's GENE_PCTILE percentile in that section, capped at 1; % right of
                   each panel = cells with the gene detected.
  OUT_FIG_CURVE    every section: left, the raw L2/3 cells (um) with the points and 2nd-order fit of the top of
                   L2/3; right, the flattened cells without depth outliers on the relative scale.

Reads:
  results/1_mouse_human/04.l23_merfish_{genes,sections,top_fit_points,rel_depth_profiles}.tsv
  results/1_mouse_human/04.l23_merfish_cells.parquet
Outputs:
  figures/1_mouse_human/06.merfish_gene_depth_profiles.pdf
  figures/1_mouse_human/06.merfish_gene_in_situ_maps.pdf
  figures/1_mouse_human/06.merfish_l23_top_fit.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RES_DIR      = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
FIG_DIR      = os.path.join(REPO_ROOT, 'figures', '1_mouse_human')
IN_GENES     = os.path.join(RES_DIR, '04.l23_merfish_genes.tsv')
IN_SECTIONS  = os.path.join(RES_DIR, '04.l23_merfish_sections.tsv')
IN_POINTS    = os.path.join(RES_DIR, '04.l23_merfish_top_fit_points.tsv')
IN_CELLS     = os.path.join(RES_DIR, '04.l23_merfish_cells.parquet')
IN_PROFILES  = os.path.join(RES_DIR, '04.l23_merfish_rel_depth_profiles.tsv')
OUT_FIG_PROFILE = os.path.join(FIG_DIR, '06.merfish_gene_depth_profiles.pdf')
OUT_FIG_ALLSEC  = os.path.join(FIG_DIR, '06.merfish_gene_in_situ_maps.pdf')
OUT_FIG_CURVE   = os.path.join(FIG_DIR, '06.merfish_l23_top_fit.pdf')

ALLSEC_YLIM = (-0.35, 1.35)   # all-sections figure: same relative-depth window in every panel
GENE_PCTILE = 98              # all-sections figure: per gene and section, log expression / this percentile, capped at 1
GENE_CMAP   = sns.color_palette('rocket_r', as_cmap=True)
SPECIES_COLORS = {'mouse': 'C3', 'human': 'k'}

os.makedirs(FIG_DIR, exist_ok=True)

genes = pd.read_csv(IN_GENES, sep='\t')
gi = genes.set_index('mouse_symbol', drop=False)
insitu = list(genes.loc[genes['insitu'], 'mouse_symbol'])
sections = pd.read_csv(IN_SECTIONS, sep='\t')
cells = pd.read_parquet(IN_CELLS)
points = pd.read_csv(IN_POINTS, sep='\t')
prof = pd.read_csv(IN_PROFILES, sep='\t')
rows = [(sp, sec, cells[(cells['section'] == sec) & ~cells['outlier']])
        for sp, sec in zip(sections['species'], sections['section'])]
n_sec = sections['species'].value_counts()

# --- depth profiles on the relative L2/3 scale, mean +- SEM over sections ---
n_col = genes.groupby('side').size().max()
sides = list(dict.fromkeys(genes['side']))
fig, axes = plt.subplots(len(sides), n_col, figsize=(2.6 * n_col, 3.4 * len(sides)), squeeze=False)
gax = []
for i, side in enumerate(sides):
    gl = genes.loc[genes['side'] == side, 'mouse_symbol']
    gax += list(zip(axes[i, :len(gl)], gl))
    for ax in axes[i, len(gl):]:
        ax.set_visible(False)
for ax, mg in gax:
    for sp in ['mouse', 'human']:
        t = prof[(prof['gene'] == mg) & (prof['species'] == sp)]
        per = t.pivot(index='bin_mid', columns='section', values='rel_expr')[list(dict.fromkeys(t['section']))]
        ax.errorbar(per.mean(axis=1), per.index, xerr=per.sem(axis=1), color=SPECIES_COLORS[sp], lw=1.8, marker='o',
                    ms=3.5, capsize=2, elinewidth=1, label=f'{sp} (n = {per.shape[1]} sections)')
    side = gi.loc[mg, 'side']
    ax.set_title(f"{gi.loc[mg, 'human_symbol']} / {mg} ({side} side)", fontsize=9, color='C0' if side == 'A' else 'C2')
    ax.set_xlabel('mean expression / peak bin', fontsize=8)
    ax.set_xlim(0, 1.1)
    ax.set_ylim(1, 0)   # top of L2/3 on top
    ax.spines[['top', 'right']].set_visible(False)
gax[0][0].legend(fontsize=7, frameon=False, loc='lower left')
for ax in axes[:, 0]:
    ax.set_ylabel('relative L2/3 depth (0 = top, 1 = bottom)', fontsize=8)
fig.suptitle(f"L2/3 depth profiles on one relative scale: {n_sec['mouse']} mouse NR vs "
             f"{n_sec['human']} human sections\n(human: {' and '.join(gi.loc[insitu, 'human_symbol'])} each detected in > "
             f"{sections['detect_threshold'].iloc[0]:.0%} of L2/3 cells)\nper section: bin means divided by the section's "
             f"peak bin; mean ± SEM over sections; {prof['bin_mid'].nunique()} bins", fontsize=9)
fig.tight_layout()
fig.savefig(OUT_FIG_PROFILE)
print(f'wrote {OUT_FIG_PROFILE}')

# --- in situ genes in every section, to inspect depth alignment ---
fig, axes = plt.subplots(len(rows), len(insitu), figsize=(4.2 * len(insitu) + 1.6, 1.25 * len(rows) + 0.8),
                         squeeze=False)
fig.subplots_adjust(left=0.24, right=0.86, top=1 - 0.9 / (1.25 * len(rows) + 0.8), bottom=0.02, hspace=0.35, wspace=0.16)
for i, (sp, sec, c) in enumerate(rows):
    for j, mg in enumerate(insitu):
        ax = axes[i, j]
        v = c[mg]
        scale = np.percentile(v, GENE_PCTILE)
        if scale <= 0:
            scale = v.max() if v.max() > 0 else 1.0   # gene detected in < 2% of cells
        o = v.sort_values().index   # zero counts = lowest expression (0); high expressers drawn on top
        sc_ = ax.scatter(c.loc[o, 'rel_x'], c.loc[o, 'rel_depth'], s=2.5, c=np.clip(v[o] / scale, 0, 1),
                         cmap=GENE_CMAP, vmin=0, vmax=1, lw=0, rasterized=True)
        for yl in (0, 1):
            ax.axhline(yl, color='k', lw=0.6, ls='--')
        ax.set_ylim(ALLSEC_YLIM[1], ALLSEC_YLIM[0])   # top of L2/3 on top
        ax.axis('off')
        if i == 0:
            ax.set_title(f"{gi.loc[mg, 'human_symbol']} / {mg} ({gi.loc[mg, 'side']} side)", fontsize=9)
        ax.text(1.01, 0.5, f'{(v > 0).mean():.0%}', transform=ax.transAxes, ha='left', va='center', fontsize=6)
    axes[i, 0].text(-0.03, 0.5, f'{sp}\n{sec}\n(n = {len(c)})', transform=axes[i, 0].transAxes,
                    ha='right', va='center', fontsize=7)
fig.colorbar(sc_, cax=fig.add_axes([0.9, 0.45, 0.012, 0.1]), label=f'log expr / {GENE_PCTILE}th pct')
fig.suptitle('Depth alignment check: L2/3 cells of every profile section on the relative scale\n'
             '(dashed = top / bottom of L2/3; % right of each panel = cells with the gene detected)', fontsize=8)
fig.savefig(OUT_FIG_ALLSEC)
print(f'wrote {OUT_FIG_ALLSEC}')

# --- every section: fitted top of L2/3 (left: raw, right: flattened, relative scale) ---
fig, axes = plt.subplots(len(rows), 2, figsize=(11, 2.6 * len(rows)), squeeze=False)
for i, (sp, sec, c) in enumerate(rows):
    raw = cells[cells['section'] == sec]
    pts = points[points['section'] == sec]
    top = np.poly1d(sections.set_index('section').loc[sec, ['top_c2', 'top_c1', 'top_c0']].to_numpy(float))
    xs = np.linspace(raw['x'].min(), raw['x'].max(), 200)
    ax = axes[i, 0]
    ax.scatter(raw['x'], raw['depth'], s=2, c='0.6', lw=0, rasterized=True)
    ax.scatter(pts['x'], pts['depth'], s=14, c='C3', zorder=3, label=f'{pts["depth_pctile"].iloc[0]}th pct depth per bin')
    ax.plot(xs, top(xs), color='C3', lw=1.2, label='2nd-order fit (top of L2/3)')
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_title(f'{sp} {sec}: raw L2/3 cells (µm)', fontsize=8, loc='left')
    ax = axes[i, 1]
    ax.scatter(c['rel_x'], c['rel_depth'], s=2, c='0.6', lw=0, rasterized=True)
    for yl in (0, 1):
        ax.axhline(yl, color='k', lw=0.6, ls='--')
    ax.set_aspect('equal')
    ax.invert_yaxis()
    ax.set_title('flattened, outliers removed, relative scale', fontsize=8, loc='left')
axes[0, 0].legend(fontsize=7, frameon=False, loc='lower left')
fig.tight_layout()
fig.savefig(OUT_FIG_CURVE)
print(f'wrote {OUT_FIG_CURVE}')
