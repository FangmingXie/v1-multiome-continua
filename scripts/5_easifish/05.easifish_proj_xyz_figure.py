"""Figures of 04: EASI-FISH LM vs RL retrograde labels in 3 samples (scomics revision/e00).

  05.easifish_cells_xyz.pdf        rows = samples; x-y, x-z, y-z views of every cell (gray) with LM, RL and
                                   both cells on top (n in the legend); depth (z) down.
  05.easifish_depth_profiles.pdf   rows: cells per 20 um depth bin (LM, RL, both); fraction of the group's
                                   binned cells per bin (LM, RL); columns = samples; values at bin centers;
                                   dashed lines = L2/3 bounds (L23_BOUNDS).
  05.easifish_xy_proj_bias.pdf     rows: log2((#RL + 1) / (#LM + 1)) per 100 um x-y bin; mean depth of RL minus
                                   LM cells per bin (blank where either is absent); columns = samples.
  05.easifish_xy_dz_null.pdf       per sample and pooled: density of the per-bin mean-depth difference
                                   (RL - LM), data vs the label-shuffle null of 04.
Display orientation per sample as e00 (ROTATE: sample 1 is shown rotated by 180 degrees in x-y); display only.

Reads:
  data/raw/easifish/easifish_sample{1,2,3}_res.csv   (cell positions and labels, for the dot maps)
  results/5_easifish/04.easifish_depth_hist.tsv
  results/5_easifish/04.easifish_xy_bins.tsv
  results/5_easifish/04.easifish_xy_dz_null.tsv
Outputs:
  figures/5_easifish/05.easifish_cells_xyz.pdf
  figures/5_easifish/05.easifish_depth_profiles.pdf
  figures/5_easifish/05.easifish_xy_proj_bias.pdf
  figures/5_easifish/05.easifish_xy_dz_null.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SAMPLES        = (1, 2, 3)
INPUT_EASIFISH = {s: os.path.join(REPO_ROOT, 'data', 'raw', 'easifish', f'easifish_sample{s}_res.csv')
                  for s in SAMPLES}
RES_DIR        = os.path.join(REPO_ROOT, 'results', '5_easifish')
IN_DEPTH       = os.path.join(RES_DIR, '04.easifish_depth_hist.tsv')
IN_XY          = os.path.join(RES_DIR, '04.easifish_xy_bins.tsv')
IN_NULL        = os.path.join(RES_DIR, '04.easifish_xy_dz_null.tsv')
FIG_DIR        = os.path.join(REPO_ROOT, 'figures', '5_easifish')
OUT_CELLS      = os.path.join(FIG_DIR, '05.easifish_cells_xyz.pdf')
OUT_DEPTH      = os.path.join(FIG_DIR, '05.easifish_depth_profiles.pdf')
OUT_XY         = os.path.join(FIG_DIR, '05.easifish_xy_proj_bias.pdf')
OUT_NULL       = os.path.join(FIG_DIR, '05.easifish_xy_dz_null.pdf')

GROUP_COLORS = {'LM': 'magenta', 'RL': 'green', 'both': 'orange'}
ROTATE       = {1: True, 2: False, 3: False}   # x-y view rotated by 180 degrees (as e00)
VIEWS        = [('x', 'y'), ('x', 'z'), ('y', 'z')]
L23_BOUNDS   = (100, 400)                       # um
LOG2_LIM     = 3
DZ_LIM       = 100                              # um
NULL_BINS    = {'sample': np.linspace(-150, 150, 50), 'pooled': np.linspace(-150, 150, 100)}
FS_TITLE     = 10
FS_LABEL     = 9
FS_TICK      = 8
CMAP_LOG2    = LinearSegmentedColormap.from_list('rl_lm', ['magenta', 'white', 'green'])
CMAP_DZ      = LinearSegmentedColormap.from_list('dz', ['black', 'white', 'C1'])

os.makedirs(FIG_DIR, exist_ok=True)
plt.rcParams['pdf.fonttype'] = 42


def style(ax):
    ax.tick_params(labelsize=FS_TICK)
    ax.spines[['top', 'right']].set_visible(False)


def orient(ax, sample, view):
    """Depth down; the e00 per-sample rotation (x-y by 180 degrees, i.e. both axes reversed)."""
    if view[1] == 'z':
        ax.invert_yaxis()
    if ROTATE[sample]:
        ax.invert_xaxis()
        if view[1] != 'z':
            ax.invert_yaxis()


def xy_map(fig, axs, xy, col, cmap, lim, label):
    """One row of x-y heatmaps (one per sample) of `col`, shared color range, one colorbar."""
    for ax, s in zip(axs, SAMPLES):
        t = xy[xy['sample'] == s]
        grid = t.pivot(index='ybin', columns='xbin', values=col)
        xe = t.drop_duplicates('xbin').sort_values('xbin')['x_lo'].to_numpy()
        ye = t.drop_duplicates('ybin').sort_values('ybin')['y_lo'].to_numpy()
        xe, ye = np.r_[xe, xe[-1] + np.diff(xe[:2])], np.r_[ye, ye[-1] + np.diff(ye[:2])]
        m = ax.pcolormesh(xe, ye, np.ma.masked_invalid(grid.to_numpy()), cmap=cmap, vmin=-lim, vmax=lim)
        ax.set_aspect('equal')
        ax.set_title(f'sample {s}', fontsize=FS_TITLE)
        ax.set_xlabel('x (um)', fontsize=FS_LABEL)
        style(ax)
        orient(ax, s, ('x', 'y'))
    axs[0].set_ylabel('y (um)', fontsize=FS_LABEL)
    cb = fig.colorbar(m, ax=axs, shrink=0.8)
    cb.set_label(label, fontsize=FS_LABEL)
    cb.ax.tick_params(labelsize=FS_TICK)


cells = {s: pd.read_csv(INPUT_EASIFISH[s]) for s in SAMPLES}
depth = pd.read_csv(IN_DEPTH, sep='\t')
xy = pd.read_csv(IN_XY, sep='\t')
null = pd.read_csv(IN_NULL, sep='\t')

# cells in x-y, x-z, y-z
fig, axs = plt.subplots(len(SAMPLES), len(VIEWS), figsize=(15, 4.2 * len(SAMPLES)), layout='constrained')
for row, s in zip(axs, SAMPLES):
    df = cells[s]
    groups = {'LM': df['LM'] == 1, 'RL': df['RL'] == 1, 'both': (df['LM'] == 1) & (df['RL'] == 1)}
    for ax, view in zip(row, VIEWS):
        ax.scatter(df[view[0]], df[view[1]], c='lightgray', s=1, linewidths=0, rasterized=True)
        for g, m in groups.items():
            ax.scatter(df.loc[m, view[0]], df.loc[m, view[1]], c=GROUP_COLORS[g], s=4, linewidths=0,
                       rasterized=True, label=f'{g} (n={m.sum():,})')
        ax.set_aspect('equal')
        ax.set_xlabel(f'{view[0]} (um)', fontsize=FS_LABEL)
        ax.set_ylabel(f'{view[1]} (um)', fontsize=FS_LABEL)
        style(ax)
        orient(ax, s, view)
    row[0].set_title(f'sample {s}', fontsize=FS_TITLE, loc='left')
    row[0].legend(fontsize=FS_TICK, markerscale=3, frameon=False, loc='upper left', bbox_to_anchor=(1, 1))
fig.savefig(OUT_CELLS, dpi=300)
print(f'wrote {OUT_CELLS}')

# depth profiles
fig, axs = plt.subplots(2, len(SAMPLES), figsize=(4 * len(SAMPLES), 8), sharey=True, layout='constrained')
for j, s in enumerate(SAMPLES):
    d = depth[depth['sample'] == s]
    for i, (col, groups, xlabel) in enumerate([('n', ['LM', 'RL', 'both'], 'number of cells'),
                                               ('frac', ['LM', 'RL'], 'fraction of cells')]):
        ax = axs[i, j]
        for g in groups:
            t = d[d['group'] == g]
            ax.plot(t[col], (t['z_lo'] + t['z_hi']) / 2, color=GROUP_COLORS[g], label=f'{g} (n={t["n"].sum():,})')
        for b in L23_BOUNDS:
            ax.axhline(b, color='gray', linestyle='--', lw=0.8)
        ax.set_xlabel(xlabel, fontsize=FS_LABEL)
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()
        ax.tick_params(labelsize=FS_TICK)
        ax.spines[['bottom', 'right']].set_visible(False)
        ax.legend(fontsize=FS_TICK, frameon=False, loc='lower right')
    axs[0, j].set_title(f'sample {s}', fontsize=FS_TITLE)
for ax in axs[:, 0]:
    ax.set_ylabel('cortical depth (um)', fontsize=FS_LABEL)
axs[0, 0].invert_yaxis()
fig.savefig(OUT_DEPTH)
print(f'wrote {OUT_DEPTH}')

# x-y projection bias and mean-depth difference
xy['log2_rl_lm'] = xy['log2_rl_lm'].where(xy['n_cells'] > 0)
fig, axs = plt.subplots(2, len(SAMPLES), figsize=(4.2 * len(SAMPLES), 7), layout='constrained')
xy_map(fig, axs[0], xy, 'log2_rl_lm', CMAP_LOG2, LOG2_LIM, 'log2((#RL + 1) / (#LM + 1))')
xy_map(fig, axs[1], xy, 'dz', CMAP_DZ, DZ_LIM, 'mean depth RL − LM (um)')
fig.savefig(OUT_XY)
print(f'wrote {OUT_XY}')

# mean-depth difference per bin: data vs shuffled labels
panels = [(f'sample {s}', xy.loc[xy['sample'] == s, 'dz'], null.loc[null['sample'] == s, 'dz'], 'sample')
          for s in SAMPLES] + [('pooled', xy['dz'], null['dz'], 'pooled')]
fig, axs = plt.subplots(1, len(panels), figsize=(4 * len(panels), 3.4), layout='constrained')
for ax, (title, data, shuff, kind) in zip(axs, panels):
    ax.hist(shuff, bins=NULL_BINS[kind], density=True, histtype='step', color='gray', label='shuffled labels')
    ax.hist(data.dropna(), bins=NULL_BINS[kind], density=True, histtype='step', color='C1', lw=1.5, label='data')
    ax.set_title(title, fontsize=FS_TITLE)
    ax.set_xlabel('mean depth diff RL − LM (um)', fontsize=FS_LABEL)
    style(ax)
axs[0].set_ylabel('density', fontsize=FS_LABEL)
axs[-1].legend(fontsize=FS_TICK, frameon=False)
fig.savefig(OUT_NULL)
print(f'wrote {OUT_NULL}')
