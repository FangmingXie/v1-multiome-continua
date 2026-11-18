"""EASI-FISH projection bias along cortical depth (LM vs RL retrograde labels), 3 samples.

Per sample, the depth distribution of LM+ and of RL+ cells, each normalized to sum 1 over the depth
bins (shape only, not abundance); one panel: mean +/- SD (ddof 0) across the 3 samples.
Depth = column z (um from pia). Z_BINS: 20 um bins from 0 to 440 um (the deepest cell of sample 3,
439 um), so sample 1 cells deeper than 440 um fall outside every bin (L2/3 is drawn at L23_RANGE).
Lines are drawn at each bin's upper edge.

Reads:
  data/raw/easifish/easifish_sample{1,2,3}_res.csv
Outputs:
  figures/5_easifish/01.easifish_proj_depth.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_SAMPLES = [os.path.join(REPO_ROOT, 'data', 'raw', 'easifish', f'easifish_sample{i}_res.csv') for i in (1, 2, 3)]
OUT_FIG       = os.path.join(REPO_ROOT, 'figures', '5_easifish', '01.easifish_proj_depth.pdf')

Z_BINS    = np.arange(0, 440 + 20, 20)   # um
L23_RANGE = (100, 400)                   # um; dashed lines
GROUPS    = {'LM': 'magenta', 'RL': 'green'}   # group -> color

os.makedirs(os.path.dirname(OUT_FIG), exist_ok=True)

dfs = [pd.read_csv(p) for p in INPUT_SAMPLES]
for p, df in zip(INPUT_SAMPLES, dfs):
    bad = {c: sorted(df[c].unique()) for c in GROUPS if not set(df[c].unique()) <= {0, 1}}
    if bad:
        raise ValueError(f'expected 0/1 projection labels in {p}, got {bad}')


def hist(df, group):
    return np.histogram(df.loc[df[group] > 0, 'z'], bins=Z_BINS)[0]


fig, ax = plt.subplots(1, 1, figsize=(1 * 5, 1 * 5))
for g, color in GROUPS.items():
    norm = np.array([hist(df, g) / hist(df, g).sum() for df in dfs])
    mean, sd = norm.mean(axis=0), norm.std(axis=0)
    ax.plot(mean, Z_BINS[1:], '-', color=color, linewidth=3, label=g)
    ax.fill_betweenx(Z_BINS[1:], mean - sd, mean + sd, color=color, alpha=0.3, linewidth=0)
    print(f'{g}: cells per sample {[int((df[g] > 0).sum()) for df in dfs]}')
ax.legend(bbox_to_anchor=(0, -0.05), loc='upper left', frameon=False)
for d in L23_RANGE:
    ax.axhline(d, color='gray', linestyle='--')
ax.grid(False)
ax.set_ylabel('cortical depth (um)')
ax.set_xlabel('fraction of cells per group (mean ± SD, 3 samples)')
ax.invert_yaxis()
ax.xaxis.set_label_position('top')
sns.despine(ax=ax, bottom=True, top=False)
fig.tight_layout()
fig.savefig(OUT_FIG)
print(f'wrote {OUT_FIG}')
