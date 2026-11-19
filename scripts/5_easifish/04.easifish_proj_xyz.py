"""EASI-FISH LM vs RL retrograde labels in 3 samples: depth profiles, x-y projection bias, and the per x-y bin
mean-depth difference with a label-shuffle null. Ported from scomics revision/e00/e00.easifish_proj_3samples.

Cells: every segmented cell of each sample; x, y, z in um (z = depth from pia); LM, RL = 0/1 labels.
Groups: LM (LM = 1), RL (RL = 1), both (LM = 1 and RL = 1).
Depth profiles: per sample and group, cells per Z_BINS bin (20 um, 0-440 um, the bins of e00 and of
5_easifish/02; cells deeper than 440 um, in sample 1 only, fall outside and are counted) and the fraction
of the group's binned cells per bin (depth shape).
x-y bins: per sample, XY_STEP squares from the sample's min x and min y (bin = floor((x - x_min) / XY_STEP));
e00 cut with pd.cut, which left out the one cell at the min x and the one at the min y.
Per x-y bin (the full grid of each sample): cell, LM and RL counts, log2((n_RL + 1) / (n_LM + 1)), the mean depth of LM and of RL cells and
dz = mean depth RL - mean depth LM (NaN where either group has no cell).
Null: per sample, N_PERM times, the RL labels and the LM labels are each permuted across the sample's cells
(independently; one seeded stream over the samples) and dz recomputed in the same x-y bins;
p_perm = (1 + #(|dz_null| >= |dz|)) / (1 + N_PERM) over the reps where dz_null is defined.

Reads:
  data/raw/easifish/easifish_sample{1,2,3}_res.csv
Outputs:
  results/5_easifish/04.easifish_depth_hist.tsv   (sample, group, z_lo, z_hi, n, frac)
  results/5_easifish/04.easifish_xy_bins.tsv      (sample, xbin, ybin, x_lo, y_lo, n_cells, n_LM, n_RL, log2_rl_lm,
                                                   mean_z_LM, mean_z_RL, dz, dz_null_mean, p_perm)
  results/5_easifish/04.easifish_xy_dz_null.tsv   (sample, rep, xbin, ybin, dz; bins with dz defined)
"""

import os
import warnings
import numpy as np
import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SAMPLES        = (1, 2, 3)
INPUT_EASIFISH = {s: os.path.join(REPO_ROOT, 'data', 'raw', 'easifish', f'easifish_sample{s}_res.csv')
                  for s in SAMPLES}
RES_DIR        = os.path.join(REPO_ROOT, 'results', '5_easifish')
OUT_DEPTH      = os.path.join(RES_DIR, '04.easifish_depth_hist.tsv')
OUT_XY         = os.path.join(RES_DIR, '04.easifish_xy_bins.tsv')
OUT_NULL       = os.path.join(RES_DIR, '04.easifish_xy_dz_null.tsv')

PROJ_GROUPS = ['LM', 'RL']
Z_BINS      = np.arange(0, 440 + 20, 20)   # um
XY_STEP     = 100                          # um
N_PERM      = 1000
SEED        = 0


def load(path):
    df = pd.read_csv(path)
    missing = {'x', 'y', 'z', *PROJ_GROUPS} - set(df.columns)
    if missing:
        raise ValueError(f'{path} lacks columns {sorted(missing)}')
    bad = {c: sorted(df[c].unique()) for c in PROJ_GROUPS if not set(df[c].unique()) <= {0, 1}}
    if bad:
        raise ValueError(f'expected 0/1 projection labels in {path}, got {bad}')
    df['xbin'] = ((df['x'] - df['x'].min()) // XY_STEP).astype(int)
    df['ybin'] = ((df['y'] - df['y'].min()) // XY_STEP).astype(int)
    return df


def mean_dz(z, xbin, ybin, rl, lm, shape):
    """Mean depth of RL minus LM cells per x-y bin (NaN where either group is empty), as an array of `shape`."""
    flat = xbin * shape[1] + ybin
    n = np.prod(shape)
    means = []
    for lab in (rl, lm):
        cnt = np.bincount(flat[lab], minlength=n)
        tot = np.bincount(flat[lab], weights=z[lab], minlength=n)
        with np.errstate(invalid='ignore', divide='ignore'):
            means.append(np.where(cnt > 0, tot / cnt, np.nan))
    return (means[0] - means[1]).reshape(shape)


def per_bin(values, t, fill=np.nan):
    """A (xbin, ybin)-indexed series aligned to the rows of t."""
    return values.reindex(pd.MultiIndex.from_frame(t[['xbin', 'ybin']]), fill_value=fill).to_numpy()


os.makedirs(RES_DIR, exist_ok=True)
rng = np.random.default_rng(SEED)
depth, xy, null = [], [], []
for s in SAMPLES:
    df = load(INPUT_EASIFISH[s])
    groups = {'LM': df['LM'] == 1, 'RL': df['RL'] == 1, 'both': (df['LM'] == 1) & (df['RL'] == 1)}
    print(f'sample {s}: {len(df)} cells; ' + ', '.join(f'{g} {m.sum()}' for g, m in groups.items())
          + f'; {(df["z"] > Z_BINS[-1]).sum()} cells deeper than {Z_BINS[-1]} um')

    # depth profiles
    for g, m in groups.items():
        n = np.histogram(df.loc[m, 'z'], bins=Z_BINS)[0]
        depth.append(pd.DataFrame({'sample': s, 'group': g, 'z_lo': Z_BINS[:-1], 'z_hi': Z_BINS[1:],
                                   'n': n, 'frac': n / n.sum()}))

    # x-y bins
    shape = (df['xbin'].max() + 1, df['ybin'].max() + 1)
    z, xb, yb = df['z'].to_numpy(), df['xbin'].to_numpy(), df['ybin'].to_numpy()
    rl, lm = groups['RL'].to_numpy(), groups['LM'].to_numpy()
    dz = mean_dz(z, xb, yb, rl, lm, shape)
    dz_null = np.stack([mean_dz(z, xb, yb, rng.permutation(rl), rng.permutation(lm), shape)
                        for _ in range(N_PERM)])
    defined = ~np.isnan(dz_null)
    exceed = (np.abs(dz_null) >= np.abs(dz)) & defined
    with np.errstate(invalid='ignore'):
        p_perm = np.where(np.isnan(dz), np.nan, (1 + exceed.sum(0)) / (1 + defined.sum(0)))
    t = pd.DataFrame({'sample': s,
                      'xbin': np.repeat(np.arange(shape[0]), shape[1]),
                      'ybin': np.tile(np.arange(shape[1]), shape[0])})
    t['x_lo'] = df['x'].min() + t['xbin'] * XY_STEP
    t['y_lo'] = df['y'].min() + t['ybin'] * XY_STEP
    t['n_cells'] = per_bin(df.groupby(['xbin', 'ybin']).size(), t, fill=0)
    for g in PROJ_GROUPS:
        t[f'n_{g}'] = per_bin(df[groups[g]].groupby(['xbin', 'ybin']).size(), t, fill=0)
    for g in PROJ_GROUPS:
        t[f'mean_z_{g}'] = per_bin(df[groups[g]].groupby(['xbin', 'ybin'])['z'].mean(), t)
    t['log2_rl_lm'] = np.log2((t['n_RL'] + 1) / (t['n_LM'] + 1))
    t['dz'] = dz.ravel()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)   # all-NaN bins
        t['dz_null_mean'] = np.nanmean(dz_null, axis=0).ravel()
    t['p_perm'] = p_perm.ravel()
    xy.append(t)
    print(f'  {shape[0]} x {shape[1]} x-y bins, {int(np.isfinite(dz).sum())} with dz; '
          f'median |dz| data {np.nanmedian(np.abs(dz)):.1f} um, null {np.nanmedian(np.abs(dz_null)):.1f} um')

    r, i, j = np.nonzero(defined)
    null.append(pd.DataFrame({'sample': s, 'rep': r, 'xbin': i, 'ybin': j, 'dz': dz_null[r, i, j]}))

depth = pd.concat(depth, ignore_index=True)
depth.to_csv(OUT_DEPTH, sep='\t', index=False)
xy = pd.concat(xy, ignore_index=True)[['sample', 'xbin', 'ybin', 'x_lo', 'y_lo', 'n_cells', 'n_LM', 'n_RL', 'log2_rl_lm',
                                       'mean_z_LM', 'mean_z_RL', 'dz', 'dz_null_mean', 'p_perm']]
xy.to_csv(OUT_XY, sep='\t', index=False)
pd.concat(null, ignore_index=True).to_csv(OUT_NULL, sep='\t', index=False)
for p in (OUT_DEPTH, OUT_XY, OUT_NULL):
    print(f'wrote {p}')
