"""A, B, C archetype scores along the L2/3 A-C axis (Yoo25 P21 NR snRNA) and along L2/3 depth (P21 NR
MERFISH), and the EASI-FISH LM vs RL retrograde labels in 3 samples: depth profiles, x-y projection bias and
the per x-y bin mean-depth difference with a label-shuffle null (ported from scomics
revision/e00/e00.easifish_proj_3samples). Figures: 02 (score and depth panels), 03 (EASI-FISH).

Archetype scores use the MARKER_CONDITION (P21) archetype gene sets (= Table S4; 2_development/01).
Score per cell = mean over the genes of a set of expression scaled per gene to [0, 1]: 0 = no
expression, 1 = the gene's SCORE_MAX_PCTILE percentile across cells (values above capped at 1):
  snRNA:   log2(CP10k+1), all archetype genes (the 0_prep/01 extract holds every one; checked).
  MERFISH: layer MERFISH_LAYER, the archetype genes on the panel (a subset; counts printed).
MERFISH percentiles are per dataset.

A-C axis (Yoo25 P21 NR cells) = score_C - score_A (negative = A end, positive = C end), with
EXCLUDE_GENES (Cdh13 is an A gene, Chrm2 a C gene, Trpc6 is B) removed from the A/C sets so that
these marker genes, when plotted against the axis, do not help define it.
Caveat, snRNA panel: the axis is built from (almost) the same A/C genes as the A and C scores, so
those curves are close to circular by construction; B is independent of the axis. The MERFISH panels
(depth) are independent of the axis.

Score panels: cells are cut into N_BINS equal-count bins along the axis / depth (MERFISH cells kept
within DEPTH_RANGE); mean +/- SEM (shaded band). Plotted values are z-scored per archetype score across
the cells of each panel (snRNA: all axis cells; MERFISH: the depth-filtered cells of that dataset) before
binning, so panels share one unit (SD of the score across cells within the dataset). MERFISH detects far
fewer transcripts per gene, which shrinks raw score levels and range; z-scoring removes that scale
difference. The bins TSV keeps both the raw (mean, sem) and z-scored (z_mean, z_sem) bin values.

EASI-FISH cells: every segmented cell of each sample; x, y, z in um (z = depth from pia); LM, RL = 0/1
labels. Groups: LM (LM = 1), RL (RL = 1), both (LM = 1 and RL = 1).
Depth profiles: per sample and group, cells per Z_BINS bin (20 um, 0-440 um, the bins of e00; cells
deeper than 440 um, in sample 1 only, fall outside and are counted) and the fraction of the group's
binned cells per bin (depth shape, not abundance).
Depth panel (next to the score panels): LM and RL fractions, mean +/- SD (ddof 0) across the 3 samples.
Only the bins whose centers fall in EASI_SHOW are drawn (and written); the fractions still use all bins.
Reconciled with the score panels: (1) all panels share the normalized y below; (2) values are at bin
centers, to line up with the MERFISH bin medians; (3) its band is SD across samples, while the score
bands are SEM across cells.
x-y bins: per sample, XY_STEP squares from the sample's min x and min y (bin = floor((x - x_min) / XY_STEP));
e00 cut with pd.cut, which left out the one cell at the min x and the one at the min y.
Per x-y bin (the full grid of each sample): cell, LM and RL counts, log2((n_RL + 1) / (n_LM + 1)), the mean
depth of LM and of RL cells and dz = mean depth RL - mean depth LM (NaN where either group has no cell).
Null: per sample, N_PERM times, the RL labels and the LM labels are each permuted across the sample's cells
(independently; one seeded stream over the samples) and dz recomputed in the same x-y bins;
p_perm = (1 + #(|dz_null| >= |dz|)) / (1 + N_PERM) over the reps where dz_null is defined.

Normalized y (all panels, 0 to 1; 0 = A end / superficial, 1 = C end / deep), linear in each
panel's position variable: snRNA and MERFISH map their cells' Y_NORM_PCTILE percentiles to 0 and 1
(snRNA: all axis cells; MERFISH: the depth-filtered L2/3 cells); EASI-FISH maps EASI_SHOW (50-400 um,
the drawn depth window) to 0 and 1. Linear, so curve shapes are unchanged; it removes the different
units (axis score vs um) and layer thickness / depth-offset differences between datasets.
The um / axis bounds used are printed and in the panels TSV (y0, y1).

Reads:
  data/input/01.yoo25_l23.h5ad
  results/2_development/01.yoo25_l23_abc_genes.tsv
  data/input/02.merfish_p21nr_l23.h5ad
  data/input/easifish/easifish_sample{1,2,3}_res.csv
Outputs:
  results/5_easifish/01.l23_archetype_scores_bins.tsv     (snRNA and MERFISH score bins)
  results/5_easifish/01.l23_archetype_scores_panels.tsv   (per panel: cells, genes, y = 0 / 1 bounds)
  results/5_easifish/01.easifish_proj_depth_bins.tsv      (LM, RL depth panel: mean, sd across samples)
  results/5_easifish/01.easifish_depth_hist.tsv           (sample, group, z_lo, z_hi, n, frac)
  results/5_easifish/01.easifish_xy_bins.tsv              (sample, xbin, ybin, x_lo, y_lo, n_cells, n_LM, n_RL,
                                                           log2_rl_lm, mean_z_LM, mean_z_RL, dz, dz_null_mean, p_perm)
  results/5_easifish/01.easifish_xy_dz_null.tsv           (sample, rep, xbin, ybin, dz; bins with dz defined)
"""

import os
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

from l23.constants import ARCH_COLORS
from l23.expr import log2_cp10k, pctile_score

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SAMPLES        = (1, 2, 3)
INPUT_YOO25    = os.path.join(REPO_ROOT, 'data', 'input', '01.yoo25_l23.h5ad')
INPUT_MARKERS  = os.path.join(REPO_ROOT, 'results', '2_development', '01.yoo25_l23_abc_genes.tsv')
INPUT_MERFISH  = {   # panel label -> L2/3 MERFISH cells (0_prep/02)
    'P21 NR': os.path.join(REPO_ROOT, 'data', 'input', '02.merfish_p21nr_l23.h5ad'),
}
INPUT_EASIFISH = {s: os.path.join(REPO_ROOT, 'data', 'input', 'easifish', f'easifish_sample{s}_res.csv')
                  for s in SAMPLES}
RES_DIR        = os.path.join(REPO_ROOT, 'results', '5_easifish')
OUT_BINS       = os.path.join(RES_DIR, '01.l23_archetype_scores_bins.tsv')
OUT_PANELS     = os.path.join(RES_DIR, '01.l23_archetype_scores_panels.tsv')
OUT_EASI_BINS  = os.path.join(RES_DIR, '01.easifish_proj_depth_bins.tsv')
OUT_DEPTH      = os.path.join(RES_DIR, '01.easifish_depth_hist.tsv')
OUT_XY         = os.path.join(RES_DIR, '01.easifish_xy_bins.tsv')
OUT_NULL       = os.path.join(RES_DIR, '01.easifish_xy_dz_null.tsv')

AGE              = 'P21'   # Yoo25 P21 NR
MARKER_CONDITION = 'P21'
EXCLUDE_GENES    = ['Cdh13', 'Trpc6', 'Chrm2']   # left out of the A-C axis only
SCORE_MAX_PCTILE = 98    # per-gene scale: 0 = no expression, this percentile = 1
N_BINS           = 15
MERFISH_LAYER    = 'ljnorm'
DEPTH_RANGE      = (50, 450)    # um; MERFISH cell filter
Y_NORM_PCTILE    = (1, 99)      # per-panel cell percentiles mapped to y = 0 and 1
PROJ_GROUPS      = ['LM', 'RL']
Z_BINS           = np.arange(0, 440 + 20, 20)   # um; EASI-FISH depth bins (the deepest cell is at 439 um)
EASI_SHOW        = (50, 400)    # um; EASI-FISH depth bins drawn (by bin center), mapped linearly to y = 0-1
XY_STEP          = 100          # um
N_PERM           = 1000
SEED             = 0

os.makedirs(RES_DIR, exist_ok=True)


def unit_bounds(v):
    return tuple(np.percentile(v, Y_NORM_PCTILE))


def to_unit(v, bounds):
    return (np.asarray(v) - bounds[0]) / (bounds[1] - bounds[0])


def dense(m):
    return m.toarray() if sp.issparse(m) else np.asarray(m)


def bin_means(pos, v, n_bins):
    """Equal-count bins along pos: median pos, mean and SEM of v and of z-scored v (over all cells
    passed in), n cells."""
    edges = np.quantile(pos, np.linspace(0, 1, n_bins + 1))
    b = np.clip(np.searchsorted(edges, pos, side='right') - 1, 0, n_bins - 1)
    z = (v - v.mean()) / v.std()
    df = pd.DataFrame({'pos': pos, 'v': v, 'z': z, 'bin': b}).groupby('bin')
    return pd.DataFrame({'pos_median': df['pos'].median(), 'mean': df['v'].mean(), 'sem': df['v'].sem(),
                         'z_mean': df['z'].mean(), 'z_sem': df['z'].sem(), 'n_cells': df['v'].size()})


def load_easifish(path):
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


mk = pd.read_csv(INPUT_MARKERS, sep='\t', usecols=['cond', 'archetype', 'gene'])
mk = mk[mk['cond'] == MARKER_CONDITION]
sets = {k: sorted(set(g['gene'])) for k, g in mk.groupby('archetype')}
if sorted(sets) != sorted(ARCH_COLORS):
    raise ValueError(f'expected archetypes {sorted(ARCH_COLORS)} in {INPUT_MARKERS}, got {sorted(sets)}')

# --- Yoo25 P21 NR: A-C axis, and the archetype scores along it ---
a = ad.read_h5ad(INPUT_YOO25)
a = a[a.obs['Age'] == AGE].copy()
for k, gs in sets.items():
    absent = sorted(set(gs) - set(a.var_names))
    if absent:
        raise ValueError(f'{k} genes not in {INPUT_YOO25}: {absent}')
counts = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
lib = np.asarray(counts.sum(axis=1)).ravel()


def score(genes):
    idx = [a.var_names.get_loc(g) for g in genes]
    return pctile_score(log2_cp10k(counts[:, idx].toarray(), lib), SCORE_MAX_PCTILE)


axis_sets = {k: sorted(set(sets[k]) - set(EXCLUDE_GENES)) for k in ['A', 'C']}
pos = score(axis_sets['C']) - score(axis_sets['A'])
print(f'Yoo25 P21 NR: {a.n_obs} cells; axis genes { {k: len(v) for k, v in axis_sets.items()} } '
      f'(excluded {EXCLUDE_GENES})')
print(pd.Series(pos, name='ac_axis').groupby(a.obs['Type'].astype(str).to_numpy()).describe().round(3).to_string())
yoo = {k: bin_means(pos, score(gs), N_BINS) for k, gs in sets.items()}
y_bounds = {'yoo25': unit_bounds(pos)}
for t in yoo.values():
    t['y_norm'] = to_unit(t['pos_median'], y_bounds['yoo25'])
print(f'genes per score: { {k: len(v) for k, v in sets.items()} }')

# --- MERFISH NR: scores along depth, one panel per dataset ---
mer, mer_n, mer_sets = {}, {}, {}
for label, path in INPUT_MERFISH.items():
    mf = ad.read_h5ad(path)
    o = mf.obs
    depth = o['depth'].to_numpy()
    in_range = (depth >= DEPTH_RANGE[0]) & (depth <= DEPTH_RANGE[1])
    lay = mf.layers[MERFISH_LAYER]
    mer_sets[label] = {k: [g for g in gs if g in mf.var_names] for k, gs in sets.items()}
    mer[label] = {}
    for k, gs in mer_sets[label].items():
        if not gs:
            raise ValueError(f'no {k} genes on the MERFISH panel of {path}')
        x = dense(lay[:, [mf.var_names.get_loc(g) for g in gs]])[in_range]
        mer[label][k] = bin_means(depth[in_range], pctile_score(x, SCORE_MAX_PCTILE), N_BINS)
    mer_n[label] = in_range.sum()
    y_bounds[label] = unit_bounds(depth[in_range])
    for t in mer[label].values():
        t['y_norm'] = to_unit(t['pos_median'], y_bounds[label])
    print(f'MERFISH {label} L2/3: {mf.n_obs} cells; {in_range.sum()} in depth {DEPTH_RANGE} um; samples '
          f'{sorted(o["sample"].astype(str).unique())}; panel genes per score: '
          f'{ {k: f"{len(mer_sets[label][k])}/{len(sets[k])}" for k in sets} }')

panels = {'yoo25_p21nr_ac_axis': yoo,
          **{f"merfish_{l.replace(' ', '').lower()}_depth": d for l, d in mer.items()}}
tbl = pd.concat({(panel, k): t for panel, d in panels.items() for k, t in d.items()},
                names=['panel', 'archetype', 'bin'])
tbl.reset_index().to_csv(OUT_BINS, sep='\t', index=False)
print(tbl.to_string(float_format=lambda v: f'{v:.3g}'))

# --- EASI-FISH: depth profiles and x-y bins per sample ---
rng = np.random.default_rng(SEED)
easi, depth_hist, xy, null = {}, [], [], []
for s in SAMPLES:
    df = easi[s] = load_easifish(INPUT_EASIFISH[s])
    groups = {'LM': df['LM'] == 1, 'RL': df['RL'] == 1, 'both': (df['LM'] == 1) & (df['RL'] == 1)}
    print(f'EASI-FISH sample {s}: {len(df)} cells; ' + ', '.join(f'{g} {m.sum()}' for g, m in groups.items())
          + f'; {(df["z"] > Z_BINS[-1]).sum()} cells deeper than {Z_BINS[-1]} um')

    # depth profiles
    for g, m in groups.items():
        n = np.histogram(df.loc[m, 'z'], bins=Z_BINS)[0]
        depth_hist.append(pd.DataFrame({'sample': s, 'group': g, 'z_lo': Z_BINS[:-1], 'z_hi': Z_BINS[1:],
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

depth_hist = pd.concat(depth_hist, ignore_index=True)
depth_hist.to_csv(OUT_DEPTH, sep='\t', index=False)
xy = pd.concat(xy, ignore_index=True)[['sample', 'xbin', 'ybin', 'x_lo', 'y_lo', 'n_cells', 'n_LM', 'n_RL', 'log2_rl_lm',
                                       'mean_z_LM', 'mean_z_RL', 'dz', 'dz_null_mean', 'p_perm']]
xy.to_csv(OUT_XY, sep='\t', index=False)
pd.concat(null, ignore_index=True).to_csv(OUT_NULL, sep='\t', index=False)

# depth panel: LM / RL fractions per bin, mean and SD across samples, bins in EASI_SHOW
z_centers = (Z_BINS[:-1] + Z_BINS[1:]) / 2
y_bounds['EASI-FISH'] = EASI_SHOW
proj = {}
for g in PROJ_GROUPS:
    frac = np.stack([depth_hist.loc[(depth_hist['sample'] == s) & (depth_hist['group'] == g), 'frac'].to_numpy()
                     for s in SAMPLES])
    t = pd.DataFrame({'z_center': z_centers, 'y_norm': to_unit(z_centers, y_bounds['EASI-FISH']),
                      'mean': frac.mean(axis=0), 'sd': frac.std(axis=0)})
    proj[g] = t[t['z_center'].between(*EASI_SHOW)].reset_index(drop=True)
pd.concat(proj, names=['group', 'bin']).reset_index().to_csv(OUT_EASI_BINS, sep='\t', index=False)
print('y = 0 / 1 bounds (axis units or um):', {k: tuple(round(float(x), 3) for x in v) for k, v in y_bounds.items()})

# per panel: cells, genes per score, the y = 0 / 1 bounds (axis units or um) and how they were set
y_rule = {'marker_condition': MARKER_CONDITION, 'y_pctile_lo': Y_NORM_PCTILE[0], 'y_pctile_hi': Y_NORM_PCTILE[1]}
pan = [{'panel': 'yoo25_p21nr_ac_axis', 'label': 'Yoo25 P21 NR snRNA', 'n_cells': a.n_obs,
        **{f'n_genes_{k}': len(gs) for k, gs in sets.items()}, 'y0': y_bounds['yoo25'][0],
        'y1': y_bounds['yoo25'][1], **y_rule}]
pan += [{'panel': f"merfish_{l.replace(' ', '').lower()}_depth", 'label': f'MERFISH {l}', 'n_cells': mer_n[l],
         **{f'n_genes_{k}': len(gs) for k, gs in mer_sets[l].items()}, 'y0': y_bounds[l][0], 'y1': y_bounds[l][1],
         **y_rule} for l in mer]
pan += [{'panel': f'easifish_{g}', 'label': 'EASI-FISH', 'n_cells': int(sum((df[g] > 0).sum() for df in easi.values())),
         'n_cells_per_sample': '/'.join(f'{int((df[g] > 0).sum()):,}' for df in easi.values()),
         'y0': EASI_SHOW[0], 'y1': EASI_SHOW[1]} for g in PROJ_GROUPS]
pd.DataFrame(pan).to_csv(OUT_PANELS, sep='\t', index=False)
for p in (OUT_BINS, OUT_PANELS, OUT_EASI_BINS, OUT_DEPTH, OUT_XY, OUT_NULL):
    print(f'wrote {p}')
