"""A', B', C' archetype scores along the L2/3 A'-C' axis (Yoo25 P21 NR snRNA) and along L2/3 depth (P21 NR
MERFISH), and the EASI-FISH LM vs RL retrograde-label depth profiles: binned panels on one normalized y.
Score panels use the MARKER_CONDITION (P21) archetype gene sets (= Table S4; 2_development/01). Figure: 03.

Score per cell = mean over all MARKER_CONDITION genes of that archetype of expression
scaled per gene to [0, 1]: 0 = no expression, 1 = the gene's SCORE_MAX_PCTILE percentile across
cells (values above capped at 1); same scoring as 5_easifish/01:
  snRNA:   log2(CP10k+1), all archetype genes (the 0_prep/01 extract holds every one; checked).
  MERFISH: layer MERFISH_LAYER, the archetype genes on the panel (a subset; counts printed).
Caveat, snRNA panel: the 5_easifish/01 axis is score_C' - score_A' built from (almost) the same A'/C' genes, so
the A' and C' curves there are close to circular by construction; B' is independent of the axis.
The MERFISH panels (depth) are independent of the axis. MERFISH percentiles are per dataset.
Cells are cut into N_BINS equal-count bins along the axis / depth (DEPTH_RANGE); mean +/- SEM (shaded band).
Plotted values are z-scored per archetype score across the cells of each panel (snRNA: all axis cells;
MERFISH: the depth-filtered cells of that dataset) before binning, so panels share one unit (SD of the
score across cells within the dataset). MERFISH detects far fewer transcripts per gene, which shrinks
raw score levels and range; z-scoring removes that scale difference. The TSV keeps both the raw
(mean, sem) and z-scored (z_mean, z_sem) bin values.

EASI-FISH panel: per sample, LM+ and RL+ cell counts in Z_BINS (20 um, 0-440 um) of
column z (um from pia), each normalized to sum 1 over all bins (depth shape only, not
abundance); mean +/- SD (ddof 0) across the 3 samples. Only the bins whose centers fall in EASI_SHOW
are drawn (and written to the TSV); the normalization still uses all bins. Reconciled with the other
panels: (1) all panels share the normalized y below; (2) values are at bin centers, to line up with the
MERFISH bin medians; (3) its band is SD across samples, while the score bands are SEM across cells. MERFISH cells are kept within DEPTH_RANGE.

Normalized y (all panels, 0 to 1; 0 = A' end / superficial, 1 = C' end / deep), linear in each
panel's position variable: snRNA and MERFISH map their cells' Y_NORM_PCTILE percentiles to 0 and 1
(snRNA: all axis cells; MERFISH: the depth-filtered L2/3 cells); EASI-FISH maps EASI_SHOW (50-400 um,
the drawn depth window) to 0 and 1. Linear, so curve shapes are unchanged; it
removes the different units (axis score vs um) and layer thickness / depth-offset differences
between datasets.
The um / axis bounds used are printed and in the panels TSV (y0, y1).

Reads:
  results/5_easifish/01.yoo25_p21nr_ac_axis.tsv
  results/0_prep/01.yoo25_l23.h5ad
  results/2_development/01.yoo25_l23_abc_genes.tsv
  results/0_prep/02.merfish_p21nr_l23.h5ad
  data/raw/easifish/easifish_sample{1,2,3}_res.csv
Outputs:
  results/5_easifish/02.l23_archetype_scores_and_easifish_proj_depth_bins.tsv            (snRNA and MERFISH)
  results/5_easifish/02.l23_archetype_scores_and_easifish_proj_depth_easifish_bins.tsv   (EASI-FISH)
  results/5_easifish/02.l23_archetype_scores_and_easifish_proj_depth_panels.tsv          (per panel: cells,
                                                                                         genes, y = 0 / 1 bounds)
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

from l23.constants import ARCH_COLORS
from l23.expr import log2_cp10k, pctile_score

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_AXIS    = os.path.join(REPO_ROOT, 'results', '5_easifish', '01.yoo25_p21nr_ac_axis.tsv')
INPUT_YOO25   = os.path.join(REPO_ROOT, 'results', '0_prep', '01.yoo25_l23.h5ad')
INPUT_MARKERS = os.path.join(REPO_ROOT, 'results', '2_development', '01.yoo25_l23_abc_genes.tsv')
INPUT_MERFISH = {   # panel label -> L2/3 MERFISH cells (0_prep/02)
    'P21 NR': os.path.join(REPO_ROOT, 'results', '0_prep', '02.merfish_p21nr_l23.h5ad'),
}
INPUT_EASIFISH = [os.path.join(REPO_ROOT, 'data', 'raw', 'easifish', f'easifish_sample{i}_res.csv') for i in (1, 2, 3)]
OUT_PREFIX    = os.path.join(REPO_ROOT, 'results', '5_easifish', '02.l23_archetype_scores_and_easifish_proj_depth')
OUT_TSV       = f'{OUT_PREFIX}_bins.tsv'
OUT_EASI_TSV  = f'{OUT_PREFIX}_easifish_bins.tsv'
OUT_PANELS    = f'{OUT_PREFIX}_panels.tsv'

MARKER_CONDITION = 'P21'
SCORE_MAX_PCTILE = 98    # per-gene scale: 0 = no expression, this percentile = 1
N_BINS           = 15
MERFISH_LAYER    = 'ljnorm'
DEPTH_RANGE      = (50, 450)    # um; MERFISH cell filter
Z_BINS           = np.arange(0, 440 + 20, 20)   # um; EASI-FISH depth bins (the deepest cell is at 439 um)
EASI_SHOW        = (50, 400)    # um; EASI-FISH bins drawn (by bin center), mapped linearly to y = 0-1
Y_NORM_PCTILE    = (1, 99)      # per-panel cell percentiles mapped to y = 0 and 1
PROJ_GROUPS      = ['LM', 'RL']

os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)


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


mk = pd.read_csv(INPUT_MARKERS, sep='\t', usecols=['cond', 'archetype', 'gene'])
mk = mk[mk['cond'] == MARKER_CONDITION]
sets = {f"{k}'": sorted(set(g['gene'])) for k, g in mk.groupby('archetype')}
if sorted(sets) != sorted(ARCH_COLORS):
    raise ValueError(f'expected archetypes {sorted(ARCH_COLORS)} in {INPUT_MARKERS}, got {sorted(sets)}')

# --- Yoo25 P21 NR: scores along the 5_easifish/01 A'-C' axis ---
axis = pd.read_csv(INPUT_AXIS, sep='\t', index_col='cell')
a = ad.read_h5ad(INPUT_YOO25, backed='r')
missing = axis.index.difference(a.obs_names)
if len(missing):
    raise ValueError(f'{len(missing)} axis cells not in {INPUT_YOO25}')
for k, gs in sets.items():
    absent = sorted(set(gs) - set(a.var_names))
    if absent:
        raise ValueError(f'{k} genes not in {INPUT_YOO25}: {absent}')
rows = np.sort(a.obs_names.get_indexer(axis.index))
cells = a.obs_names[rows]
counts = a.X[rows]
counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)
lib = np.asarray(counts.sum(axis=1)).ravel()
pos = axis.loc[cells, 'ac_axis'].to_numpy()
yoo = {}
for k, gs in sets.items():
    idx = [a.var_names.get_loc(g) for g in gs]
    yoo[k] = bin_means(pos, pctile_score(log2_cp10k(counts[:, idx].toarray(), lib), SCORE_MAX_PCTILE), N_BINS)
y_bounds = {'yoo25': unit_bounds(pos)}
for t in yoo.values():
    t['y_norm'] = to_unit(t['pos_median'], y_bounds['yoo25'])
print(f'Yoo25 P21 NR: {len(cells)} cells; genes per score: { {k: len(v) for k, v in sets.items()} }')

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

# --- EASI-FISH: LM / RL depth distributions, each normalized to sum 1 per sample ---
easi = []
for p in INPUT_EASIFISH:
    df = pd.read_csv(p)
    bad = {c: sorted(df[c].unique()) for c in PROJ_GROUPS if not set(df[c].unique()) <= {0, 1}}
    if bad:
        raise ValueError(f'expected 0/1 projection labels in {p}, got {bad}')
    easi.append(df)
z_centers = (Z_BINS[:-1] + Z_BINS[1:]) / 2
y_bounds['EASI-FISH'] = EASI_SHOW
proj = {}
for g in PROJ_GROUPS:
    h = np.array([np.histogram(df.loc[df[g] > 0, 'z'], bins=Z_BINS)[0] for df in easi])
    norm = h / h.sum(axis=1, keepdims=True)
    t = pd.DataFrame({'z_center': z_centers, 'y_norm': to_unit(z_centers, y_bounds['EASI-FISH']),
                      'mean': norm.mean(axis=0), 'sd': norm.std(axis=0)})
    proj[g] = t[t['z_center'].between(*EASI_SHOW)].reset_index(drop=True)
    print(f'EASI-FISH {g}: cells per sample {[int((df[g] > 0).sum()) for df in easi]}')
pd.concat(proj, names=['group', 'bin']).reset_index().to_csv(OUT_EASI_TSV, sep='\t', index=False)

panels = {'yoo25_p21nr_ac_axis': yoo,
          **{f"merfish_{l.replace(' ', '').lower()}_depth": d for l, d in mer.items()}}
tbl = pd.concat({(panel, k): t for panel, d in panels.items() for k, t in d.items()},
                names=['panel', 'archetype', 'bin'])
tbl.reset_index().to_csv(OUT_TSV, sep='\t', index=False)
print(tbl.to_string(float_format=lambda v: f'{v:.3g}'))
print(f'wrote {OUT_TSV}')
print(f'wrote {OUT_EASI_TSV}')
print('y = 0 / 1 bounds (axis units or um):', {k: tuple(round(float(x), 3) for x in v) for k, v in y_bounds.items()})

# per panel: cells, genes per score, the y = 0 / 1 bounds (axis units or um) and how they were set
y_rule = {'marker_condition': MARKER_CONDITION, 'y_pctile_lo': Y_NORM_PCTILE[0], 'y_pctile_hi': Y_NORM_PCTILE[1]}
pan = [{'panel': 'yoo25_p21nr_ac_axis', 'label': 'Yoo25 P21 NR snRNA', 'n_cells': len(cells),
        **{f'n_genes_{k}': len(gs) for k, gs in sets.items()}, 'y0': y_bounds['yoo25'][0],
        'y1': y_bounds['yoo25'][1], **y_rule}]
pan += [{'panel': f"merfish_{l.replace(' ', '').lower()}_depth", 'label': f'MERFISH {l}', 'n_cells': mer_n[l],
         **{f'n_genes_{k}': len(gs) for k, gs in mer_sets[l].items()}, 'y0': y_bounds[l][0], 'y1': y_bounds[l][1],
         **y_rule} for l in mer]
pan += [{'panel': f'easifish_{g}', 'label': 'EASI-FISH', 'n_cells': int(sum((df[g] > 0).sum() for df in easi)),
         'n_cells_per_sample': '/'.join(f'{int((df[g] > 0).sum()):,}' for df in easi),
         'y0': EASI_SHOW[0], 'y1': EASI_SHOW[1]} for g in PROJ_GROUPS]
pd.DataFrame(pan).to_csv(OUT_PANELS, sep='\t', index=False)
print(f'wrote {OUT_PANELS}')
