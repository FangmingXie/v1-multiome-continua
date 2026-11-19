"""Conserved L2/3 A/C genes in situ, mouse vs human MERFISH: L2/3 cells on one relative depth scale, and depth profiles.

Genes: PROFILE_GENES (mouse symbol, a row per side; A = superficial, C = deep), mapped to human by 1-to-1
orthologs (data/raw/human_mouse_orthologs.tsv); every gene must be on both MERFISH panels (checked).
INSITU_GENES (a subset) are drawn in situ by 08.
Sections: every NR mouse section (PROFILE_MOUSE_SECTIONS) and the Fang22 human sections with >= 300
eL2/3.IT cells in which every INSITU_GENES gene is detected in > H_MIN_DETECT of the L2/3 cells
(PROFILE_HUMAN_SECTIONS; checked). Three other 4000-gene sections with >= 300 eL2/3.IT cells barely detect
KCNH5 (< 7% of L2/3 cells) and are not used: H18.06.006.MTG.4000.expand.rep3, H22.26.401.MTG.4000.expand.rep1,
H20.30.001.STG.4000.expand.rep1.
Cells and expression:
  mouse  P28 NR sections: the L2/3 MERFISH cells of the 0_prep/02 extract with depth from pia in
         MOUSE_DEPTH_RANGE (um; drops label-transfer outliers); expression log2(jnorm + 1) (the section's
         500-gene panel); position = width x depth from pia (um).
  human  Fang22 section (250 or 4000 genes): eL2/3.IT cells; counts scaled to the median total per cell, log1p;
         position = pre-expansion adjusted.x/y (um), rotated so pia is up (superficial eL2/3.IT above deep cells).
Every section: depth is re-measured from a 2nd-order polynomial fit to the top of L2/3 (CURVE_PCTILE depth
percentile in CURVE_N_BINS equal-count bins along the section); cells whose depth then lies outside
Tukey's fences (Q1 - TUKEY_K IQR, Q3 + TUKEY_K IQR) are removed. Relative L2/3 depth = (depth - p2) /
(p98 - p2) over the kept cells (0 = top, 1 = bottom of the L2/3 band), and position along the section in
the same unit (L2/3 thicknesses), so mouse (~300 um band) and human (~1500 um band) share one scale.
Profiles: per gene and section, mean log expression in N_REL_BINS equal-width relative-depth bins (cells
outside [0, 1] left out), divided by its highest bin (peak = 1). Sections, not cells, are the replicates
(08 plots mean +- SEM over sections): pooling cells first would weight sections by cell count and mix panels
of different sensitivity (250 vs 4000 genes), whereas per-section peak normalization removes those offsets.

Reads:
  data/raw/human_mouse_orthologs.tsv
  results/0_prep/02.merfish_p28nr_l23.h5ad                                     (L2/3 cells, depth, expression)
  data/raw/fang22_human_merfish/{section}.{features,genes,matrix}.csv          (every PROFILE_HUMAN_SECTIONS)
Outputs:
  results/1_mouse_human/06.l23_merfish_genes.tsv              (side, human_symbol, mouse_symbol, insitu)
  results/1_mouse_human/06.l23_merfish_sections.tsv           (every section: cells, outliers, L2/3 thickness,
                                                               top-of-L2/3 fit, INSITU_GENES detection)
  results/1_mouse_human/06.l23_merfish_top_fit_points.tsv     (every section: the points the top was fit to)
  results/1_mouse_human/06.l23_merfish_cells.parquet          (every section, every L2/3 cell: raw position,
                                                               outlier flag, relative position, log expression)
  results/1_mouse_human/06.l23_merfish_rel_depth_profiles.tsv (every section: gene x bin profiles)
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import coo_matrix

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RES_DIR    = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
IN_ORTHO   = os.path.join(REPO_ROOT, 'data', 'raw', 'human_mouse_orthologs.tsv')
IN_MOUSE   = os.path.join(REPO_ROOT, 'results', '0_prep', '02.merfish_p28nr_l23.h5ad')
HUMAN_DIR  = os.path.join(REPO_ROOT, 'data', 'raw', 'fang22_human_merfish')
OUT_GENES    = os.path.join(RES_DIR, '06.l23_merfish_genes.tsv')
OUT_SECTIONS = os.path.join(RES_DIR, '06.l23_merfish_sections.tsv')
OUT_POINTS   = os.path.join(RES_DIR, '06.l23_merfish_top_fit_points.tsv')
OUT_CELLS    = os.path.join(RES_DIR, '06.l23_merfish_cells.parquet')
OUT_PROFILES = os.path.join(RES_DIR, '06.l23_merfish_rel_depth_profiles.tsv')

PROFILE_GENES = {   # mouse symbol, by side (A = superficial, C = deep); all on both MERFISH panels
    'A': ['Cdh13', 'Robo1', 'Epha6'],
    'C': ['Kcnh5', 'Tox', 'Lrp1b'],
}
INSITU_GENES = ['Cdh13', 'Kcnh5']   # drawn in situ (08); the human-section detection filter uses them
H_MIN_DETECT = 0.10   # every human section must detect every INSITU_GENES gene in > this fraction of its
                      # L2/3 cells (checked)
MOUSE_DEPTH_RANGE = (50, 450)      # um from pia; ~ the 5-95% depth range of L2/3 MERFISH cells
MOUSE_LAYER    = 'jnorm'
H_SUBCLASS, H_SUBCLASS_COL = 'eL2/3.IT', 'cluster_L2'
H_DEEP = ['eL6.IT', 'eL6b', 'eL6.CT', 'eL6.IT.CAR3', 'eL5/6.NP', 'eL6.CAR3', 'eNP']   # orient the section
H_SRC_X, H_SRC_Y = 'adjusted.x', 'adjusted.y'   # pre-expansion um
H_MIN_ANCHOR = 5
PROFILE_MOUSE_SECTIONS = ['P28NRa_ant', 'P28NRa_pos', 'P28NRb_ant', 'P28NRb_pos']   # all NR
PROFILE_HUMAN_SECTIONS = [   # Fang22: >= 300 eL2/3.IT cells and INSITU_GENES detected (see the docstring)
    'H18.06.006.MTG.250.expand.rep1', 'H18.06.006.MTG.4000.expand.rep1',
    'H22.26.401.MTG.250.unexpand.rep1',
    'H20.30.001.STG.250.unexpand.rep1', 'H20.30.001.STG.250.unexpand.rep2',
    'H20.30.001.STG.4000.expand.rep2',
]
CURVE_N_BINS  = 10    # top-of-L2/3 fit (every section): equal-count bins along the section ...
CURVE_PCTILE  = 5     # ... each giving (median x, this depth percentile); a 2nd-order polynomial fit to these
TUKEY_K     = 1.5
REL_PCTILE  = (2, 98)  # relative L2/3 depth: these depth percentiles of the section's kept L2/3 cells -> 0, 1
N_REL_BINS  = 5        # profile: equal-width bins over relative depth [0, 1]

os.makedirs(RES_DIR, exist_ok=True)

# --- genes: sides and 1-to-1 human orthologs, on both panels ---
ortho = pd.read_csv(IN_ORTHO, sep='\t').drop_duplicates('human_symbol').drop_duplicates('mouse_symbol')
SIDE = {g: s for s, gs in PROFILE_GENES.items() for g in gs}
missing = [g for g in INSITU_GENES if g not in SIDE]
if missing:
    raise ValueError(f'INSITU_GENES {missing} not in PROFILE_GENES')
genes = ortho[ortho['mouse_symbol'].isin(SIDE)].copy()
missing = set(SIDE) - set(genes['mouse_symbol'])
if missing:
    raise ValueError(f'no 1-to-1 human ortholog for {sorted(missing)} in {IN_ORTHO}')
genes['side'] = genes['mouse_symbol'].map(SIDE)
genes['insitu'] = genes['mouse_symbol'].isin(INSITU_GENES)
genes = genes.set_index('mouse_symbol').loc[list(SIDE)].reset_index()[['side', 'human_symbol', 'mouse_symbol', 'insitu']]
genes.to_csv(OUT_GENES, sep='\t', index=False)
print(f'wrote {OUT_GENES}')


def tukey_keep(depth):
    q1, q3 = np.percentile(depth, [25, 75])
    return (depth >= q1 - TUKEY_K * (q3 - q1)) & (depth <= q3 + TUKEY_K * (q3 - q1))


def flatten_top(pos, section):
    """Depth measured from a 2nd-order polynomial fit to the top of L2/3: per CURVE_N_BINS equal-count bins along
    x, the CURVE_PCTILE depth percentile; depth -> depth - top(x) (vertical offset; the curvature is gentle).
    Returns (flattened positions, fit points, polynomial)."""
    b = pd.qcut(pos['x'], CURVE_N_BINS, labels=False)
    pts = pd.DataFrame({'x': pos['x'].groupby(b).median(),
                        'depth': pos['depth'].groupby(b).quantile(CURVE_PCTILE / 100)})
    top = np.poly1d(np.polyfit(pts['x'], pts['depth'], 2))
    print(f'{section}: top of L2/3 fit {top.coeffs[0]:.2e} x^2 + {top.coeffs[1]:.2e} x + {top.coeffs[2]:.1f}; '
          f'range of the fit over the section {np.ptp(top(pos["x"])):.0f} um')
    return pos.assign(depth=pos['depth'] - top(pos['x'])), pts, top


def finish(species, section, pos, ex):
    """Flatten, drop depth outliers, relative scale. Returns (cells table, fit points, section summary)."""
    flat, pts, top = flatten_top(pos, section)
    keep = tukey_keep(flat['depth'].to_numpy())
    lo, hi = np.percentile(flat.loc[keep, 'depth'], REL_PCTILE)
    cells = pd.DataFrame({'species': species, 'section': section, 'cell': pos.index, 'x': pos['x'].to_numpy(),
                          'depth': pos['depth'].to_numpy(), 'outlier': ~keep,
                          'rel_depth': np.where(keep, (flat['depth'] - lo) / (hi - lo), np.nan),
                          'rel_x': np.where(keep, (flat['x'] - flat.loc[keep, 'x'].min()) / (hi - lo), np.nan)})
    cells = pd.concat([cells, ex.set_axis(genes['mouse_symbol'], axis=1).reset_index(drop=True)], axis=1)
    info = {'species': species, 'section': section, 'n_cells': int(keep.sum()), 'n_outliers': int((~keep).sum()),
            'l23_thickness_um': hi - lo, 'top_c2': top.coeffs[0], 'top_c1': top.coeffs[1], 'top_c0': top.coeffs[2]}
    print(f'{species} {section}: {len(pos)} L2/3 cells; {(~keep).sum()} depth outliers removed')
    return cells, pts.assign(species=species, section=section), info


# --- mouse: L2/3 NR cells (0_prep/02, in section-file order) in the depth range ---
mad = ad.read_h5ad(IN_MOUSE)
absent = [g for g in genes['mouse_symbol'] if g not in mad.var_names]
if absent:
    raise ValueError(f'{absent} not on the panel of {IN_MOUSE}')
depth = mad.obs['depth'].to_numpy()
mad = mad[(depth >= MOUSE_DEPTH_RANGE[0]) & (depth <= MOUSE_DEPTH_RANGE[1])].copy()
mad.obs_names.name = 'cell'
print(f'mouse P28 NR L2/3 MERFISH cells in {MOUSE_DEPTH_RANGE} um: {mad.n_obs}')


def load_mouse(section):
    """L2/3 cells of one NR section: (positions, log2(jnorm+1) of the genes)."""
    sd = mad[mad.obs['sample'] == section]
    if sd.n_obs == 0:
        raise ValueError(f'no {section} cells in {IN_MOUSE}')
    lay = np.asarray(sd.layers[MOUSE_LAYER])
    ex = pd.DataFrame(np.log2(lay[:, [sd.var_names.get_loc(g) for g in genes['mouse_symbol']]] + 1),
                      index=sd.obs_names, columns=genes['mouse_symbol'])
    pos = pd.DataFrame({'x': sd.obs['width'].to_numpy(), 'depth': sd.obs['depth'].to_numpy()}, index=sd.obs_names)
    return pos, ex


def load_human(section):
    """eL2/3.IT cells of one Fang22 section: (positions, log1p median-normalized expression of the genes);
    rotated so pia is up (superficial eL2/3.IT above the deep H_DEEP cells)."""
    feat = pd.read_csv(os.path.join(HUMAN_DIR, f'{section}.features.csv'))
    hg = pd.read_csv(os.path.join(HUMAN_DIR, f'{section}.genes.csv'))['name'].tolist()
    absent = [g for g in genes['human_symbol'] if g not in hg]
    if absent:
        raise ValueError(f'{section}: {absent} not on the panel')
    trip = pd.read_csv(os.path.join(HUMAN_DIR, f'{section}.matrix.csv'))
    counts = coo_matrix((trip['val'].values, (trip['col'].values - 1, trip['row'].values - 1)),
                        shape=(len(feat), len(hg))).tocsr()
    total = np.asarray(counts.sum(axis=1)).ravel()
    if (total == 0).any():
        raise ValueError(f'{section}: {(total == 0).sum()} cells with zero total counts')
    norm = counts.multiply((np.median(total) / total)[:, None]).tocsr()
    norm.data = np.log1p(norm.data)
    sup, deep = feat[feat[H_SUBCLASS_COL] == H_SUBCLASS], feat[feat[H_SUBCLASS_COL].isin(H_DEEP)]
    if len(sup) < H_MIN_ANCHOR or len(deep) < H_MIN_ANCHOR:
        raise ValueError(f'{section}: too few anchor cells to orient ({len(sup)} {H_SUBCLASS}, {len(deep)} deep)')
    deg = (np.degrees(np.arctan2(sup[H_SRC_Y].mean() - deep[H_SRC_Y].mean(),
                                 sup[H_SRC_X].mean() - deep[H_SRC_X].mean())) - 90) % 360
    rad = np.deg2rad(deg)
    hx = feat[H_SRC_X] * np.cos(rad) + feat[H_SRC_Y] * np.sin(rad)
    hy = -feat[H_SRC_X] * np.sin(rad) + feat[H_SRC_Y] * np.cos(rad)
    m = (feat[H_SUBCLASS_COL] == H_SUBCLASS).to_numpy()
    ex = pd.DataFrame(norm[m][:, [hg.index(g) for g in genes['human_symbol']]].toarray(),
                      index=feat.loc[m, 'name'], columns=genes['human_symbol'])
    pos = pd.DataFrame({'x': hx[m].to_numpy(), 'depth': -hy[m].to_numpy()}, index=ex.index)
    print(f'human {section}: rotated {deg:.1f} deg')
    return pos, ex


done = [finish('mouse', sec, *load_mouse(sec)) for sec in PROFILE_MOUSE_SECTIONS]
done += [finish('human', sec, *load_human(sec)) for sec in PROFILE_HUMAN_SECTIONS]

# human sections: detection of INSITU_GENES in the kept L2/3 cells
sections = pd.DataFrame([d[2] for d in done])
for g in INSITU_GENES:
    sections[f'detect_{g}'] = [(c.loc[~c['outlier'], g] > 0).mean() for c, _, _ in done]
h = sections['species'] == 'human'
sections['detect_threshold'] = H_MIN_DETECT
print(f'human sections, fraction of L2/3 cells detecting {INSITU_GENES} (each must be > {H_MIN_DETECT}):')
print(sections[h].set_index('section')[[f'detect_{g}' for g in INSITU_GENES]].round(3).to_string())
low = sections.loc[h & (sections[[f'detect_{g}' for g in INSITU_GENES]] <= H_MIN_DETECT).any(axis=1), 'section']
if len(low):
    raise ValueError(f'human sections {list(low)} detect one of {INSITU_GENES} in <= {H_MIN_DETECT} of L2/3 cells')
sections.to_csv(OUT_SECTIONS, sep='\t', index=False)
(pd.concat([d[1] for d in done])[['species', 'section', 'x', 'depth']].assign(depth_pctile=CURVE_PCTILE)
 .to_csv(OUT_POINTS, sep='\t', index=False))
cells = pd.concat([d[0] for d in done], ignore_index=True)
cells.to_parquet(OUT_CELLS, index=False)
print(f'wrote {OUT_SECTIONS}\nwrote {OUT_POINTS}\nwrote {OUT_CELLS}')

# --- depth profiles: per gene and section, bin means over relative depth / peak bin ---
edges = np.linspace(0, 1, N_REL_BINS + 1)
mids = (edges[:-1] + edges[1:]) / 2
rows = []
for mg in genes['mouse_symbol']:
    for c, _, info in done:
        sp, sec, c = info['species'], info['section'], c[~c['outlier']]
        inb = c['rel_depth'].between(0, 1)
        b = np.clip(np.searchsorted(edges, c.loc[inb, 'rel_depth'], side='right') - 1, 0, N_REL_BINS - 1)
        prof = c.loc[inb, mg].groupby(b).mean().reindex(range(N_REL_BINS))
        rows.append(pd.DataFrame({'gene': mg, 'species': sp, 'section': sec, 'bin_mid': mids,
                                  'rel_expr': (prof / prof.max()).to_numpy()}))
pd.concat(rows).to_csv(OUT_PROFILES, sep='\t', index=False)
print(f'wrote {OUT_PROFILES}')
