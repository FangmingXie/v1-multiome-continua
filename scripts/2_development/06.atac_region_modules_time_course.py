"""A/B/C archetype chromatin region modules (M1-M5) over postnatal day in Yoo25 L2/3 ATAC, NR vs DR, per
cell group of 2_development/03 (all L2/3, A, B, C cells). The ATAC counterpart of 2_development/05.

Modules: L23{A,B,C}_regions_M{1-5}.bed; every region is one consensus peak of the peak matrices (checked).
A region in modules of two archetypes is in each.
Pseudobulk: peak counts summed per (cell group, Sample); cells, Sample, condition, day and A/B/C
group from the 2_development/03 cell table (barcodes are the same strings; every ATAC cell must be there).
lib = the unit's counts over all consensus peaks (CPM in peaks). Per region and sample:
log2(counts / lib * 1e6 + 1) minus the region's reference REF (mean over the P21 NR samples, all L2/3
cells), as in 2_development/05. Module score per sample = mean of that value over the module's regions that pass
the filter in the cell group (>= MIN_PB_COUNTS counts in >= MIN_PB_SAMPLES of the group's samples,
all Ages).
Statistics per module and cell group (l23.timecourse.dr_ols): OLS of the module score over the P12-P21 samples (9 NR, 8 DR),
score ~ age (factor) + DR + DR*(day - D0), D0 = 16: b_DR = average DR effect (at day 16), b_int =
change in the DR effect per day; t-test p-values, BH padj over the 15 modules of the cell group.
Sensitivity: the same scores and OLS with lib = summed obs n_fragment (*_nfrag columns; printed
side by side). Figures: 2_development/08.

Reads:
  data/input/atac/pmat_l23concensus_{Age}.h5ad
  data/input/atac_modules/L23{A,B,C}_regions_M{1-5}.bed
  results/2_development/03.yoo25_l23_cell_groups.tsv
Outputs:
  results/2_development/06.yoo25_l23_atac_region_modules_regions.tsv   (region -> module, pass filter per cell group)
  results/2_development/06.yoo25_l23_atac_region_modules_samples.tsv   (module score per sample)
  results/2_development/06.yoo25_l23_atac_region_modules_stats.tsv     (OLS DR and slope per module)
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

from l23.constants import AGES, TEST_DAYS
from l23.timecourse import dr_ols, bh

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ARCHETYPES   = ['A', 'B', 'C']
MODULES      = ['M1', 'M2', 'M3', 'M4', 'M5']
INPUT_PMAT   = {a: os.path.join(REPO_ROOT, 'data', 'input', 'atac', f'pmat_l23concensus_{a}.h5ad') for a in AGES}
INPUT_BED    = {(k, m): os.path.join(REPO_ROOT, 'data', 'input', 'atac_modules', f'L23{k}_regions_{m}.bed')
                for k in ARCHETYPES for m in MODULES}
INPUT_CELLS  = os.path.join(REPO_ROOT, 'results', '2_development', '03.yoo25_l23_cell_groups.tsv')
OUT_REGIONS  = os.path.join(REPO_ROOT, 'results', '2_development', '06.yoo25_l23_atac_region_modules_regions.tsv')
OUT_SAMPLES  = os.path.join(REPO_ROOT, 'results', '2_development', '06.yoo25_l23_atac_region_modules_samples.tsv')
OUT_STATS    = os.path.join(REPO_ROOT, 'results', '2_development', '06.yoo25_l23_atac_region_modules_stats.tsv')

GROUPS         = {'all L2/3': None, 'A cells': 'A', 'B cells': 'B', 'C cells': 'C'}   # name -> 2_development/03 archetype
REF            = ('all L2/3', 'NR', 21)   # (cells, condition, day) of the per-region reference
MIN_PB_COUNTS  = 10
MIN_PB_SAMPLES = 2

os.makedirs(os.path.dirname(OUT_STATS), exist_ok=True)

# modules
mod = pd.concat([pd.read_csv(p, sep='\t', header=None, usecols=[0, 1, 2], names=['chr', 'start', 'end'])
                 .assign(archetype=k, module=m) for (k, m), p in INPUT_BED.items()], ignore_index=True)
mod['region'] = mod['chr'] + ':' + mod['start'].astype(str) + '-' + mod['end'].astype(str)
if mod.duplicated(['archetype', 'region']).any():
    raise ValueError('a region is listed twice within one archetype region set')
regions = sorted(mod['region'].unique())
print('regions per module:', mod.groupby(['archetype', 'module']).size().unstack().to_string(), sep='\n')
print(f'{len(mod)} module regions, {len(regions)} unique')

# pseudobulk per (cell group, Sample): module-region counts, lib (all peaks), n_fragment
cells = pd.read_csv(INPUT_CELLS, sep='\t', index_col='cell')
peaks, counts, meta = None, [], []
for age, path in INPUT_PMAT.items():
    a = ad.read_h5ad(path)
    if peaks is None:
        peaks = a.var_names
        absent = sorted(set(regions) - set(peaks))
        if absent:
            raise ValueError(f'{len(absent)} BED regions not peaks of {path}: {absent[:5]}')
        col = peaks.get_indexer(regions)
    elif not a.var_names.equals(peaks):
        raise ValueError(f'peaks of {path} differ from those of {next(iter(INPUT_PMAT.values()))}')
    absent = a.obs_names.difference(cells.index)
    if len(absent):
        raise ValueError(f'{len(absent)} cells of {path} not in {INPUT_CELLS}: {list(absent[:3])}')
    c = cells.loc[a.obs_names]
    if (c['Age'] != age).any():
        raise ValueError(f'{path}: cells of Ages {sorted(c["Age"].unique())} in {INPUT_CELLS}, expected {age}')
    X = sp.csr_matrix(a.X)
    lib_cell = np.asarray(X.sum(axis=1)).ravel()
    for grp, k in GROUPS.items():
        in_grp = np.ones(len(c), bool) if k is None else (c['archetype'] == k).to_numpy()
        for s in sorted(c['Sample'].unique()):
            w = in_grp & (c['Sample'] == s).to_numpy()
            counts.append(np.asarray(X[w][:, col].sum(axis=0)).ravel())
            meta.append({'cells': grp, 'Sample': s, 'Age': age, 'condition': c['condition'].iloc[0],
                         'day': int(c['day'].iloc[0]), 'n_cells': int(w.sum()), 'lib': lib_cell[w].sum(),
                         'n_fragment': a.obs['n_fragment'].to_numpy()[w].sum()})
    print(f'{age}: {a.n_obs} cells, samples {sorted(c["Sample"].unique())}')
meta = pd.DataFrame(meta).set_index(['cells', 'Sample'])
pb = pd.DataFrame(np.vstack(counts), index=meta.index, columns=regions)
print(meta.to_string())

# per-region log2(CPM+1) minus the P21 NR reference, for each library
is_ref = ((meta.index.get_level_values('cells') == REF[0]) & (meta['condition'] == REF[1])
          & (meta['day'] == REF[2])).to_numpy()
if not is_ref.any():
    raise ValueError(f'no reference samples {REF}')
rel = {}
for key, lib_col in [('', 'lib'), ('_nfrag', 'n_fragment')]:
    logcpm = np.log2(pb.div(meta[lib_col], axis=0) * 1e6 + 1)
    rel[key] = logcpm - logcpm[is_ref].mean(axis=0)

# region filter per cell group, module scores per sample
passes = pd.DataFrame({g: (pb.loc[g] >= MIN_PB_COUNTS).sum(axis=0) >= MIN_PB_SAMPLES for g in GROUPS})
mod[['archetype', 'module', 'region']].join(passes, on='region').to_csv(OUT_REGIONS, sep='\t', index=False)
print(f'wrote {OUT_REGIONS}')
rows = []
for grp in GROUPS:
    for (k, m), t in mod.groupby(['archetype', 'module']):
        keep = [r for r in t['region'] if passes.loc[r, grp]]
        if not keep:
            raise ValueError(f'{grp}: no region of {k} {m} passes the filter')
        rows.append(meta.loc[grp].assign(score=rel[''].loc[grp, keep].mean(axis=1),
                                         score_nfrag=rel['_nfrag'].loc[grp, keep].mean(axis=1),
                                         n_regions=len(keep), n_regions_set=len(t),
                                         cells=grp, archetype=k, module=m).reset_index())
samp = pd.concat(rows, ignore_index=True)[['cells', 'archetype', 'module', 'Sample', 'Age', 'condition', 'day',
                                           'n_cells', 'lib', 'n_fragment', 'n_regions', 'n_regions_set',
                                           'score', 'score_nfrag']]
samp.to_csv(OUT_SAMPLES, sep='\t', index=False)
print(f'wrote {OUT_SAMPLES}')

# OLS per module and cell group, P12-P21, for each library
stats = []
for (grp, k, m), t in samp[samp['day'].isin(TEST_DAYS)].groupby(['cells', 'archetype', 'module'], sort=False):
    stats.append({'cells': grp, 'archetype': k, 'module': m, 'n_regions': int(t['n_regions'].iloc[0]),
                  'n_NR': int((t['condition'] == 'NR').sum()), 'n_DR': int((t['condition'] == 'DR').sum()),
                  **dr_ols(t, 'score'), **dr_ols(t, 'score_nfrag', suffix='_nfrag')})
stats = pd.DataFrame(stats)
for key in ['DR_d16', 'slope', 'DR_d16_nfrag', 'slope_nfrag']:
    stats[f'padj_{key}'] = stats.groupby('cells')[f'pvalue_{key}'].transform(bh)
stats.to_csv(OUT_STATS, sep='\t', index=False)
print(stats[['cells', 'archetype', 'module', 'n_regions', 'DR_d16', 'padj_DR_d16', 'DR_d16_nfrag', 'padj_DR_d16_nfrag',
             'slope_per_day', 'padj_slope', 'slope_per_day_nfrag', 'padj_slope_nfrag']]
      .to_string(index=False, float_format=lambda v: f'{v:.3g}'))
print('sensitivity (lib = n_fragment vs counts in peaks), Pearson r over the 60 modules: '
      f'DR {np.corrcoef(stats["DR_d16"], stats["DR_d16_nfrag"])[0, 1]:.3f}, '
      f'slope {np.corrcoef(stats["slope_per_day"], stats["slope_per_day_nfrag"])[0, 1]:.3f}')
print(f'wrote {OUT_STATS}')
