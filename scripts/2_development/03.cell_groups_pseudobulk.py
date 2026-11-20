"""Yoo25 L2/3 neurons over postnatal development, normally reared (NR) vs dark reared (DR): archetype
cell groups and the pseudobulk counts (all genes, all Ages) used by 2_development/04-06 and 4_dr/02.

Cells: the Yoo25 L2/3 neurons of every Age (0_prep/01). NR is sampled at P6-P21; DR only at P12-P21
(P6-P10 NR are the shared pre-DR baseline).
Expression: raw counts, per cell log2(CP10k+1) (library = sum over all genes).
Archetype cells (obs Type has A/B/C only at P17-P21, so it is not used): every cell of every Age
is scored with the same A, B and C gene sets: the union over all conditions of that archetype's
Table S4 genes (2_development/01) (a gene listed under different archetypes in different conditions is in each of
those sets, overlaps printed), as in 5_easifish/01 but within each Age (each NR / DR condition separately):
per gene log2(CP10k+1) / its SCORE_MAX_PCTILE percentile over the Age's cells, capped at 1, mean over
the set. A / B / C cells of an Age = its cells in the top TOP_FRAC of that score (>= the Age's
1 - TOP_FRAC quantile); the three groups are picked independently, and cells in the top of more
than one score are dropped from all of them (counts printed), so the groups are disjoint.
Pseudobulk: raw counts of every gene summed per sample (the replicate) and cell group (all L2/3,
A, B, C cells), all Ages; lib = the unit's summed library sizes (= row sum of the counts).

Reads:
  data/input/01.yoo25_l23.h5ad
  results/2_development/01.yoo25_l23_abc_genes.tsv
Outputs:
  results/2_development/03.yoo25_l23_cell_groups.tsv            (per cell: archetype group, '' = none)
  results/2_development/03.yoo25_l23_pseudobulk_counts.tsv      ((cells, Sample) rows x all genes)
  results/2_development/03.yoo25_l23_pseudobulk_samples.tsv     ((cells, Sample): Age, condition, day, n_cells, lib)
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

from l23.constants import AGES
from l23.expr import log2_cp10k, pctile_score

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_YOO25   = os.path.join(REPO_ROOT, 'data', 'input', '01.yoo25_l23.h5ad')
INPUT_MARKERS = os.path.join(REPO_ROOT, 'results', '2_development', '01.yoo25_l23_abc_genes.tsv')
OUT_CELLS      = os.path.join(REPO_ROOT, 'results', '2_development', '03.yoo25_l23_cell_groups.tsv')
OUT_PB         = os.path.join(REPO_ROOT, 'results', '2_development', '03.yoo25_l23_pseudobulk_counts.tsv')
OUT_PB_SAMPLES = os.path.join(REPO_ROOT, 'results', '2_development', '03.yoo25_l23_pseudobulk_samples.tsv')

ARCHETYPES       = ['A', 'B', 'C']
TOP_FRAC         = 0.20   # per Age, the top fraction of cells by each score
SCORE_MAX_PCTILE = 98     # per-gene scale: 0 = no expression, this percentile = 1

os.makedirs(os.path.dirname(OUT_PB), exist_ok=True)

src = ad.read_h5ad(INPUT_YOO25)
counts = sp.csr_matrix(src.X)
lib = src.obs['lib'].to_numpy()
cells = src.obs[['Age', 'Sample', 'condition', 'day']].astype({'Age': str, 'Sample': str, 'condition': str})

mk = pd.read_csv(INPUT_MARKERS, sep='\t', usecols=['cond', 'archetype', 'gene'])
missing = sorted(set(AGES) - set(mk['cond']))
if missing:
    raise ValueError(f'Conditions {missing} absent from {INPUT_MARKERS}')
sets = {k: sorted(set(t['gene'])) for k, t in mk.groupby('archetype')}   # union over conditions
if sorted(sets) != ARCHETYPES:
    raise ValueError(f'expected archetypes {ARCHETYPES} in {INPUT_MARKERS}, got {sorted(sets)}')
print(f'archetype genes per score (union of all conditions): { {k: len(v) for k, v in sets.items()} }; '
      f'shared: { {f"{a}{b}": len(set(sets[a]) & set(sets[b])) for i, a in enumerate(ARCHETYPES) for b in ARCHETYPES[i + 1:]} }')
genes = sorted({x for gs in sets.values() for x in gs})
absent = [x for x in genes if x not in src.var_names]
if absent:
    raise ValueError(f'{len(absent)} genes not in {INPUT_YOO25}: {absent}')
col = {x: i for i, x in enumerate(genes)}
L = log2_cp10k(counts[:, src.var_names.get_indexer(genes)].toarray(), lib)   # cells x score genes


top = pd.DataFrame(False, index=cells.index, columns=ARCHETYPES)
for age in AGES:
    m = (cells['Age'] == age).to_numpy()
    for k in ARCHETYPES:
        sc = pctile_score(L[m][:, [col[x] for x in sets[k]]], SCORE_MAX_PCTILE)
        top.loc[m, k] = sc >= np.quantile(sc, 1 - TOP_FRAC)
n_groups = top.sum(axis=1)
print(f'top {TOP_FRAC:.0%} cells per Age and score:')
print(pd.concat([top.groupby(cells['Age']).sum(),
                 (n_groups > 1).groupby(cells['Age']).sum().rename('in >1 group, dropped')], axis=1).to_string())

cells['archetype'] = ''
for k in ARCHETYPES:
    cells.loc[top[k] & (n_groups == 1), 'archetype'] = k
cells[['Age', 'Sample', 'condition', 'day', 'archetype']].rename_axis('cell').to_csv(OUT_CELLS, sep='\t')
print(f'wrote {OUT_CELLS}')

# pseudobulk units = (group, sample); indicator matrix units x cells (cells in the extract's row order)
GROUPS = {'all L2/3': None, **{f'{k} cells': k for k in ARCHETYPES}}   # -> archetype
units = [(grp, smp) for grp, k in GROUPS.items()
         for smp in sorted(cells.loc[cells['archetype'] == k if k else slice(None), 'Sample'].unique())]
unit_idx = {u: i for i, u in enumerate(units)}
r, c = [], []
for j, (smp, k) in enumerate(zip(cells['Sample'], cells['archetype'])):
    for grp, gk in GROUPS.items():
        if gk is None or gk == k:
            r.append(unit_idx[(grp, smp)])
            c.append(j)
M = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(len(units), len(cells)))
pb = (M @ counts).toarray()
if not np.allclose(pb, np.round(pb)):
    raise ValueError(f'{INPUT_YOO25} .X is not integer counts')
index = pd.MultiIndex.from_tuples(units, names=['cells', 'Sample'])
pb = pd.DataFrame(pb.astype(int), index=index, columns=src.var_names)
smp = cells.groupby('Sample')[['Age', 'condition', 'day']].first()
meta = pd.DataFrame({'Age': smp.loc[index.get_level_values('Sample'), 'Age'].to_numpy(),
                     'condition': smp.loc[index.get_level_values('Sample'), 'condition'].to_numpy(),
                     'day': smp.loc[index.get_level_values('Sample'), 'day'].to_numpy(),
                     'n_cells': np.asarray(M.sum(axis=1)).ravel().astype(int),
                     'lib': pb.sum(axis=1).to_numpy()}, index=index)
if not np.allclose(meta['lib'], M @ lib):
    raise ValueError('pseudobulk library sizes differ from the summed per-cell library sizes')
pb.to_csv(OUT_PB, sep='\t')
meta.to_csv(OUT_PB_SAMPLES, sep='\t')
print(meta.groupby(['cells', 'condition', 'day'])['n_cells'].agg(['size', 'sum']).to_string())
print(f'wrote {OUT_PB} {pb.shape}')
print(f'wrote {OUT_PB_SAMPLES}')
