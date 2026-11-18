"""Yoo25 L2/3 neurons over postnatal development, normally reared (NR) vs dark reared (DR): archetype
cell groups and the pseudobulk counts (all genes, all Ages) used by 2_development/03-06 and 4_dr/01.

Cells: Study == '2023 Multiome' (Yoo25), Subclass == 'L2/3' (as in 0_prep), all ages in AGES. NR is
sampled at P6-P21; DR only at P12-P21 (P6-P10 NR are the shared pre-DR baseline). The 2022 RNA
(Cheng22) samples in the same object are left out (other study and chemistry).
Expression: raw counts from .raw, per cell log2(CP10k+1) (library = sum over all genes).
Archetype cells (obs Type has A/B/C only at P17-P21, so it is not used): every cell of every Age
is scored with the same A, B and C gene sets: the union over all conditions of that archetype's
Table S4 genes (0_prep/02) (a gene listed under different archetypes in different conditions is in each of
those sets, overlaps printed), as in 5_easifish/02 but within each Age (each NR / DR condition separately):
per gene log2(CP10k+1) / its SCORE_MAX_PCTILE percentile over the Age's cells, capped at 1, mean over
the set. A' / B' / C' cells of an Age = its cells in the top TOP_FRAC of that score (>= the Age's
1 - TOP_FRAC quantile); the three groups are picked independently, and cells in the top of more
than one score are dropped from all of them (counts printed), so the groups are disjoint.
Pseudobulk: raw counts of every gene summed per sample (the replicate) and cell group (all L2/3,
A', B', C' cells), all Ages; lib = the unit's summed library sizes (= row sum of the counts).

Reads:
  data/raw/superdupermegaRNA_hasraw.h5ad
  results/0_prep/02.yoo25_l23_abc_genes.tsv
Outputs:
  results/2_development/02.yoo25_l23_cell_groups.tsv            (per cell: archetype group, '' = none)
  results/2_development/02.yoo25_l23_pseudobulk_counts.tsv      ((cells, Sample) rows x all genes)
  results/2_development/02.yoo25_l23_pseudobulk_samples.tsv     ((cells, Sample): Age, condition, day, n_cells, lib)
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_SOURCE  = os.path.join(REPO_ROOT, 'data', 'raw', 'superdupermegaRNA_hasraw.h5ad')
INPUT_MARKERS = os.path.join(REPO_ROOT, 'results', '0_prep', '02.yoo25_l23_abc_genes.tsv')
OUT_CELLS      = os.path.join(REPO_ROOT, 'results', '2_development', '02.yoo25_l23_cell_groups.tsv')
OUT_PB         = os.path.join(REPO_ROOT, 'results', '2_development', '02.yoo25_l23_pseudobulk_counts.tsv')
OUT_PB_SAMPLES = os.path.join(REPO_ROOT, 'results', '2_development', '02.yoo25_l23_pseudobulk_samples.tsv')

STUDY    = '2023 Multiome'
SUBCLASS = 'L2/3'
# obs Age -> (postnatal day, condition)
AGES = {'P6': (6, 'NR'), 'P8': (8, 'NR'), 'P10': (10, 'NR'), 'P12': (12, 'NR'), 'P14': (14, 'NR'),
        'P17': (17, 'NR'), 'P21': (21, 'NR'),
        'P12DR': (12, 'DR'), 'P14DR': (14, 'DR'), 'P17DR': (17, 'DR'), 'P21DR': (21, 'DR')}
ARCHETYPES       = ["A'", "B'", "C'"]
TOP_FRAC         = 0.20   # per Age, the top fraction of cells by each score
SCORE_MAX_PCTILE = 98     # per-gene scale: 0 = no expression, this percentile = 1
CHUNK            = 5000   # cells per .raw read

os.makedirs(os.path.dirname(OUT_PB), exist_ok=True)

src = ad.read_h5ad(INPUT_SOURCE, backed='r')
if src.raw is None:
    raise ValueError(f'{INPUT_SOURCE} has no .raw; raw counts are expected there')
obs = src.obs[['Age', 'Study', 'Subclass', 'Sample']].astype(str)
sel = ((obs['Study'] == STUDY) & (obs['Subclass'] == SUBCLASS)).to_numpy()
ages = sorted(obs.loc[sel, 'Age'].unique())
if sorted(AGES) != ages:
    raise ValueError(f'expected {STUDY} L2/3 ages {sorted(AGES)}, found {ages}')
idx = np.flatnonzero(sel)

mk = pd.read_csv(INPUT_MARKERS, sep='\t', usecols=['cond', 'archetype', 'gene'])
missing = sorted(set(AGES) - set(mk['cond']))
if missing:
    raise ValueError(f'Conditions {missing} absent from {INPUT_MARKERS}')
sets = {f"{k}'": sorted(set(t['gene'])) for k, t in mk.groupby('archetype')}   # union over conditions
if sorted(sets) != ARCHETYPES:
    raise ValueError(f'expected archetypes {ARCHETYPES} in {INPUT_MARKERS}, got {sorted(sets)}')
print(f'archetype genes per score (union of all conditions): { {k: len(v) for k, v in sets.items()} }; '
      f'shared: { {f"{a}{b}": len(set(sets[a]) & set(sets[b])) for i, a in enumerate(ARCHETYPES) for b in ARCHETYPES[i + 1:]} }')
genes = sorted({x for gs in sets.values() for x in gs})
absent = [x for x in genes if x not in src.raw.var_names]
if absent:
    raise ValueError(f'{len(absent)} genes not in .raw of {INPUT_SOURCE}: {absent}')
gidx = src.raw.var_names.get_indexer(genes)
col = {x: i for i, x in enumerate(genes)}

lib, X = np.empty(len(idx)), np.empty((len(idx), len(genes)), dtype=np.float32)
for s in range(0, len(idx), CHUNK):
    x = sp.csr_matrix(src.raw.X[idx[s:s + CHUNK]])
    lib[s:s + CHUNK] = np.asarray(x.sum(axis=1)).ravel()
    X[s:s + CHUNK] = x[:, gidx].toarray()
if (lib <= 0).any():
    raise ValueError('cells with zero library size')
L = np.log2(X / lib[:, None] * 1e4 + 1)   # log2(CP10k+1), cells x score genes

cells = obs.iloc[idx].copy()
cells['day'] = cells['Age'].map(lambda a: AGES[a][0])
cells['condition'] = cells['Age'].map(lambda a: AGES[a][1])
if (cells.groupby('Sample')['Age'].nunique() != 1).any():
    raise ValueError('a sample spans more than one Age')


def score(x):
    """Per gene x / (SCORE_MAX_PCTILE percentile), capped at 1; mean over genes (cells x genes in)."""
    hi = np.percentile(x, SCORE_MAX_PCTILE, axis=0)
    hi = np.where(hi > 0, hi, 1.0)   # gene zero in >= 98% of cells: raw value, capped at 1
    return np.clip(x / hi, 0.0, 1.0).mean(axis=1)


top = pd.DataFrame(False, index=cells.index, columns=ARCHETYPES)
for age in AGES:
    m = (cells['Age'] == age).to_numpy()
    for k in ARCHETYPES:
        sc = score(L[m][:, [col[x] for x in sets[k]]])
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

# pseudobulk units = (group, sample); indicator matrix units x cells (cells in .raw row order = idx)
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
pb = np.zeros((len(units), len(src.raw.var_names)))
for s in range(0, len(idx), CHUNK):
    pb += np.asarray((M[:, s:s + CHUNK] @ sp.csr_matrix(src.raw.X[idx[s:s + CHUNK]])).todense())
    print(f'  summed {min(s + CHUNK, len(idx))}/{len(idx)} cells', flush=True)
if not np.allclose(pb, np.round(pb)):
    raise ValueError('.raw.X is not integer counts')
index = pd.MultiIndex.from_tuples(units, names=['cells', 'Sample'])
pb = pd.DataFrame(pb.astype(int), index=index, columns=src.raw.var_names)
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
