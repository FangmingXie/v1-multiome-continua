"""Yoo25 L2/3 over postnatal development, NR vs DR: cells in P8 and P21 PC space of the global HVGs, and the
spread of each Age along the A-C and B axes of the P21 archetype triangle.
Cleaned from scomics revision/imported/analysis_multiome_py/31.multiome_basics_pcs_global_hvgs.

Cells: the Yoo25 L2/3 neurons of every Age (0_prep/01).
Genes: drop ^mt- and Xist, keep genes detected in > MIN_CELLS cells over all Ages (the rule of 2_development/01
step 3). Expression: CP10k with lib = sum over all genes (obs lib); log2(CP10k+1).
HVGs: genes split into N_BIN quantile bins of mean CP10k; per bin, the top HVG_FRAC by CP10k variance / mean.
PCs: PCA of the HVG log2(CP10k+1) fit on the P_EARLY cells (PC1) and on the P_LATER cells (PC1, PC2),
then applied to every cell (signs flipped as in the original; randomized SVD, seeded). Null: PCA fit
on the P_LATER cells with each gene shuffled across cells (seeded), applied to the real data of every cell.
Archetypes: PCHA (NOC, DELTA, seeded) on the P_LATER cells in their (PC1, PC2). A-C axis = unit vector
between the two vertices with the lowest and highest PC1, B axis = its perpendicular; signs set so that
the A-C axis increases with score C - A and the B axis with score B (2_development/01 scores); the vertices are
labeled A, C (ends of the A-C axis) and B (the third). Per cell: ac, b = (PC1, PC2) projected on the axes
(real and null PCs). Variance (ddof 1) of pc_early1, ac, b, ac_shuff, b_shuff per Age and per Sample.

Reads:
  data/input/01.yoo25_l23.h5ad
  results/2_development/01.yoo25_l23_abc_scores.tsv
Outputs:
  results/2_development/02.yoo25_l23_hvg_pcs_cells.tsv        (per cell: PCs, ac, b, null ac, b)
  results/2_development/02.yoo25_l23_hvg_pcs_archetypes.tsv   (A, B, C vertices in P_LATER (PC1, PC2))
  results/2_development/02.yoo25_l23_hvg_pcs_variance.tsv     (variances per Age and per Sample)
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from sklearn.decomposition import PCA

from l23.constants import AGES
from l23.expr import log2_cp10k
from l23.utils import pcha, shuffle_rows_per_col

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_YOO25  = os.path.join(REPO_ROOT, 'data', 'input', '01.yoo25_l23.h5ad')
INPUT_SCORES = os.path.join(REPO_ROOT, 'results', '2_development', '01.yoo25_l23_abc_scores.tsv')
OUT_CELLS    = os.path.join(REPO_ROOT, 'results', '2_development', '02.yoo25_l23_hvg_pcs_cells.tsv')
OUT_ARCH     = os.path.join(REPO_ROOT, 'results', '2_development', '02.yoo25_l23_hvg_pcs_archetypes.tsv')
OUT_VAR      = os.path.join(REPO_ROOT, 'results', '2_development', '02.yoo25_l23_hvg_pcs_variance.tsv')

MIN_CELLS = 10       # keep genes detected in more than this many cells
N_BIN     = 20       # mean-expression bins for HVG selection
HVG_FRAC  = 0.3      # fraction of genes kept per bin
P_EARLY   = 'P8'     # Age the early PCs are fit on
P_LATER   = 'P21'    # Age the later PCs and the archetypes are fit on
NOC       = 3        # archetypes
DELTA     = 0.2      # PCHA relaxation
SEED      = 0
VAR_COLS  = ['pc_early1', 'ac', 'b', 'ac_shuff', 'b_shuff']


def select_hvgs(cp10k):
    """Indices of the HVGs: per quantile bin of mean CP10k, the top HVG_FRAC genes by variance / mean."""
    mean = np.ravel(cp10k.mean(axis=0))
    var = np.ravel(cp10k.multiply(cp10k).mean(axis=0)) - mean ** 2
    g = pd.DataFrame({'bin': pd.qcut(mean, N_BIN, labels=False), 'ratio': var / mean})
    n_per_bin = int(HVG_FRAC * len(g) / N_BIN)
    top = g.groupby('bin')['ratio'].nlargest(n_per_bin)
    return np.sort(top.index.get_level_values(1).to_numpy())


def fit_project(x_fit, x_all, n_pcs, flip):
    """PCA fit on x_fit, applied to x_all; columns in `flip` (0-based) change sign."""
    pcs = PCA(n_components=n_pcs, random_state=SEED).fit(x_fit).transform(x_all)
    pcs[:, flip] *= -1
    return pcs


def unit(v):
    return v / np.linalg.norm(v)


os.makedirs(os.path.dirname(OUT_CELLS), exist_ok=True)

# cells and genes
src = ad.read_h5ad(INPUT_YOO25)
obs = src.obs[['Age', 'Sample', 'condition', 'day']].astype({'Age': str, 'Sample': str})
unknown = sorted(set(obs['Age']) - set(AGES))
if unknown:
    raise ValueError(f'Ages {unknown} not in l23.constants.AGES')
counts = sp.csr_matrix(src.X)
lib = src.obs['lib'].to_numpy()
genes = src.var_names
keep = ~genes.str.contains(r'^mt-|^Xist$') & (np.ravel((counts > 0).sum(axis=0)) > MIN_CELLS)
counts, genes = counts[:, keep], genes[keep]
print(f'{counts.shape[0]} cells, {counts.shape[1]} genes after the mt/Xist and > {MIN_CELLS} cells filter')

scores = pd.read_csv(INPUT_SCORES, sep='\t', index_col='cell')
if not scores.index.sort_values().equals(obs.index.sort_values()):
    raise ValueError(f'{INPUT_SCORES} cells differ from {INPUT_YOO25}')
scores = scores.reindex(obs.index)

# HVGs
cp10k = sp.diags(1 / lib) @ counts * 1e4
hvg = select_hvgs(cp10k)
print(f'{len(hvg)} HVGs')
x = log2_cp10k(counts[:, hvg].toarray(), lib)

# PCs fit on the early / later Age, applied to every cell
early = (obs['Age'] == P_EARLY).to_numpy()
later = (obs['Age'] == P_LATER).to_numpy()
pc_early = fit_project(x[early], x, 1, flip=[0])
pc_later = fit_project(x[later], x, 2, flip=[0, 1])
pc_shuff = fit_project(shuffle_rows_per_col(x[later], seed=SEED), x, 2, flip=[])

# archetypes of the later Age, A-C and B axes
np.random.seed(SEED)                    # py_pcha initializes with np.random
aa, varexpl = pcha(pc_later[later].T, noc=NOC, delta=DELTA)   # vertices sorted by PC1
print(f'PCHA on {P_LATER} cells: variance explained {varexpl:.3f}')
t1 = unit(aa[:, 2] - aa[:, 0])
t2 = np.array([t1[1], -t1[0]])
if np.corrcoef(pc_later[later] @ t1, (scores['score_c'] - scores['score_a'])[later])[0, 1] < 0:
    t1 = -t1
if np.corrcoef(pc_later[later] @ t2, scores['score_b'][later])[0, 1] < 0:
    t2 = -t2
proj_t1 = aa.T @ t1
vertex = {'A': int(np.argmin(proj_t1[[0, 2]]) * 2), 'C': int(np.argmax(proj_t1[[0, 2]]) * 2), 'B': 1}
arch = pd.DataFrame({'archetype': list(vertex),
                     'pc_later1': [aa[0, i] for i in vertex.values()],
                     'pc_later2': [aa[1, i] for i in vertex.values()]})
arch.to_csv(OUT_ARCH, sep='\t', index=False)
print(arch.to_string(index=False))

# per cell
cells = obs[['Age', 'Sample']].copy()
cells['pc_early1'] = pc_early[:, 0]
cells['pc_later1'] = pc_later[:, 0]
cells['pc_later2'] = pc_later[:, 1]
cells['ac'] = pc_later @ t1
cells['b'] = pc_later @ t2
cells['ac_shuff'] = pc_shuff @ t1
cells['b_shuff'] = pc_shuff @ t2
cells.index.name = 'cell'
cells.to_csv(OUT_CELLS, sep='\t')
print(f'wrote {OUT_CELLS} {cells.shape}')
for k, s in [('score C - A', scores['score_c'] - scores['score_a']), ('score B', scores['score_b'])]:
    print(f'{P_LATER}: r(ac, {k}) = {np.corrcoef(cells["ac"][later], s[later])[0, 1]:+.2f}, '
          f'r(b, {k}) = {np.corrcoef(cells["b"][later], s[later])[0, 1]:+.2f}')

# variance per Age and per Sample
var = []
for level in ['Age', 'Sample']:
    v = cells.groupby(level, observed=True)[VAR_COLS].var()
    v.insert(0, 'unit', v.index)
    v.insert(0, 'level', level)
    v['Age'] = cells.groupby(level, observed=True)['Age'].first().to_numpy() if level == 'Sample' else v['unit']
    var.append(v)
var = pd.concat(var, ignore_index=True)
var['condition'] = var['Age'].map(lambda a: AGES[a][1])
var['day'] = var['Age'].map(lambda a: AGES[a][0])
var['Age'] = pd.Categorical(var['Age'], categories=list(AGES))
var = var.sort_values(['level', 'Age', 'unit'])[['level', 'unit', 'Age', 'condition', 'day'] + VAR_COLS]
var.to_csv(OUT_VAR, sep='\t', index=False)
print(var[var['level'] == 'Age'].to_string(index=False))
print(f'wrote {OUT_VAR}')
