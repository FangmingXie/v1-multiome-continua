"""Yoo25 L2/3 archetype (A/B/C) genes per condition (= Table S4), from raw counts.

Cells: Yoo25 ('2023 Multiome') L2/3, all 11 conditions (cond = Age: NR P6-P21, DR P12-P21).
Steps (per condition unless noted):
  1) scores. Cheng22 L2/3 archetype genes (P17on = A/B/C) present in Yoo25. log2(CP10k+1) with
     lib = sum over ALL genes; each scored gene min-max scaled over the condition's cells (+EPS);
     score_X = mean over the X genes.
  2) archetype cells. A = the N_CELLS lowest (score_c - score_a), C = the highest, B = the highest score_b;
     a cell in more than one group is dropped from all.
  3) archetype genes. Test genes: drop ^mt- and Xist, keep genes detected in > MIN_CELLS cells
     (over all conditions); CP10k with lib = sum over the TEST genes (not step 1's lib). For each pair of
     archetypes: Student's t-test on log2(CP10k+1) (NaN p = 1), BH per pair; log2FC = log2(mean CP10k +
     OFFSET) difference. Gene of X: log2FC > 0 vs both others, max log2FC > L2FC_TH, min q < ALPHA;
     for B also min log2FC > L2FC_TH_B_MIN.

Reads:
  data/raw/superdupermegaRNA_hasraw.h5ad        (Yoo25 counts in .raw)
  data/raw/Cheng22_Xie25_L23_ABC_genes_n286.csv     (Cheng22 adult L2/3 A/B/C genes)
Outputs:
  results/0_prep/02.yoo25_l23_abc_scores.tsv    (step 1; cell, cond, score_a, score_b, score_c)
  results/0_prep/02.yoo25_l23_abc_cells.tsv     (step 2; cell, cond, archetype: A/B/C, '' = none)
  results/0_prep/02.yoo25_l23_abc_genes.tsv     (step 3, = Table S4; cond, archetype, gene, l2fc_max/min/avg, q_min)
  results/0_prep/02.yoo25_l23_abc_gene_stats.parquet  (step 3, every test gene: l2fc_vs_X and q_vs_X
                                                for each other archetype X, l2fc_max/min/avg, q_min, is_gene)
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_H5AD  = os.path.join(REPO_ROOT, 'data', 'raw', 'superdupermegaRNA_hasraw.h5ad')
INPUT_GENES = os.path.join(REPO_ROOT, 'data', 'raw', 'Cheng22_Xie25_L23_ABC_genes_n286.csv')
OUT_SCORES  = os.path.join(REPO_ROOT, 'results', '0_prep', '02.yoo25_l23_abc_scores.tsv')
OUT_CELLS   = os.path.join(REPO_ROOT, 'results', '0_prep', '02.yoo25_l23_abc_cells.tsv')
OUT_GENES   = os.path.join(REPO_ROOT, 'results', '0_prep', '02.yoo25_l23_abc_genes.tsv')
OUT_STATS   = os.path.join(REPO_ROOT, 'results', '0_prep', '02.yoo25_l23_abc_gene_stats.parquet')

STUDY    = '2023 Multiome'
SUBCLASS = 'L2/3'
CONDS = ['P6', 'P8', 'P10', 'P12', 'P14', 'P17', 'P21', 'P12DR', 'P14DR', 'P17DR', 'P21DR']
ARCHS = ['A', 'B', 'C']
# step 1
EPS = 1e-10                    # min-max denominator
# step 2
N_CELLS = 100                  # archetype cells per group and condition
# step 3
MIN_CELLS     = 10             # keep test genes detected in more than this many cells
OFFSET        = 1              # log2(mean CP10k + OFFSET)
ALPHA         = 0.05
L2FC_TH       = 1
L2FC_TH_B_MIN = np.log2(1.2)   # B only

os.makedirs(os.path.dirname(OUT_GENES), exist_ok=True)

# ---- inputs: counts and gene sets ----
src = ad.read_h5ad(INPUT_H5AD, backed='r')
if src.raw is None:
    raise ValueError(f'{INPUT_H5AD} has no .raw; raw counts are expected there')
obs = src.obs[['Age', 'Study', 'Subclass']].astype(str)
idx = np.flatnonzero(((obs['Study'] == STUDY) & (obs['Subclass'] == SUBCLASS)).to_numpy())
cell_ids = obs.index[idx]
cond = obs['Age'].iloc[idx].to_numpy()
if sorted(set(cond)) != sorted(CONDS):
    raise ValueError(f'expected conditions {CONDS}, found {sorted(set(cond))}')
X = sp.csr_matrix(src.raw.X[idx]).astype(np.float64)
genes = src.raw.var_names
if not np.allclose(X.data, np.round(X.data)):
    raise ValueError(f'{INPUT_H5AD} .raw.X is not integer counts')
print(f'{len(idx)} cells x {len(genes)} genes')

anno = pd.read_csv(INPUT_GENES)
score_sets = {}
for k in ARCHS:
    g = anno.loc[anno['P17on'] == k, 'gene'].astype(str)
    score_sets[k] = genes.get_indexer(np.intersect1d(g, genes))
    print(f'score genes {k}: {len(g)}, {len(score_sets[k])} in Yoo25')

# ---- step 1: scores ----
lib_all = np.asarray(X.sum(1)).ravel()
scores = pd.DataFrame(index=pd.Index(cell_ids, name='cell'))
scores['cond'] = cond
for k, gi in score_sets.items():
    L = np.log2(1 + X[:, gi].toarray() / lib_all[:, None] * 1e4)
    s = np.empty(len(idx))
    for c in CONDS:
        m = cond == c
        lo, hi = L[m].min(0), L[m].max(0)
        s[m] = ((L[m] - lo) / (hi - lo + EPS)).mean(1)
    scores[f'score_{k.lower()}'] = s
scores.to_csv(OUT_SCORES, sep='\t')

# ---- step 2: archetype cells ----
axis_ca = (scores['score_c'] - scores['score_a']).to_numpy()
score_b = scores['score_b'].to_numpy()
arch = np.full(len(idx), '', dtype=object)
for c in CONDS:
    m = np.flatnonzero(cond == c)
    r_ca, r_b = pd.Series(axis_ca[m]).rank().to_numpy(), pd.Series(score_b[m]).rank().to_numpy()
    top = {'A': r_ca <= N_CELLS, 'B': r_b > len(m) - N_CELLS, 'C': r_ca > len(m) - N_CELLS}
    multi = np.sum(list(top.values()), axis=0) > 1
    for k, v in top.items():
        arch[m[v & ~multi]] = k
cells = scores[['cond']].assign(archetype=arch)
cells.to_csv(OUT_CELLS, sep='\t')
print(cells[cells['archetype'] != ''].groupby(['archetype', 'cond']).size().unstack()[CONDS].to_string())

# ---- step 3: archetype genes ----
keep = ~genes.str.contains(r'^mt-|^Xist$') & (np.asarray((X > 0).sum(0)).ravel() > MIN_CELLS)
Xt, test_genes = X[:, keep], genes[keep].to_numpy()
cp10k = sp.diags(1 / np.asarray(Xt.sum(1)).ravel()) @ Xt * 1e4
print(f'test genes: {len(test_genes)}')

stats_rows = []
for c in CONDS:
    rows = np.flatnonzero((cond == c) & (arch != ''))   # only the archetype cells are needed
    xn, grp = cp10k[rows].toarray(), arch[rows]
    cp = {k: xn[grp == k] for k in ARCHS}
    mu = {k: np.log2(v.mean(0) + OFFSET) for k, v in cp.items()}
    lg = {k: np.log2(v + 1) for k, v in cp.items()}
    q = {}
    for i, a in enumerate(ARCHS):
        for b in ARCHS[i + 1:]:
            q[a, b] = q[b, a] = multipletests(np.nan_to_num(stats.ttest_ind(lg[a], lg[b]).pvalue, nan=1),
                                              method='fdr_bh')[1]
    for k in ARCHS:
        others = [j for j in ARCHS if j != k]
        fc = np.stack([mu[k] - mu[j] for j in others])
        fmax, fmin = fc.max(0), fc.min(0)
        qmin = np.minimum(*[q[k, j] for j in others])
        sig = (fmin > 0) & (fmax > L2FC_TH) & (qmin < ALPHA)
        if k == 'B':
            sig &= fmin > L2FC_TH_B_MIN
        stats_rows.append(pd.DataFrame({
            'cond': c, 'archetype': k, 'gene': test_genes,
            **{f'l2fc_vs_{j}': mu[k] - mu[j] if j != k else np.nan for j in ARCHS},
            **{f'q_vs_{j}': q[k, j] if j != k else np.nan for j in ARCHS},
            'l2fc_max': fmax, 'l2fc_min': fmin, 'l2fc_avg': (fmax + fmin) / 2, 'q_min': qmin, 'is_gene': sig}))
gene_stats = pd.concat(stats_rows, ignore_index=True)
gene_stats.to_parquet(OUT_STATS, index=False)
res = gene_stats.loc[gene_stats['is_gene'], ['cond', 'archetype', 'gene', 'l2fc_max', 'l2fc_min', 'l2fc_avg', 'q_min']]
res.to_csv(OUT_GENES, sep='\t', index=False)
print(res.groupby(['archetype', 'cond']).size().unstack()[CONDS].to_string())
print(f'wrote {OUT_SCORES}\nwrote {OUT_CELLS}\nwrote {OUT_GENES}\nwrote {OUT_STATS}')
