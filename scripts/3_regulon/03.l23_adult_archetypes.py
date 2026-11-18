"""Adult L2/3 archetypes (n=3) on the two-dataset Harmony embedding of 02, and their marker genes.

PCHA with NOC=3 on the cheng22 P28 NR + yoo25 P21 NR L2/3 embedding, then the genes associated
with each archetype. Per-cell expression is reconstructed across both datasets and aligned to
the Harmony coords' dataset:barcode cell order; the marker test runs over ALL expressed genes
with nonzero variance (the shared two-dataset gene set, mito removed), not just the HVGs.

Archetype labels: archetype_1, _2, _3 = C', B', A' (laminar-depth order A' superficial).

Procedure:
  1) load Harmony coords; use H1..H5 as the PCHA feature space
  2) fit PCHA (NDIM=5, NOC=3); back-project archetypes to H-space
  3) take the N_TOP_CELLS cells closest to each archetype
  4) reconstruct CP10k->log2 expression for the full shared gene set across both datasets,
     filter to expressed (nonzero total) + nonzero-variance genes
  5) conservative one-vs-each Wilcoxon per archetype -> markers (worst-case p / log2FC)

Reads:
  results/3_regulon/02.l23_harmony_coords.tsv   (H1..H10 + Type + Sample + Dataset)
  data/raw/superdupermegaRNA_cheng22_IT_P28NR.h5ad
  data/raw/superdupermegaRNA_yoo25_IT_P21.h5ad
Outputs:
  results/3_regulon/03.l23_adult_archetype_markers.tsv
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import scipy.stats
from statsmodels.stats.multitest import multipletests

from l23.archetype import SCA

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# --- file paths ---
IN_COORDS      = os.path.join(REPO_ROOT, 'results', '3_regulon', '02.l23_harmony_coords.tsv')
INPUT_CHENG22  = os.path.join(REPO_ROOT, 'data', 'raw', 'superdupermegaRNA_cheng22_IT_P28NR.h5ad')
INPUT_YOO25    = os.path.join(REPO_ROOT, 'data', 'raw', 'superdupermegaRNA_yoo25_IT_P21.h5ad')
OUT_MARKERS    = os.path.join(REPO_ROOT, 'results', '3_regulon', '03.l23_adult_archetype_markers.tsv')

# --- per-dataset config (as 02) ---
DATASETS = [
    dict(tag='cheng22', path=INPUT_CHENG22, subclass_col='Subclass',
         subclass_val='L2/3', sample_col='Sample', depth_col='n_counts'),
    dict(tag='yoo25', path=INPUT_YOO25, subclass_col='Subclass',
         subclass_val='L2/3', sample_col='Sample', depth_col='total_counts'),
]

# --- parameters ---
CLUSTER_COL     = 'Type'
SAMPLE_COL      = 'Sample'
DATASET_COL     = 'Dataset'
H_COLS          = ['H1', 'H2', 'H3', 'H4', 'H5']   # top-5 Harmony PCs
NOC             = 3
NDIM            = len(H_COLS)   # 5
DROP_PCS        = []            # use all 5 Harmony PCs directly (no PC dropped)
N_TOP_CELLS     = 300
FRAC_IN_THRESH  = 0.25
FDR_THRESH      = 0.001
LOG2FC_THRESH   = np.log2(1.5)

os.makedirs(os.path.dirname(OUT_MARKERS), exist_ok=True)

# ===================== load Harmony coords =====================

print(f'Loading Harmony coords from {IN_COORDS}...')
coords  = pd.read_csv(IN_COORDS, sep='\t', index_col=0)
xn      = coords[H_COLS].values                       # (n_cells, 5) PCHA feature space
types   = coords[CLUSTER_COL].astype(str).values
samples = coords[SAMPLE_COL].astype(str).values
cell_index = coords.index.values
n_cells = xn.shape[0]
print(f'  {n_cells} cells, {len(np.unique(samples))} dataset:sample groups')

# ===================== fit PCHA =====================

print(f'Fitting PCHA: NOC={NOC}, NDIM={NDIM}...')
sca = SCA(xn, types)
sca.setup_feature_matrix(method='data')
sca.proj_and_pcha(NDIM, NOC, drop_pcs=DROP_PCS)

# back-project archetype coords from PCHA space to the H feature space
aa_feat = sca.aa.T @ sca.pca_.components_[:NDIM] + sca.pca_.mean_   # (NOC, NDIM)

# --- cell-to-archetype distances in H feature space ---
dists     = np.stack([np.linalg.norm(xn - aa_feat[k], axis=1) for k in range(NOC)], axis=1)
top_cells = [np.argsort(dists[:, k])[:N_TOP_CELLS] for k in range(NOC)]
print(f'Top {N_TOP_CELLS} cells selected per archetype.')

# ===================== reconstruct expression across both datasets =====================

# --- 1. load each dataset, filter to L2/3 ---
print('Loading datasets and filtering to L2/3...')
adatas = {}
for d in DATASETS:
    a = ad.read_h5ad(d['path'])
    a = a[a.obs[d['subclass_col']] == d['subclass_val']].copy()
    adatas[d['tag']] = a
    print(f"  {d['tag']:8s}: {a.n_obs} cells ({d['subclass_val']}), {a.n_vars} genes")

# --- 2. shared gene set (inner join), ordered by cheng22, mito removed ---
common_set = set(adatas['cheng22'].var_names)
for tag in [d['tag'] for d in DATASETS if d['tag'] != 'cheng22']:
    common_set &= set(adatas[tag].var_names)
common_genes = [g for g in adatas['cheng22'].var_names if g in common_set]
n_before = len(common_genes)
common_genes = [g for g in common_genes if not g.lower().startswith('mt-')]
print(f"  Shared genes: {n_before} (intersection); "
      f"{len(common_genes)} after removing {n_before - len(common_genes)} mito genes")

# --- 3. CP10k->log2 expression per dataset (validated depth + raw counts), then merge ---
expr_frames = []
for d in DATASETS:
    a = adatas[d['tag']]
    depth = a.obs[d['depth_col']].values.astype(np.float64)
    assert np.all(np.isfinite(depth)) and np.all(depth > 0), \
        f"{d['tag']}: invalid depth column '{d['depth_col']}' (NaN or <=0)"
    xc = a.raw[:, common_genes].X
    xc = xc.toarray() if sp.issparse(xc) else np.asarray(xc, dtype=np.float64)
    xc = np.log2(xc / depth[:, None] * 1e4 + 1).astype(np.float32)   # CP10k -> log2(1+x)
    idx = np.array([f"{d['tag']}:{n}" for n in a.obs_names])
    expr_frames.append(pd.DataFrame(xc, index=idx, columns=common_genes))
    print(f"  {d['tag']:8s}: depth median {np.median(depth):.0f}  ->  {xc.shape[0]} cells")

expr_df = pd.concat(expr_frames)
# align to coords cell order (fail-fast if any coord cell is missing from the merged expression)
missing = np.setdiff1d(cell_index, expr_df.index.values)
assert missing.size == 0, f'{missing.size} coord cells missing from reconstructed expression, e.g. {missing[:5]}'
expr_df = expr_df.reindex(cell_index)

# --- 4. filter to expressed (nonzero total) + nonzero-variance genes ---
X_norm    = expr_df.values
expressed = X_norm.sum(axis=0) > 0
nonzero_var = X_norm.var(axis=0) > 0
keep = expressed & nonzero_var
X_norm    = X_norm[:, keep]
gene_names = expr_df.columns.values[keep]
n_genes = X_norm.shape[1]   # gene universe size
print(f'Gene universe for marker test: {n_genes} expressed nonzero-variance genes '
      f'(dropped {(~keep).sum()} of {keep.size})')
print(f'X_norm shape: {X_norm.shape}')

# ===================== conservative one-vs-each Wilcoxon per archetype =====================

print('Running conservative one-vs-each Wilcoxon per archetype (all expressed genes)...')
all_markers = []

for k in range(NOC):
    archetype_label = f'archetype_{k + 1}'
    print(f'  {archetype_label}...')
    in_idx  = top_cells[k]
    X_in    = X_norm[in_idx]
    frac_in = (X_in > 0).mean(axis=0)

    log2fc_mat, pval_mat, frac_out_list = [], [], []
    for j in range(NOC):
        if j == k:
            continue
        X_out    = X_norm[top_cells[j]]
        log2fc_j = (X_in.mean(axis=0) - X_out.mean(axis=0)) / np.log(2)
        pvals_j  = np.array([scipy.stats.ranksums(X_in[:, g], X_out[:, g])[1]
                              for g in range(n_genes)])
        log2fc_mat.append(log2fc_j)
        pval_mat.append(pvals_j)
        frac_out_list.append((X_out > 0).mean(axis=0))

    # worst-case across pairwise comparisons
    log2fc   = np.stack(log2fc_mat).min(axis=0)
    pvals    = np.stack(pval_mat).max(axis=0)
    frac_out = np.stack(frac_out_list).mean(axis=0)

    _, fdr, _, _ = multipletests(pvals, method='fdr_bh')

    df = pd.DataFrame({
        'gene':      gene_names,
        'archetype': archetype_label,
        'log2FC':    log2fc,
        'pval':      pvals,
        'fdr':       fdr,
        'frac_in':   frac_in,
        'frac_out':  frac_out,
    })
    df = df[(df['frac_in'] >= FRAC_IN_THRESH) & (df['fdr'] < FDR_THRESH) & (df['log2FC'] > LOG2FC_THRESH)]
    df = df.sort_values('log2FC', ascending=False)
    all_markers.append(df)
    print(f'    {len(df)} markers after filtering')

markers_df = pd.concat(all_markers, ignore_index=True)
markers_df.to_csv(OUT_MARKERS, sep='\t', index=False)
print(f'Saved {OUT_MARKERS}  ({len(markers_df)} total markers)')

print('Done.')
