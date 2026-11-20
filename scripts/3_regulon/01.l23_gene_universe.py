"""Gene universe of the regulon enrichment: genes expressed in adult L2/3 IT neurons (Cheng22 P28 NR + Yoo25 P21 NR).

Cells: the L2/3 cells of the 0_prep/01 extracts, Cheng22 (all, P28 NR) first, then Yoo25 Age == 'P21' (NR),
each in file order.
Genes: the shared gene set (inner join, Cheng22 order, mito removed); expression log2(CP10k+1) from the
raw counts with each dataset's validated library-size column (cheng22 'n_counts', yoo25
'total_counts'). Universe = genes with nonzero total and nonzero variance over all cells; their
mean / var / detection rate of log2(CP10k+1) are what 02's expression strata are built from.

Reads:
  data/input/01.cheng22_p28nr_l23.h5ad
  data/input/01.yoo25_l23.h5ad
Outputs:
  results/3_regulon/01.l23_gene_universe_stats.tsv   (gene, mean, var, detect)
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_H5AD = {'cheng22': os.path.join(REPO_ROOT, 'data', 'input', '01.cheng22_p28nr_l23.h5ad'),
              'yoo25': os.path.join(REPO_ROOT, 'data', 'input', '01.yoo25_l23.h5ad')}
OUT_STATS  = os.path.join(REPO_ROOT, 'results', '3_regulon', '01.l23_gene_universe_stats.tsv')

YOO25_AGE    = 'P21'   # NR
DEPTH_COL    = {'cheng22': 'n_counts', 'yoo25': 'total_counts'}

os.makedirs(os.path.dirname(OUT_STATS), exist_ok=True)

subs = {}
for tag, path in INPUT_H5AD.items():
    a = ad.read_h5ad(path)
    subs[tag] = a[a.obs['Age'] == YOO25_AGE].copy() if tag == 'yoo25' else a
    print(f'  {tag:8s}: {subs[tag].n_obs} L2/3 cells, {subs[tag].n_vars} genes')

# shared gene set (inner join), ordered by cheng22, mito removed
common_set = set(subs['cheng22'].var_names) & set(subs['yoo25'].var_names)
common_genes = [g for g in subs['cheng22'].var_names if g in common_set]
common_genes = [g for g in common_genes if not g.lower().startswith('mt-')]

# log2(CP10k+1) per dataset (validated depth + raw counts), then merged
frames = []
for tag, a in subs.items():
    depth = a.obs[DEPTH_COL[tag]].values.astype(np.float64)
    if not (np.all(np.isfinite(depth)) and np.all(depth > 0)):
        raise ValueError(f"{tag}: invalid depth column '{DEPTH_COL[tag]}' (NaN or <=0)")
    xc = a[:, common_genes].X
    xc = xc.toarray() if sp.issparse(xc) else np.asarray(xc, dtype=np.float64)
    frames.append(pd.DataFrame(np.log2(xc / depth[:, None] * 1e4 + 1).astype(np.float32), columns=common_genes))
X = pd.concat(frames).values   # (the float32 sums below depend on this memory layout)

keep = (X.sum(axis=0) > 0) & (X.var(axis=0) > 0)
Xk = X[:, keep]
stats = pd.DataFrame({'mean': Xk.mean(axis=0), 'var': Xk.var(axis=0), 'detect': (Xk > 0).mean(axis=0)},
                     index=pd.Index(np.asarray(common_genes)[keep], name='gene'))
stats.to_csv(OUT_STATS, sep='\t')
print(f'gene universe: {keep.sum()} expressed nonzero-variance genes over {X.shape[0]} cells '
      f'(dropped {(~keep).sum()} of {keep.size})')
print(f'wrote {OUT_STATS}')
