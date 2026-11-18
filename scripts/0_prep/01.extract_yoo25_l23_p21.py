"""Extract Yoo25 P21 L2/3 neurons, normally reared (NR) and dark reared (DR), with raw counts.

Cells: Subclass == 'L2/3' (not Type: type labels differ by condition) and
  NR: Age == 'P21' and Study == '2023 Multiome' (the other P21 samples are Cheng22)
  DR: Age == 'P21DR'
Raw integer counts (source .raw.X) are written to .X. Two obs columns are added:
  condition  NR / DR
  replicate  the pseudobulk unit (= Sample)

Reads:
  data/raw/superdupermegaRNA_hasraw.h5ad
Outputs:
  results/0_prep/01.yoo25_l23_p21nrdr.h5ad
  results/0_prep/01.cell_counts.tsv
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_SOURCE = os.path.join(REPO_ROOT, 'data', 'raw', 'superdupermegaRNA_hasraw.h5ad')
OUT_H5AD     = os.path.join(REPO_ROOT, 'results', '0_prep', '01.yoo25_l23_p21nrdr.h5ad')
OUT_COUNTS   = os.path.join(REPO_ROOT, 'results', '0_prep', '01.cell_counts.tsv')

SUBCLASS = 'L2/3'
STUDY    = '2023 Multiome'
GROUPS = {   # condition -> Age and the expected samples
    'NR': {'Age': 'P21', 'samples': ['P21a', 'P21b']},
    'DR': {'Age': 'P21DR', 'samples': ['P21DRa', 'P21DRb']},
}

os.makedirs(os.path.dirname(OUT_H5AD), exist_ok=True)

print(f'Reading {INPUT_SOURCE} (backed)', flush=True)
src = ad.read_h5ad(INPUT_SOURCE, backed='r')
if src.raw is None:
    raise ValueError(f'{INPUT_SOURCE} has no .raw; raw counts are expected there')
if not src.raw.var_names.equals(src.var_names):
    raise ValueError('.raw.var_names differ from .var_names')
obs = src.obs[['Age', 'Study', 'Sample', 'Subclass']].astype(str)

masks, cond = [], np.full(len(obs), '', dtype=object)
for c, g in GROUPS.items():
    m = ((obs['Age'] == g['Age']) & (obs['Study'] == STUDY) & (obs['Subclass'] == SUBCLASS)).to_numpy()
    found = sorted(obs.loc[m, 'Sample'].unique())
    if found != sorted(g['samples']):
        raise ValueError(f'{c}: expected samples {g["samples"]}, found {found}')
    cond[m] = c
    masks.append(m)
idx = np.flatnonzero(np.logical_or.reduce(masks))

print(f'loading {len(idx)} L2/3 cells', flush=True)
X = sp.csr_matrix(src.raw.X[idx])
if X.data.size and not np.allclose(X.data, np.round(X.data)):
    raise ValueError('.raw.X is not integer counts')
sub_obs = src.obs.iloc[idx].copy()
sub_obs['condition'] = pd.Categorical(cond[idx], categories=list(GROUPS))
sub_obs['replicate'] = sub_obs['Sample'].astype(str)

out = ad.AnnData(X=X.astype(np.float32), obs=sub_obs, var=src.var.copy())
out.write_h5ad(OUT_H5AD)
print(f'wrote {OUT_H5AD} {out.shape}', flush=True)

counts = sub_obs.groupby(['condition', 'replicate'], observed=True).size().rename('n_cells').reset_index()
counts.to_csv(OUT_COUNTS, sep='\t', index=False)
print(counts.to_string(index=False))
print(f'wrote {OUT_COUNTS}')
