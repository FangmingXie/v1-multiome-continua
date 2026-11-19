"""Extract the L2/3 neurons with raw counts: Yoo25 of every Age (NR and DR), and Cheng22 P28 NR.

Every later step that needs Yoo25 or Cheng22 L2/3 counts reads these extracts instead of the combined object.
Both: Subclass == 'L2/3' (not Type: type labels differ by condition), in the source row order, raw integer
counts (source .raw.X) written to .X, every source obs column kept.
  Yoo25    Study == '2023 Multiome', all 11 Ages of l23.constants.AGES (NR P6-P21, DR P12-P21).
           obs columns added: day, condition (NR / DR) and lib (library size = sum over all genes).
  Cheng22  Study missing (Cheng22 has no Study label in the combined object) and Age == 'P28' (NR; DR and
           dark-light P28 cells have their own Ages, P28_dr and P28_dl).

Reads:
  data/raw/superdupermegaRNA_hasraw.h5ad
Outputs:
  results/0_prep/01.yoo25_l23.h5ad
  results/0_prep/01.cheng22_p28nr_l23.h5ad
  results/0_prep/01.cell_counts.tsv     (Yoo25 cells per Age and Sample)
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

from l23.constants import YOO25_STUDY, CHENG22_AGE, L23_SUBCLASS, AGES

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_SOURCE = os.path.join(REPO_ROOT, 'data', 'raw', 'superdupermegaRNA_hasraw.h5ad')
OUT_H5AD     = os.path.join(REPO_ROOT, 'results', '0_prep', '01.yoo25_l23.h5ad')
OUT_CHENG22  = os.path.join(REPO_ROOT, 'results', '0_prep', '01.cheng22_p28nr_l23.h5ad')
OUT_COUNTS   = os.path.join(REPO_ROOT, 'results', '0_prep', '01.cell_counts.tsv')

os.makedirs(os.path.dirname(OUT_H5AD), exist_ok=True)


def raw_counts(src, idx):
    """CSR raw integer counts (source .raw.X) of rows idx."""
    X = sp.csr_matrix(src.raw.X[idx])
    if X.data.size and not np.allclose(X.data, np.round(X.data)):
        raise ValueError('.raw.X is not integer counts')
    return X


print(f'Reading {INPUT_SOURCE} (backed)', flush=True)
src = ad.read_h5ad(INPUT_SOURCE, backed='r')
if src.raw is None:
    raise ValueError(f'{INPUT_SOURCE} has no .raw; raw counts are expected there')
if not src.raw.var_names.equals(src.var_names):
    raise ValueError('.raw.var_names differ from .var_names')
obs = src.obs[['Age', 'Study', 'Sample', 'Subclass']].astype(str)

# --- Yoo25, every Age ---
sel = ((obs['Study'] == YOO25_STUDY) & (obs['Subclass'] == L23_SUBCLASS)).to_numpy()
ages = sorted(obs.loc[sel, 'Age'].unique())
if ages != sorted(AGES):
    raise ValueError(f'expected {YOO25_STUDY} {L23_SUBCLASS} ages {sorted(AGES)}, found {ages}')
if (obs[sel].groupby('Sample')['Age'].nunique() != 1).any():
    raise ValueError('a sample spans more than one Age')
idx = np.flatnonzero(sel)

print(f'loading {len(idx)} Yoo25 L2/3 cells', flush=True)
X = raw_counts(src, idx)
sub_obs = src.obs.iloc[idx].copy()
age = obs['Age'].iloc[idx]
sub_obs['day'] = age.map(lambda a: AGES[a][0]).to_numpy()
sub_obs['condition'] = age.map(lambda a: AGES[a][1]).to_numpy()
sub_obs['lib'] = np.asarray(X.sum(axis=1), dtype=np.float64).ravel()
if (sub_obs['lib'] <= 0).any():
    raise ValueError('cells with zero library size')

out = ad.AnnData(X=X.astype(np.float32), obs=sub_obs, var=src.var.copy())
out.write_h5ad(OUT_H5AD)
print(f'wrote {OUT_H5AD} {out.shape}', flush=True)

# --- Cheng22, P28 NR ---
idx = np.flatnonzero((src.obs['Study'].isna() & (obs['Age'] == CHENG22_AGE)
                      & (obs['Subclass'] == L23_SUBCLASS)).to_numpy())
if len(idx) == 0:
    raise ValueError(f'no Cheng22 (Study missing) {CHENG22_AGE} {L23_SUBCLASS} cells')
print(f'loading {len(idx)} Cheng22 {CHENG22_AGE} L2/3 cells', flush=True)
c_obs = src.obs.iloc[idx].copy()
for col in c_obs.select_dtypes('category'):   # read whole (no AnnData subsetting), so drop the unused levels here
    c_obs[col] = c_obs[col].cat.remove_unused_categories()
out = ad.AnnData(X=raw_counts(src, idx).astype(np.float32), obs=c_obs, var=src.var.copy())
out.write_h5ad(OUT_CHENG22)
print(f'wrote {OUT_CHENG22} {out.shape}', flush=True)

counts = (sub_obs.astype({'Age': str, 'Sample': str}).groupby(['Age', 'Sample']).size().rename('n_cells')
          .reset_index())
counts['Age'] = pd.Categorical(counts['Age'], categories=list(AGES))
counts = counts.sort_values(['Age', 'Sample'])
counts.to_csv(OUT_COUNTS, sep='\t', index=False)
print(counts.to_string(index=False))
print(f'wrote {OUT_COUNTS}')
