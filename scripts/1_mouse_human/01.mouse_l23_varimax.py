"""Mouse L2/3 IT (Cheng22, P28 NR): PCA + varimax rotation + variance partitioning, top-2000 HVGs.

Normalization: cheng22 .X is log1p(CP10k), so the integer raw counts in .raw are re-normalized with
log2(CP10k + 1). HVGs = top N_HVG by raw variance of that matrix; StandardScaler -> PCA(N_PCS) ->
Kaiser varimax; VX components ordered by descending score variance. The post-varimax variance
partition (cell_type vs sample vs library_size) is what the Gate-A VX selection (VX dominated by
cell type) of scripts 03 and 05 is read from.

Reads:
  data/raw/superdupermegaRNA_cheng22_IT_P28NR.h5ad
Outputs:
  results/1_mouse_human/01.mouse_l23_varimax_coords.tsv
  results/1_mouse_human/01.mouse_l23_varimax_loadings.tsv
  results/1_mouse_human/01.mouse_l23_pc_variance_partition.tsv
  results/1_mouse_human/01.mouse_l23_vx_variance_partition.tsv
"""

import os
import numpy as np
import pandas as pd
import anndata as ad

from l23.crossspecies import pca_varimax, variance_partition, dummies, library_size_factor, mouse_log2_cp10k

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_MOUSE     = os.path.join(REPO_ROOT, 'data', 'raw', 'superdupermegaRNA_cheng22_IT_P28NR.h5ad')
OUT_DIR         = os.path.join(REPO_ROOT, 'results', '1_mouse_human')
OUT_VX_COORDS   = os.path.join(OUT_DIR, '01.mouse_l23_varimax_coords.tsv')
OUT_VX_LOADINGS = os.path.join(OUT_DIR, '01.mouse_l23_varimax_loadings.tsv')
OUT_PC_VAR_TSV  = os.path.join(OUT_DIR, '01.mouse_l23_pc_variance_partition.tsv')
OUT_VX_VAR_TSV  = os.path.join(OUT_DIR, '01.mouse_l23_vx_variance_partition.tsv')

MOUSE_SUBCLASS = 'L2/3'
CLUSTER_COL    = 'Type'
SAMPLE_COL     = 'Sample'
LIBSIZE_COL    = 'n_counts'
N_HVG          = 2000
N_PCS          = 10

os.makedirs(OUT_DIR, exist_ok=True)

print('Loading mouse cheng22 data...')
m_adata = ad.read_h5ad(INPUT_MOUSE)
m_adata = m_adata[m_adata.obs['Subclass'] == MOUSE_SUBCLASS]
print(f'  {m_adata.n_obs} cells, {m_adata.n_vars} genes')
types   = m_adata.obs[CLUSTER_COL].values
samples = m_adata.obs[SAMPLE_COL].values
factors = {'cell_type': dummies(types, 'cell_type'), 'sample': dummies(samples, 'sample'),
           'library_size': library_size_factor(m_adata.obs[LIBSIZE_COL].values)}

X_norm = mouse_log2_cp10k(m_adata)
print(f'Top {N_HVG} HVGs -> PCA -> varimax...')
hvg_names, scores, vx_scores, vx_loadings = pca_varimax(X_norm, m_adata.var_names.values, N_HVG, N_PCS)

pc_cols = [f'PC{i+1}' for i in range(N_PCS)]
vx_cols = [f'VX{i+1}' for i in range(N_PCS)]
variance_partition(scores, pc_cols, factors, 'pre-varimax PCs').to_csv(OUT_PC_VAR_TSV, sep='\t')
variance_partition(vx_scores, vx_cols, factors, 'post-varimax VX').to_csv(OUT_VX_VAR_TSV, sep='\t')

vx_df = pd.DataFrame(vx_scores, index=m_adata.obs_names, columns=vx_cols)
vx_df[CLUSTER_COL] = types
vx_df[SAMPLE_COL]  = samples
vx_df.to_csv(OUT_VX_COORDS, sep='\t')
pd.DataFrame(vx_loadings, index=hvg_names, columns=vx_cols).to_csv(OUT_VX_LOADINGS, sep='\t')
for p in (OUT_PC_VAR_TSV, OUT_VX_VAR_TSV, OUT_VX_COORDS, OUT_VX_LOADINGS):
    print(f'Saved {p}')
